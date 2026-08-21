# -*- coding: utf-8 -*-
"""
kcat_optimization.py — 适配 calculate_ETM 的酶动力学参数优化模块

功能说明
--------
该模块用于迭代优化 ETGEM 模型中反应的 kcat 和 kcat_MW 参数，以提升
酶-热约束模型（ETM）的求解性能。核心思路来源于 ECMpy 的酶用量校正方法
（change_enz_model_by_enz_usage），但直接适配到 ETGEM 框架的 Model 对象，
无需中间 JSON 文件转换。

核心流程
--------
1. 运行 calculate_ETM 获取初始解
2. 计算每个反应的酶用量 E_j = v_j * mw / (kcat * 3600000)
3. 识别酶用量占比最高的瓶颈反应
4. 根据 EC 编号查询 kcat 数据库（EC_kcat_max.json），获取更高的 kcat 值
5. 更新模型对象的 reaction.kcat 和 reaction.kcat_MW
6. 重新运行 calculate_ETM 验证性能提升
7. 迭代直至收敛或达到最大迭代次数

与 calculate_ETM 的兼容性
-------------------------
- 输入 model 为 ETGEM 的 Model 对象（core.model.Model）
- 调整后可直接传入 calculate_ETM(model, "maximize", B_value) 使用
"""

import json
import os
import numpy as np
import pandas as pd
from copy import deepcopy
from typing import Dict, List, Tuple, Optional, Any


# ─── 自包含 JSON 模型序列化（不依赖 io/json.py，避免影响其他模块）───────────────
# 原始 io/json.py 中的 _serialize_* 函数遍历 dict 时误用 keys()，
# 导致保存失败。此处内嵌修复版，仅供本模块保存模型使用。

def _serialize_genes(genes: Dict[str, Any]) -> List[Dict[str, Any]]:
    """基因对象转化为字典（修复版：遍历 .values()"""
    genelist = []
    for gene in genes.values():
        if "uniprot" in gene.annotation.keys():
            genelist.append({"id": gene.id, "name": gene.name,
                             "uniprot": gene.annotation["uniprot"],
                             "annotation": gene.annotation})
        else:
            genelist.append({"id": gene.id, "name": gene.name,
                             "uniprot": "", "annotation": gene.annotation})
    return genelist


def _serialize_reactions(reactions: Dict[str, Any]) -> List[Dict[str, Any]]:
    """反应对象转化为字典（修复版：完整字段 + 兼容 metabolites 键类型）。"""
    serialized_reactions = []
    for reaction in reactions.values():
        # reaction.metabolites 可能是 {id_str: coef} 或 {Metabolite: coef}
        metabolites = {}
        for met, value in reaction.metabolites.items():
            met_id = met.id if hasattr(met, 'id') else met
            metabolites[met_id] = value

        # kcat_MW：若对象已有此属性则用对象值，否则按公式计算
        if hasattr(reaction, 'kcat_MW') and reaction.kcat_MW not in (None, ''):
            kcat_MW = reaction.kcat_MW
        else:
            try:
                kcat_MW = reaction.kcat / reaction.mw * 3600000
            except (TypeError, ZeroDivisionError):
                kcat_MW = 0

        # 浓度上下界（默认 -14.50865774 / -3.912023005，与 io/json.py 一致）
        concentration_lb = getattr(reaction, 'concentration_lb', None)
        if concentration_lb in (None, ''):
            concentration_lb = -3.912023005
        concentration_ub = getattr(reaction, 'concentration_ub', None)
        if concentration_ub in (None, ''):
            concentration_ub = -14.50865774

        # 目标系数（默认 None；若所有反应均为空，则在含"BIOMASS"的反应中设为 1.0）
        objective_coefficient = getattr(reaction, 'objective_coefficient', None)
        if objective_coefficient in (None, ''):
            objective_coefficient = 0  # 由 write_json_model 统一补全

        notes = getattr(reaction, 'notes', {}) or {}

        serialized_reaction = {
            "id": reaction.id,
            "name": reaction.name,
            "metabolites": metabolites,
            "lower_bound": reaction.lower_bound,
            "upper_bound": reaction.upper_bound,
            "gene_reaction_rule": reaction.gene_reaction_rule,
            "notes": notes,
            "annotation": reaction.annotation,
            "kcat": reaction.kcat,
            "kcat_MW": kcat_MW,
            "mw": reaction.mw,
            "g0": reaction.g0,
            "concentration_lb": concentration_lb,
            "concentration_ub": concentration_ub,
            "objective_coefficient": objective_coefficient,
        }
        serialized_reactions.append(serialized_reaction)
    return serialized_reactions


def _serialize_metabolites(metabolites: Dict[str, Any]) -> List[Dict[str, Any]]:
    """代谢物对象转化为字典（修复版：遍历 .values()）"""
    metabolitelist = []
    for metabolite in metabolites.values():
        serialize_metabolite = {"id": metabolite.id,
                                "name": metabolite.name,
                                "compartment": metabolite.compartment,
                                "charge": metabolite.charge,
                                "formula": metabolite.formula,
                                "concentration_ub": metabolite.concentration_ub,
                                "concentration_lb": metabolite.concentration_lb,
                                "annotation": metabolite.annotation}
        metabolitelist.append(serialize_metabolite)
    return metabolitelist


def _serialize_environment(environment: Any) -> Dict[str, Any]:
    """环境对象转化为字典"""
    return {
        "p_h": environment.p_h,
        "p_mg": environment.p_mg,
        "ionic_strength": environment.ionic_strength,
        "temperature": environment.temperature
    }


def write_json_model(model: Any, filename: str,
                     template_path: str = None) -> None:
    """将模型完整写入 JSON 文件（本地实现，不依赖 io/json.py）。

    完整保留顶层 id/compartments/version/enzyme_constraint 及反应
    kcat_MW/notes/concentration_lb|ub/objective_coefficient 等所有字段，
    格式严格对齐 iZM564_SYN_etcmodel.json。

    Parameters
    ----------
    model : etgems_web.Model
        要保存的模型对象。
    filename : str
        输出 JSON 文件路径。
    template_path : str, optional
        模板文件路径，用于补充顶层结构字段（id/compartments/version/
        enzyme_constraint）。默认使用 ./model/iZM564_SYN_etcmodel.json。
    """
    # 加载模板以获取顶层结构字段
    template = {}
    try:
        if template_path is None:
            template_path = os.path.join(
                os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
                    os.path.dirname(os.path.abspath(__file__)))))),
                'model', 'iZM564_SYN_etcmodel.json'
            )
        if os.path.exists(template_path):
            with open(template_path, 'r', encoding='utf-8') as f:
                template = json.load(f)
    except Exception:
        template = {}

    data = {
        "metabolites": _serialize_metabolites(model.metabolites),
        "reactions": _serialize_reactions(model.reactions),
        "genes": _serialize_genes(model.genes),
        "environment": _serialize_environment(model.environment),
        "additional_properties": {
            "obj_name": model.obj_name,
            "medium": model.medium,
            "K_value": model.K_value,
            "E_total": model.E_total,
            "solver": model.solver,
        }
    }

    # 自动补全 objective_coefficient：若所有反应均为空，则在含"BIOMASS"的反应中设为 1.0
    if all(r.get("objective_coefficient") in (None, "", 0, 0.0)
           for r in data["reactions"]):
        biomass_assigned = False
        for r in data["reactions"]:
            if "BIOMASS" in r["id"]:
                r["objective_coefficient"] = 1.0
                biomass_assigned = True
        # 如未找到含 BIOMASS 的反应，但模型设置了 obj_name，则在 obj_name 对应反应上设
        if not biomass_assigned and model.obj_name:
            for r in data["reactions"]:
                if r["id"] == model.obj_name:
                    r["objective_coefficient"] = 1.0
                    biomass_assigned = True
                    break

    # 补充顶层结构字段（id/compartments/version/enzyme_constraint）
    for key in ("id", "compartments", "version", "enzyme_constraint"):
        if key in template:
            data[key] = template[key]

    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(data, file, indent=4, ensure_ascii=False)



def analyze_reaction_enzyme_usage(opt_etm, model) -> pd.DataFrame:
    """从 ETM 求解结果中分析每个反应的酶用量。

    计算原理（对应 flux_balance_analysis.py 第 153-154 行）：
        kcat_MW = kcat / mw * 3600000
        E_j = v_j / kcat_MW_j = v_j * mw / (kcat * 3600000)

    Parameters
    ----------
    opt_etm : Pyomo ConcreteModel
        calculate_ETM 返回的求解后的模型对象。
    model : etgems_web.Model
        原始 ETGEM 模型对象，用于读取反应的 kcat 和 mw。

    Returns
    -------
    pd.DataFrame
        包含 flux, kcat, mw, kcat_MW, enzyme_usage, enzyme_ratio 的 DataFrame，
        按 enzyme_usage 降序排列。
    """
    results = []
    total_enzyme = 0.0

    for rea_id in opt_etm.reaction:
        flux = opt_etm.reaction[rea_id].value
        if flux is None or flux <= 1e-10:
            continue

        try:
            reaction = model.get_reaction_by_id(rea_id)
        except (KeyError, TypeError):
            try:
                reaction = model.reactions.get_by_id(rea_id)
            except (KeyError, AttributeError):
                continue

        kcat = getattr(reaction, 'kcat', None)
        mw = getattr(reaction, 'mw', None)

        if kcat is None or mw is None or kcat == '' or mw == '':
            continue
        kcat = float(kcat) if not isinstance(kcat, float) else kcat
        mw = float(mw) if not isinstance(mw, float) else mw

        if kcat <= 0 or mw <= 0:
            continue

        # kcat_MW = kcat / mw * 3600000  (对应 Get_Model_Data 第154行)
        kcat_mw = kcat / mw * 3600000
        # 酶用量 E_j = v_j / kcat_MW
        enzyme_usage = flux / kcat_mw

        results.append({
            'reaction_id': rea_id,
            'flux': flux,
            'kcat': kcat,
            'mw': mw,
            'kcat_MW': kcat_mw,
            'enzyme_usage': enzyme_usage
        })

    if not results:
        return pd.DataFrame()

    result_df = pd.DataFrame(results)
    total_enzyme = result_df['enzyme_usage'].sum()
    result_df['enzyme_ratio'] = result_df['enzyme_usage'] / total_enzyme if total_enzyme > 0 else 0.0
    result_df = result_df.sort_values('enzyme_usage', ascending=False).reset_index(drop=True)

    return result_df


def validate_reaction_kcat(reaction) -> Tuple[bool, str]:
    """验证单个反应的 kcat 和 mw 参数是否有效。

    Parameters
    ----------
    reaction : etgems_web.Reaction
        待验证的反应对象。

    Returns
    -------
    (bool, str)
        (是否有效, 描述信息)
    """
    kcat = getattr(reaction, 'kcat', None)
    mw = getattr(reaction, 'mw', None)

    if kcat is None or kcat == '':
        return False, f"kcat 为空"
    try:
        kcat = float(kcat)
    except (ValueError, TypeError):
        return False, f"kcat ({kcat}) 无法转换为浮点数"

    if mw is None or mw == '':
        return False, f"mw 为空"
    try:
        mw = float(mw)
    except (ValueError, TypeError):
        return False, f"mw ({mw}) 无法转换为浮点数"

    if kcat <= 0:
        return False, f"kcat ({kcat}) 必须大于 0"
    if mw <= 0:
        return False, f"mw ({mw}) 必须大于 0"

    return True, f"kcat={kcat}, mw={mw}, kcat_MW={kcat / mw * 3600000:.2f}"


def lookup_higher_kcat(reaction, ec_kcat_max_db: dict) -> Optional[float]:
    """根据反应的 EC 编号查询 kcat 数据库，获取更高的 kcat 值。

    优先使用 E. coli 来源的 kcat 值，若无则使用全局最大值。

    Parameters
    ----------
    reaction : etgems_web.Reaction
        反应对象，需包含 annotation 中的 ec-code。
    ec_kcat_max_db : dict
        EC kcat 最大值数据库，格式为：
        {"EC_number": {"kcat_max": float, "org": str}}

    Returns
    -------
    Optional[float]
        更高的 kcat 值（单位：/h），若无可用数据则返回 None。
    """
    annotation = getattr(reaction, 'annotation', {})
    if not annotation:
        return None

    ec_code = annotation.get('ec-code', None)
    if ec_code is None:
        return None

    kcat_max_list = []

    # EC 编号可能是字符串或列表
    ec_numbers = [ec_code] if isinstance(ec_code, str) else list(ec_code)

    for ec in ec_numbers:
        if ec in ec_kcat_max_db:
            kcat_max_list.append(ec_kcat_max_db[ec]['kcat_max'])

    if not kcat_max_list:
        return None

    return max(kcat_max_list)


def _calculate_etm_objective(model, B_value: float = 0.0):
    """运行 calculate_ETM 并返回目标值。

    辅助函数，用于在调整过程中频繁调用 calculate_ETM。

    Parameters
    ----------
    model : etgems_web.Model
        ETGEM 模型对象。
    B_value : float
        热力学 MDF 约束值。

    Returns
    -------
    (Pyomo ConcreteModel, float or None)
        (求解后的模型, 目标函数值)
    """
    try:
        from .flux_balance_analysis import calculate_ETM
    except ImportError:
        from flux_balance_analysis import calculate_ETM

    try:
        opt = calculate_ETM(model, "maximize", B_value)
        obj_val = opt.obj()
        return opt, obj_val
    except Exception as e:
        return None, None


def adjust_kcat_for_ETM(model, ec_kcat_max_path: str,
                        B_value: float = 0.0,
                        max_iterations: int = 20,
                        improvement_threshold: float = 0.001,
                        top_n: int = 1,
                        target_objective: float = None,
                        save_path: str = None,
                        verbose: bool = True) -> dict:
    """迭代优化 ETM 模型中的 kcat 参数，提升模型性能。

    核心算法：
    1. 运行 calculate_ETM 获取初始解
    2. 计算酶用量 → 识别瓶颈反应（酶用量占比最高的反应）
    3. 查询 EC kcat 数据库，获取更高 kcat 值
    4. 更新模型对象的 kcat 和 kcat_MW 属性
    5. 重新求解，验证目标值提升
    6. 重复直至收敛或达到目标值

    Parameters
    ----------
    model : etgems_web.Model
        ETGEM 模型对象（将被原地修改，建议先 deepcopy）。
    ec_kcat_max_path : str
        EC kcat 最大值数据库 JSON 文件路径。
    B_value : float, default=0.0
        热力学 MDF 约束值，传递给 calculate_ETM。
    max_iterations : int, default=20
        最大迭代次数。
    improvement_threshold : float, default=0.001
        目标值相对提升阈值，低于此值视为收敛。
    top_n : int, default=1
        每次迭代调整的瓶颈反应数量（按酶用量占比排序取前 N 个）。
    target_objective : float, optional
        目标值。当优化目标（如 biomass）达到或超过此值时提前停止。
        例如 target_objective=0.4422 表示 biomass >= 0.4422 就停止。
    save_path : str, optional
        优化后模型的保存路径（JSON 格式）。提供后自动保存。
    verbose : bool, default=True
        是否打印详细日志。

    Returns
    -------
    dict
        包含以下键的调整报告：
        - 'success': bool, 是否成功完成
        - 'model': etgems_web.Model, 优化后的模型对象
        - 'initial_objective': float, 初始目标值
        - 'final_objective': float, 最终目标值
        - 'improvement': float, 相对提升比例
        - 'target_achieved': bool, 是否达到目标值
        - 'n_iterations': int, 实际迭代次数
        - 'changed_reactions': list, 被调整的反应列表
        - 'adjustment_log': list, 每次调整的详细日志
        - 'save_path': str or None, 模型保存路径
        - 'error': str, 错误信息（如有）

    Examples
    --------
    >>> model = load_etgem_model(...)
    >>> # 优化直到 biomass 达到 0.4422
    >>> result = adjust_kcat_for_ETM(model, "data/EC_kcat_max.json",
    ...                              B_value=0.0, target_objective=0.4422)
    >>> print(f"目标达成: {result['target_achieved']}")
    """
    # ─── 参数验证 ────────────────────────────────────────────────
    try:
        with open(ec_kcat_max_path, 'r') as f:
            ec_kcat_max_db = json.load(f)
    except FileNotFoundError:
        return {
            'success': False,
            'error': f"EC kcat 数据库文件未找到: {ec_kcat_max_path}"
        }
    except json.JSONDecodeError as e:
        return {
            'success': False,
            'error': f"EC kcat 数据库 JSON 解析失败: {e}"
        }

    if not isinstance(ec_kcat_max_db, dict):
        return {
            'success': False,
            'error': "EC kcat 数据库格式错误，应为 dict"
        }

    if max_iterations < 1:
        return {'success': False, 'error': "max_iterations 必须 >= 1"}
    if improvement_threshold < 0:
        return {'success': False, 'error': "improvement_threshold 必须 >= 0"}
    if top_n < 1:
        return {'success': False, 'error': "top_n 必须 >= 1"}

    # ─── 初始求解 ────────────────────────────────────────────────
    if verbose:
        print("=" * 60)
        print("kcat 参数优化 - 初始求解")
        print("=" * 60)

    opt_initial, initial_obj = _calculate_etm_objective(model, B_value)

    if initial_obj is None:
        return {
            'success': False,
            'error': "初始 calculate_ETM 求解失败"
        }

    if verbose:
        print(f"初始目标值 (biomass): {initial_obj:.6f}")

    # ─── 迭代优化 ────────────────────────────────────────────────
    changed_reactions = []
    adjustment_log = []
    current_obj = initial_obj
    best_obj = initial_obj
    target_achieved = False

    for iteration in range(1, max_iterations + 1):
        if verbose:
            print(f"\n--- 迭代 {iteration}/{max_iterations} ---")

        # 1. 分析酶用量
        opt_etm, obj_val = _calculate_etm_objective(model, B_value)

        if obj_val is None:
            if verbose:
                print("  求解失败，终止迭代")
            break

        # 2. 检查是否达到目标值（智能停止：取最接近目标值的结果）
        if target_objective is not None and obj_val >= target_objective:
            target_achieved = True
            prev_obj = current_obj  # 上一轮的目标值（当前 kcat 调整之前）
            diff_prev = abs(target_objective - prev_obj)
            diff_curr = abs(target_objective - obj_val)

            if diff_prev < diff_curr and adjustment_log:
                # 上一轮更接近目标值 → 回退最后一轮 kcat 变更
                last_changes = adjustment_log[-1]['changed_reactions']
                for change in last_changes:
                    rea_id = change['reaction_id']
                    old_kcat = change['old_kcat']
                    try:
                        reaction = model.get_reaction_by_id(rea_id)
                        reaction.set_kcat(old_kcat)
                        reaction.kcat_MW = old_kcat / float(reaction.mw) * 3600000
                    except (KeyError, TypeError, AttributeError):
                        try:
                            reaction = model.reactions.get_by_id(rea_id)
                            reaction.kcat = old_kcat
                            reaction.kcat_MW = old_kcat / float(reaction.mw) * 3600000
                        except (KeyError, AttributeError):
                            pass
                current_obj = prev_obj
                if verbose:
                    print(f"  ✓ 目标达成（回退）: final_obj={prev_obj:.6f}, "
                          f"差值 {diff_prev:.6f} < {diff_curr:.6f}")
            else:
                # 当前值更接近（或相等取较高）→ 保留当前值
                current_obj = obj_val
                if verbose:
                    print(f"  ✓ 目标达成: final_obj={obj_val:.6f}, "
                          f"差值 {diff_curr:.6f} <= {diff_prev:.6f}")
            break

        enzyme_df = analyze_reaction_enzyme_usage(opt_etm, model)

        if enzyme_df.empty:
            if verbose:
                print("  没有有效的酶用量数据，终止迭代")
            break

        # 3. 检查收敛性
        relative_improvement = (obj_val - current_obj) / abs(current_obj) if current_obj != 0 else 0
        if iteration > 1 and relative_improvement < improvement_threshold:
            if verbose:
                print(f"  收敛: 相对提升 {relative_improvement*100:.4f}% < 阈值 {improvement_threshold*100:.2f}%")
            break

        current_obj = obj_val
        if current_obj > best_obj:
            best_obj = current_obj

        # 3. 识别瓶颈反应（从酶用量排序列表中依次尝试 top_n 个）
        iteration_changed = []
        skipped_count = 0

        while len(iteration_changed) < top_n and skipped_count < min(len(enzyme_df), 50):
            start_idx = skipped_count
            candidate_pool = enzyme_df.iloc[start_idx:start_idx + top_n]

            if verbose:
                print(f"  当前目标值: {obj_val:.6f}")
                print(f"  总酶用量: {enzyme_df['enzyme_usage'].sum():.6f}")
                if skipped_count > 0:
                    print(f"  (已跳过 {skipped_count} 个无法更新的瓶颈，继续尝试...)")

            for _, row in candidate_pool.iterrows():
                if len(iteration_changed) >= top_n:
                    break

                rea_id = row['reaction_id']

                # 跳过已调整过的反应
                if rea_id in changed_reactions:
                    skipped_count += 1
                    continue

                try:
                    reaction = model.get_reaction_by_id(rea_id)
                except (KeyError, TypeError):
                    try:
                        reaction = model.reactions.get_by_id(rea_id)
                    except (KeyError, AttributeError):
                        skipped_count += 1
                        continue

                ec_code = reaction.annotation.get('ec-code', None) if hasattr(reaction, 'annotation') else None
                if ec_code is None:
                    if verbose:
                        print(f"  跳过 {rea_id}: 无 EC 编号")
                    skipped_count += 1
                    continue

                # 查询数据库中的更高 kcat
                new_kcat = lookup_higher_kcat(reaction, ec_kcat_max_db)

                if new_kcat is None:
                    if verbose:
                        print(f"  跳过 {rea_id}: 数据库中无对应 EC 编号")
                    skipped_count += 1
                    continue

                old_kcat = float(reaction.kcat) if not isinstance(reaction.kcat, float) else reaction.kcat

                if new_kcat <= old_kcat:
                    if verbose:
                        print(f"  跳过 {rea_id}: 数据库 kcat ({new_kcat:.2f}) <= 当前 kcat ({old_kcat:.2f})")
                    skipped_count += 1
                    continue

                if verbose:
                    print(f"  瓶颈: {rea_id} | "
                          f"通量={row['flux']:.4f} | "
                          f"酶占比={row['enzyme_ratio']*100:.2f}% | "
                          f"kcat={old_kcat:.2f} /h")

                # 更新 kcat（保持 mw 不变）
                try:
                    reaction.set_kcat(new_kcat)
                except (ValueError, AttributeError):
                    try:
                        reaction.kcat = new_kcat
                    except Exception as e:
                        if verbose:
                            print(f"  更新 {rea_id} 失败: {e}")
                        skipped_count += 1
                        continue

                # 更新 kcat_MW（用于 Get_Model_Data 读取）
                try:
                    reaction.kcat_MW = new_kcat / float(reaction.mw) * 3600000
                except (AttributeError, TypeError, ValueError):
                    pass

                iteration_changed.append({
                    'reaction_id': rea_id,
                    'old_kcat': old_kcat,
                    'new_kcat': new_kcat,
                    'ec_code': ec_code
                })
                changed_reactions.append(rea_id)

                if verbose:
                    print(f"  ✓ 更新 {rea_id}: kcat {old_kcat:.2f} → {new_kcat:.2f} /h "
                          f"(EC: {ec_code})")

            if iteration_changed:
                break
            else:
                # 所有候选都无法更新，跳过并继续尝试下一批
                skipped_count += top_n

        if not iteration_changed:
            if verbose:
                print("  所有瓶颈反应均无法更新，终止迭代")
            break

        adjustment_log.append({
            'iteration': iteration,
            'objective': obj_val,
            'changed_reactions': iteration_changed
        })

    # ─── 最终验证 ────────────────────────────────────────────────
    opt_final, final_obj = _calculate_etm_objective(model, B_value)

    if final_obj is None:
        return {
            'success': False,
            'error': "最终 calculate_ETM 求解失败",
            'model': model,
            'initial_objective': initial_obj,
            'changed_reactions': changed_reactions,
            'adjustment_log': adjustment_log
        }

    improvement = (final_obj - initial_obj) / abs(initial_obj) if initial_obj != 0 else 0

    # 保存模型（使用本地 write_json_model，不依赖 io/json.py）
    saved_path = None
    if save_path:
        try:
            write_json_model(model, save_path)
            saved_path = save_path
        except Exception as e:
            print(f"  保存模型失败: {e}")

    if verbose:
        print("\n" + "=" * 60)
        print("kcat 参数优化 - 完成")
        print("=" * 60)
        print(f"初始目标值:     {initial_obj:.6f}")
        print(f"最终目标值:     {final_obj:.6f}")
        print(f"相对提升:       {improvement*100:.2f}%")
        print(f"调整的反应数:   {len(set(changed_reactions))}")
        print(f"实际迭代次数:   {len(adjustment_log)}")
        if saved_path:
            print(f"模型已保存:     {saved_path}")
        print("=" * 60)

    return {
        'success': True,
        'model': model,
        'initial_objective': initial_obj,
        'final_objective': final_obj,
        'improvement': improvement,
        'target_achieved': target_achieved,
        'n_iterations': len(adjustment_log),
        'changed_reactions': list(set(changed_reactions)),
        'adjustment_log': adjustment_log,
        'save_path': saved_path
    }


def evaluate_kcat_adjustment(model_original, model_adjusted,
                             ec_kcat_max_path: str,
                             B_value: float = 0.0,
                             verbose: bool = True) -> dict:
    """评估 kcat 调整前后的性能对比。

    对原始模型和调整后的模型分别运行 calculate_ETM，
    对比目标值、酶用量分布等指标。

    Parameters
    ----------
    model_original : etgems_web.Model
        原始模型（未调整）。
    model_adjusted : etgems_web.Model
        调整后的模型。
    ec_kcat_max_path : str
        EC kcat 最大值数据库路径（用于验证）。
    B_value : float, default=0.0
        热力学 MDF 约束值。
    verbose : bool, default=True
        是否打印评估报告。

    Returns
    -------
    dict
        包含各项评估指标的字典。
    """
    # 求解原始模型
    opt_orig, obj_orig = _calculate_etm_objective(model_original, B_value)
    # 求解调整后模型
    opt_adj, obj_adj = _calculate_etm_objective(model_adjusted, B_value)

    if obj_orig is None:
        return {'error': '原始模型求解失败'}
    if obj_adj is None:
        return {'error': '调整后模型求解失败'}

    # 酶用量分析
    enzyme_orig = analyze_reaction_enzyme_usage(opt_orig, model_original)
    enzyme_adj = analyze_reaction_enzyme_usage(opt_adj, model_adjusted)

    # 汇总统计
    total_enzyme_orig = enzyme_orig['enzyme_usage'].sum() if not enzyme_orig.empty else 0
    total_enzyme_adj = enzyme_adj['enzyme_usage'].sum() if not enzyme_adj.empty else 0

    improvement = (obj_adj - obj_orig) / abs(obj_orig) if obj_orig != 0 else 0

    # 统计 kcat 变化
    kcat_changes = {}
    for rea_id in model_original.reactions:
        r_orig = model_original.get_reaction_by_id(rea_id)
        r_adj = model_adjusted.get_reaction_by_id(rea_id)
        k_orig = float(r_orig.kcat) if r_orig.kcat != '' and r_orig.kcat is not None else 0
        k_adj = float(r_adj.kcat) if r_adj.kcat != '' and r_adj.kcat is not None else 0
        if k_orig != k_adj:
            kcat_changes[rea_id] = {'old_kcat': k_orig, 'new_kcat': k_adj}

    # 酶用量占比最高的前 5 个反应变化
    top_enzyme_orig = enzyme_orig.head(5)[['reaction_id', 'enzyme_usage', 'enzyme_ratio']].to_dict('records') if not enzyme_orig.empty else []
    top_enzyme_adj = enzyme_adj.head(5)[['reaction_id', 'enzyme_usage', 'enzyme_ratio']].to_dict('records') if not enzyme_adj.empty else []

    result = {
        'original_objective': obj_orig,
        'adjusted_objective': obj_adj,
        'improvement': improvement,
        'improvement_percent': improvement * 100,
        'original_total_enzyme': total_enzyme_orig,
        'adjusted_total_enzyme': total_enzyme_adj,
        'n_reactions_changed': len(kcat_changes),
        'kcat_changes': kcat_changes,
        'top_enzyme_original': top_enzyme_orig,
        'top_enzyme_adjusted': top_enzyme_adj
    }

    if verbose:
        print("=" * 60)
        print("kcat 调整性能评估报告")
        print("=" * 60)
        print(f"目标值 (原始):     {obj_orig:.6f}")
        print(f"目标值 (调整后):   {obj_adj:.6f}")
        print(f"相对提升:          {improvement*100:.2f}%")
        print(f"总酶用量 (原始):   {total_enzyme_orig:.6f}")
        print(f"总酶用量 (调整后): {total_enzyme_adj:.6f}")
        print(f"调整的反应数量:    {len(kcat_changes)}")
        print("-" * 60)
        if kcat_changes:
            print("kcat 变化详情:")
            for rid, chg in kcat_changes.items():
                print(f"  {rid}: {chg['old_kcat']:.2f} → {chg['new_kcat']:.2f} /h")
        print("=" * 60)

    return result


def load_ec_kcat_database(ec_kcat_max_path: str) -> dict:
    """加载 EC kcat 最大值数据库并进行基础校验。

    Parameters
    ----------
    ec_kcat_max_path : str
        JSON 文件路径，格式为：
        {"EC_number": {"kcat_max": float, "org": str}}

    Returns
    -------
    dict
        校验后的数据库字典。

    Raises
    ------
    FileNotFoundError
        文件不存在时抛出。
    ValueError
        数据库格式无效时抛出。
    """
    try:
        with open(ec_kcat_max_path, 'r') as f:
            db = json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"EC kcat 数据库文件未找到: {ec_kcat_max_path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"JSON 解析失败: {e}")

    if not isinstance(db, dict):
        raise ValueError("数据库格式错误: 根对象应为 dict")

    # 验证每个条目
    invalid_entries = []
    for ec, entry in db.items():
        if not isinstance(entry, dict):
            invalid_entries.append(ec)
            continue
        if 'kcat_max' not in entry:
            invalid_entries.append(ec)
            continue
        try:
            val = float(entry['kcat_max'])
            if val <= 0:
                invalid_entries.append(ec)
        except (ValueError, TypeError):
            invalid_entries.append(ec)

    if invalid_entries:
        print(f"警告: {len(invalid_entries)} 个条目格式无效: {invalid_entries[:5]}...")

    return db


def batch_kcat_optimization_pipeline(model, ec_kcat_max_path: str,
                                     B_values: List[float] = None,
                                     max_iterations: int = 20,
                                     improvement_threshold: float = 0.001,
                                     verbose: bool = True) -> dict:
    """批量 kcat 优化流水线：对多个 B_value 分别进行优化并汇总结果。

    Parameters
    ----------
    model : etgems_web.Model
        ETGEM 模型对象。
    ec_kcat_max_path : str
        EC kcat 最大值数据库路径。
    B_values : list of float, optional
        需要测试的 B_value 列表。默认 [0.0]。
    max_iterations : int, default=20
        最大迭代次数。
    improvement_threshold : float, default=0.001
        收敛阈值。
    verbose : bool, default=True
        是否打印日志。

    Returns
    -------
    dict
        包含每个 B_value 的优化结果。
    """
    if B_values is None:
        B_values = [0.0]

    results = {}
    for B in B_values:
        if verbose:
            print(f"\n{'#' * 60}")
            print(f"# B_value = {B}")
            print(f"{'#' * 60}")

        model_copy = deepcopy(model)
        result = adjust_kcat_for_ETM(
            model_copy, ec_kcat_max_path,
            B_value=B,
            max_iterations=max_iterations,
            improvement_threshold=improvement_threshold,
            verbose=verbose
        )
        results[B] = result

    return results