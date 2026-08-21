# etiZM564：运动发酵单胞菌酶-热力学约束代谢模型构建与应用

本仓库围绕 **运动发酵单胞菌（*Zymomonas mobilis* ZM4）** 开展基因组规模代谢模型（GEM）的构建、
酶约束（ecGEM）与热力学约束（ETGEM）升级，并在此基础上进行 2,3-丁二醇（2,3-BDO）等
多种高附加值产物的代谢工程改造靶标预测。

整个工作流以 Jupyter Notebook 按编号串联（`00` → `07`），核心计算依托 COBRApy，
酶约束部分使用 [ECMpy](https://github.com/tibbdc/ECMpy) 工作流与
[AutoPACMEN](https://github.com/klamt-lab/autopacmen) 动力学参数数据库，
热力学约束使用随仓库附带的 `ETGEM/etgems_web` 工具包实现。

## 模型体系

| 模型 | 全称 | 说明 |
| --- | --- | --- |
| **iZM564** | 基因组规模代谢模型 | 由 KEGG 注释 Excel（1670 反应 / 1641 代谢物 / 564 基因）重建，FBA / pFBA 求解并通过 Biolog 底物利用测试验证 |
| **eciZM564** | 酶约束模型 | 在 iZM564 上加入酶约束（kcat / 蛋白分子量 / 酶分配比例 f），经 kcat 迭代校准，使 MM / MM0.15 培养基下的预测生长速率与实验值（0.36 / 0.1533 h⁻¹）吻合 |
| **ETiZM564** | 酶-热力学约束模型（ETGEM） | 在酶约束基础上补充热力学参数（反应吉布斯自由能 g0、代谢物浓度上下限），构建五个模型变体（基础 / 23BDO / SYN / MM / MM015），可用于 `calculate_ETM`（热力学 + 酶约束 FBA） |
| **etiZM564_* / ETiZM564_*_23BDO** | 工程/添加外源反应模型 | 在 ETGEM 上添加外源反应或产物路径（如 2,3-BDO、乙二醇、1,3-PDO 等），用于产物合成能力评估 |

## 仓库结构

```
.
├── 00_model_build&Biolog_test_optimized.ipynb     # ① iZM564 重建 + Biolog 底物利用测试（C/N/P/S）
├── 01_model_iZM564_modification.ipynb             # ② 模型格式修正 + GPR 修改 + UniProt 注释
├── 02_construct_eciZM564_optimized.ipynb          # ③ 酶约束模型 eciZM564 构建与 kcat 校准
├── 03_ETGEM模型从头构建_optimized.ipynb           # ④ 五个变体的 ETGEM（酶-热力学约束）构建
├── 03-2_ETGEM_biomass&product copy.ipynb          #     ETGEM 生物量 & 乙醇产物分析
├── 03-3_ETGEM_add_rxns.ipynb                      #     向 ETGEM 添加外源 2,3-BDO 反应（微氧条件产量最高）
├── 04_ETGEM_comparison copy.ipynb                 # ⑤ ETGEM 鲁棒性分析（底物摄取 0~70）与乙醇对比
├── 05_13C_flux_model_comparison copy.ipynb        # ⑥ 13C 通量模型对比（etiZM564 vs 已发表 eciZM547）
├── 06_23BDO_FVSEOF copy.ipynb                     # ⑦ 2,3-BDO FVSEOF 菌株设计（含氧气反馈约束）
├── 07_etiZM564_multi_product.ipynb                # ⑧ 多产物合成能力评估（C2/C3/C4/C5）
│
├── ETGEM/
│   └── etgems_web/                                # 酶-热力学约束建模工具包（随仓库提供）
│       └── script/
│           ├── core/                              # 基础类：Model / Reaction / Metabolite / Gene / Environment
│           ├── ecmpy/                             # ECMpy 酶约束工作流（ecGEM 构建、AutoPACMEN kcat）
│           ├── tcmodel/                           # 热力学建模（g0、K 值、ΔfG 热力学驱动力）
│           ├── cave/                              # SPAV 通量模拟
│           ├── io/                                # 模型格式转换与 JSON 读写
│           ├── analysis/                          # FBA / FVA / kcat 优化 / BeReTa / MDF 等分析
│           └── plotting/                          # 绘图工具
│
├── script/                                        # 顶层本地脚本（ECMpy_function.py、AutoPACMEN_function.py）
├── data/                                          # 输入数据：kcat 库、GPR 修改表、蛋白丰度、亚基数等
├── model/                                         # 输出模型（iZM564 / ETiZM564 各变体 JSON）
├── analysis/                                      # 通量与 kcat 校准过程文件
├── ZM4_get_data/                                  # AutoPACMEN kcat 数据库产物（BRENDA / SABIO-RK 合并）
├── _cache/ncbi_taxonomy/                          # NCBI 分类学本地缓存
│
├── iZM564_20260807_KEGG_met.xls                   # iZM564 模型源数据（KEGG 注释）
├── Supplementary_file_1_iZM516.sbml               # 参考模型 iZM516（SBML）
├── iZM4_478-pcbi.1008137.s002.xml                 # 参考模型 iZM4-478（SBML）
├── eciZM547.json                                  # 已发表酶约束模型（用于 13C 通量对比）
└── iZM4_478.xlsx                                  # 参考模型数据
```

## 主流程（Notebook 编号顺序）

1. **00 · iZM564 重建与验证**：从 KEGG 注释 Excel 读取代谢物 / 反应 / 基因三张工作表，
   逐项加入 COBRApy 模型；FBA / pFBA 分析（无氧生物量 ≈ 0.047）；质量 / 电荷守恒检查；
   导出 SBML / JSON；以 pFBA 逐一模拟 Biolog 表型芯片，测试碳 / 氮 / 磷 / 硫源可利用性。
2. **01 · 模型修改**：统一代谢物与反应 ID 格式（满足 AutoPACMEN / ECMpy 输入要求）、
   按 `gpr_modifications.tsv` / `reaction_modifications.tsv` 修正 GPR 并删除指定反应、
   通过 ZM4 基因 ID 转换表添加 UniProt 注释（作为动力学参数获取基础）。
3. **02 · eciZM564 构建**：基于 ECMpy + AutoPACMEN 获取反应级 kcat / 蛋白分子量；
   由蛋白组丰度计算酶分配比例 f（三次重复 ≈ 0.51）；构建酶约束模型；
   以 BRENDA + 手动覆盖校准 kcat，按 MM / MM0.15 培养基迭代拟合至目标生长速率。
4. **03 · ETGEM 构建**：将五个 iZM564 变体升级为酶-热力学约束模型
   （补充 kcat / kcat_MW / mw / g0 / 代谢物浓度上下限），并针对 MM / MM015 变体做 kcat 迭代优化。
5. **03-2 / 03-3**：ETGEM 生物量与乙醇产物分析；向 ETGEM 添加外源 2,3-BDO 反应，
   模拟表明微氧条件下产量最高。
6. **04 · 对比与鲁棒性**：厌氧条件下底物摄取 0~70 mmol/gDW/h 的生长鲁棒性分析，乙醇对比。
7. **05 · 13C 通量对比**：将 etiZM564 通量结果与已发表的 eciZM547 模型对比。
8. **06 · FVSEOF 菌株设计**：对 2,3-BDO 模型执行含氧气反馈约束的 FVSEOF 扫描，
   按通量倍变趋势识别上调（建议过表达）/ 下调（建议敲除）靶标基因。
9. **07 · 多产物评估**：在 ETGEM 上添加外源产物路径，评估
   C2（乙醇、乙二醇）、C3（乳酸、1,3-PDO、丙氨酸）、C4（异丁醇、1,4-BDO、琥珀酸）、
   C5（衣康酸、木糖醇 / 木糖酸）等产物的理论合成能力。

## 环境与依赖

开发验证环境：**Python 3.9**，建议按各 Notebook 头部单元格安装依赖。

| 包 | 用途 | 版本（已验证） |
| --- | --- | --- |
| cobra | 代谢建模与 FBA | 0.22.1 |
| pandas | 表格数据处理 | 1.5.3 |
| numpy | 数值计算 | 1.23.5 |
| gurobipy | 线性规划求解器（推荐） | 12.0.3 |
| xlrd | 读取 .xls 格式 Excel | 2.0.1 |
| pyomo | 热力学 / ETM 优化建模 | — |
| requests | UniProt API 调用 | — |

> **求解器说明**：FBA 部分建议使用 Gurobi（数值稳定、结果可复现），缺失时 COBRApy 会自动回退；
> 热力学 ETM 计算（`calculate_ETM`）基于 pyomo，默认使用 `cplex_direct` 求解器。

## 快速开始

```bash
# 1. 安装依赖（以 pip 为例）
pip install cobra pandas numpy gurobipy xlrd pyomo requests

# 2. 按编号顺序运行 Notebook：00 → 01 → 02 → 03 → ...
# 建议使用 Jupyter / VS Code 逐单元格执行
```

关键可调参数均集中在各 Notebook 的配置单元格（如 `CONFIG` 字典、`EXCEL_FILE` / `MODEL_NAME` 等），
更换数据文件时只需修改对应配置。

## 数据说明

- `data/`：GPR 修改表（`gpr_modifications.tsv`）、反应删除表（`reaction_modifications.tsv`）、
  蛋白组丰度（`ZM4_1/2/3.tsv`、`iZM564_Protein_Abundance.xlsx`）、基因-亚基映射、EC kcat 上限、手动 kcat 覆盖等。
- `ZM4_get_data/`：AutoPACMEN 产物——BRENDA / SABIO-RK / UniProt 合并的 kcat 数据库、
  反应级 kcat 与蛋白分子量表（`iZM564_mean_reaction_kcat_MW.csv` 等）。
- `model/`：各阶段产物模型 JSON（`iZM564_20260807_KEGG_met.json`、`ETiZM564_MM_optimized.json`、
  `ETiZM564_MM_optimized_23BDO.json` 等）。

## 参考工具

- [ECMpy](https://github.com/tibbdc/ECMpy)：酶约束代谢模型构建与 kcat 校准工作流
- [AutoPACMEN](https://github.com/klamt-lab/autopacmen)：基于 BRENDA / SABIO-RK / UniProt 的 kcat 与蛋白分子量数据库构建
- COBRApy：基因组规模代谢建模（FBA / pFBA / 通量分析）
- ETGEMs（`ETGEM/etgems_web`）：酶-热力学约束建模工具包（热力学驱动力、K 值、BeReTa 等）
