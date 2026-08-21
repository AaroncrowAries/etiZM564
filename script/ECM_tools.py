# -*- coding: utf-8 -*-
"""
ECM_tools - 高效的酶约束代谢模型工具
整合comparison目录中notebook的公共功能
"""

import cobra
import json
import pandas as pd
from typing import List, Dict, Optional, Union
from pathlib import Path

cobra.Configuration.solver = "glpk"


class ECMTools:
    def __init__(self, model_path: str = None):
        self.model_path = model_path
        self.model = None
        self.model_dict = None
        if model_path:
            self.load_model(model_path)

    def load_model(self, model_path: str) -> cobra.Model:
        """加载酶约束模型"""
        self.model_path = model_path
        self.model = cobra.io.load_json_model(model_path)
        self.model_dict = self._json_load(model_path)
        return self.model

    @staticmethod
    def _json_load(path: str) -> dict:
        """加载JSON文件"""
        with open(path) as f:
            return json.load(f)

    def get_ec_model(self) -> cobra.Model:
        """从JSON模型文件创建酶约束模型"""
        if not self.model_dict:
            raise ValueError("Model not loaded. Call load_model first.")
        
        model = cobra.io.json.load_json_model(self.model_path)
        coefficients = {}
        
        for rxn in model.reactions:
            for eachr in self.model_dict['reactions']:
                if rxn.id == eachr['id'] and eachr.get('kcat_MW'):
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                    break
        
        constraint = model.problem.Constraint(
            0,
            lb=self.model_dict['enzyme_constraint']['lowerbound'],
            ub=self.model_dict['enzyme_constraint']['upperbound']
        )
        model.add_cons_vars(constraint)
        model.solver.update()
        constraint.set_linear_coefficients(coefficients=coefficients)
        return model

    @staticmethod
    def create_reaction(
        rxn_id: str,
        name: str,
        lower_bound: float,
        upper_bound: float,
        mets: Dict[str, float],
        kcat: float = None,
        kcat_mw: float = None
    ) -> dict:
        """创建反应字典"""
        rxn = {
            'id': rxn_id,
            'name': name,
            'lower_bound': lower_bound,
            'upper_bound': upper_bound,
            'metabolites': mets,
            'kcat': kcat if kcat else '',
            'kcat_MW': kcat_mw if kcat_mw else ''
        }
        return rxn

    def add_reactions_to_model(
        self,
        reactions: List[dict],
        bigmodel: cobra.Model = None,
        output_path: str = None
    ) -> str:
        """添加反应列表到模型并保存
        
        Args:
            reactions: 反应字典列表
            bigmodel: 提供代谢物信息的大模型（可选）
            output_path: 输出文件路径（可选）
            
        Returns:
            更新后的模型JSON路径
        """
        existing_rxn_ids = {r['id'] for r in self.model_dict['reactions']}
        existing_met_ids = {m['id'] for m in self.model_dict.get('metabolites', [])}
        
        for rxn in reactions:
            if rxn['id'] not in existing_rxn_ids:
                self.model_dict['reactions'].append(rxn)
                existing_rxn_ids.add(rxn['id'])
                
                if bigmodel:
                    for met_id in rxn['metabolites'].keys():
                        if met_id not in existing_met_ids:
                            met_info = self._get_metabolite_info(bigmodel, met_id)
                            if met_info:
                                self.model_dict['metabolites'].append(met_info)
                                existing_met_ids.add(met_id)
        
        if output_path:
            with open(output_path, 'w') as f:
                json.dump(self.model_dict, f, indent=2)
            return output_path
        return self.model_path

    def _get_metabolite_info(self, bigmodel: cobra.Model, met_id: str) -> Optional[dict]:
        """从大模型获取代谢物信息"""
        compartment_map = {'c0': 'c', 'e0': 'e', 'p0': 'p'}
        
        try:
            rep_id = met_id.replace('e0', 'c0').replace('p0', 'c0')
            metabolite = bigmodel.metabolites.get_by_id(rep_id)
            compartment = next(
                (v for k, v in compartment_map.items() if met_id.endswith(k)),
                'c'
            )
            return {
                'id': met_id,
                'name': metabolite.name,
                'compartment': compartment,
                'charge': metabolite.charge,
                'formula': str(metabolite.formula) if metabolite.formula else '',
                'annotation': metabolite.annotation
            }
        except:
            return None

    @staticmethod
    def run_pfba(
        model: cobra.Model,
        objective: str = None,
        medium: Dict[str, float] = None,
        constraints: Dict[str, tuple] = None,
        flux_threshold: float = 1e-4
    ) -> Dict[str, float]:
        """运行pFBA分析
        
        Args:
            model: COBRA模型
            objective: 目标反应ID
            medium: 培养基条件字典
            constraints: 反应约束字典 {rxn_id: (lower, upper)}
            flux_threshold: 通量阈值
            
        Returns:
            显著通量字典 {reaction_id: flux}
        """
        model = model.copy()
        
        if medium:
            current_medium = model.medium
            current_medium.update(medium)
            model.medium = current_medium
            
        if constraints:
            for rxn_id, bounds in constraints.items():
                try:
                    model.reactions.get_by_id(rxn_id).bounds = bounds
                except:
                    pass
                    
        if objective:
            model.objective = objective
            
        solution = cobra.flux_analysis.pfba(model)
        
        return {
            k: v for k, v in dict(solution.fluxes).items()
            if abs(v) > flux_threshold
        }

    @staticmethod
    def set_medium(
        model: cobra.Model,
        substrate: str = None,
        substrate_flux: float = 40,
        anaerobic: bool = True,
        **kwargs
    ) -> cobra.Model:
        """设置培养基条件
        
        Args:
            model: COBRA模型
            substrate: 底物反应ID (如 'EX_glc_D_e')
            substrate_flux: 底物吸收速率
            anaerobic: 是否厌氧
            **kwargs: 其他培养基组分
        """
        model = model.copy()
        medium = model.medium
        
        if substrate:
            substrate_reverse = substrate + '_reverse'
            medium[substrate_reverse] = substrate_flux
            medium[substrate] = 0
            
        if anaerobic:
            medium['EX_o2_e'] = 0
            
        for k, v in kwargs.items():
            medium[k] = v
            
        model.medium = medium
        return model

    @staticmethod
    def print_fluxes(fluxes: Dict[str, float], model: cobra.Model) -> None:
        """格式化打印通量结果"""
        print(f"{'Reaction':<30} {'Flux':>12}  Reaction String")
        print("-" * 80)
        for rxn_id, flux in sorted(fluxes.items(), key=lambda x: abs(x[1]), reverse=True):
            try:
                rxn_str = model.reactions.get_by_id(rxn_id).build_reaction_string(use_metabolite_names=True)
                print(f"{rxn_id:<30} {flux:>12.4f}  {rxn_str}")
            except:
                print(f"{rxn_id:<30} {flux:>12.4f}")


def quick_add_product(
    base_model_path: str,
    product_rxn: dict,
    exchange_rxn: dict,
    transport_rxn: dict = None,
    product_id: str = None,
    bigmodel_path: str = None,
    objective: str = None,
    substrate: str = 'EX_glc_D_e',
    substrate_flux: float = 40,
    anaerobic: bool = True
) -> Dict[str, float]:
    """快速添加产品途径并运行pFBA
    
    Args:
        base_model_path: 基础模型JSON路径
        product_rxn: 产物合成反应
        exchange_rxn: 交换反应
        transport_rxn: 转运反应（可选）
        product_id: 目标产物ID（用于设置objective）
        bigmodel_path: 大模型路径（提供代谢物信息）
        objective: 目标函数，默认使用产物交换反应
        substrate: 底物
        substrate_flux: 底物通量
        anaerobic: 厌氧条件
        
    Returns:
        pFBA通量结果
    """
    tools = ECMTools(base_model_path)
    
    reactions_to_add = [product_rxn, exchange_rxn]
    if transport_rxn:
        reactions_to_add.append(transport_rxn)
    
    bigmodel = None
    if bigmodel_path:
        bigmodel = cobra.io.read_sbml_model(bigmodel_path)
    
    output_path = base_model_path.replace('.json', f'_{product_id}.json')
    tools.add_reactions_to_model(reactions_to_add, bigmodel, output_path)
    
    ec_model = tools.get_ec_model()
    
    obj = objective if objective else exchange_rxn['id']
    fluxes = tools.run_pfba(
        ec_model,
        objective=obj,
        medium={substrate: substrate_flux},
        constraints={'EX_o2_e': (0, 1000)} if anaerobic else None
    )
    
    print(f"\n{'='*60}")
    print(f"Product: {product_id}")
    print(f"Yield: {fluxes.get(exchange_rxn['id'], 0):.4f}")
    print(f"{'='*60}\n")
    
    return fluxes


def compare_models(model_paths: List[str], labels: List[str] = None) -> pd.DataFrame:
    """比较多个模型的统计数据
    
    Args:
        model_paths: 模型路径列表
        labels: 模型标签列表
        
    Returns:
        包含模型统计信息的DataFrame
    """
    if labels is None:
        labels = [f"Model_{i+1}" for i in range(len(model_paths))]
    
    stats = []
    for path, label in zip(model_paths, labels):
        model = cobra.io.load_json_model(path)
        stats.append({
            'Model': label,
            'Genes': len(model.genes),
            'Reactions': len(model.reactions),
            'Metabolites': len(model.metabolites)
        })
    
    return pd.DataFrame(stats)