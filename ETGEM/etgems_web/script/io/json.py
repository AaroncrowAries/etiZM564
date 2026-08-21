import json
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
from typing import Dict,Any,List
import re

def read_json_model_reactions(filename):
    """Read a JSON file and return a list of reactions."""
    with open(filename, 'r') as file:
        data = json.load(file)
        
    reactions = {}
    for reaction_data in data['reactions']:
        reaction = Reaction(
            id=reaction_data.get('id'),
            name=reaction_data.get('name', ''),
            subsystem=reaction_data.get('subsystem', ''),
            metabolites=reaction_data.get('metabolites', {}),
            gene_reaction_rule=reaction_data.get('gene_reaction_rule', ''),
            lower_bound=float(reaction_data.get('lower_bound', 0.0)),
            upper_bound=float(reaction_data.get('upper_bound', 1000.0)),
            kcat=float(reaction_data.get('kcat', 100000)),
            mw=float(reaction_data.get('mw', 1.0)),
            g0=float(reaction_data.get('g0', -99999.0)),
            annotation=reaction_data.get('annotation',{})
        )
        reactions[reaction.id] = reaction
    
    return reactions


def read_json_model_metabolites(filename):
    """Read a JSON file and return a list of metabolites."""
    with open(filename, 'r') as file:
        data = json.load(file)
        
    metabolites = {}
    for metabolite_data in data['metabolites']:
        metabolite = Metabolite(
            id=metabolite_data.get('id'),
            name=metabolite_data.get('name', ''),
            compartment=metabolite_data.get('compartment', ''),
            charge=int(metabolite_data.get('charge', 0)),
            formula=metabolite_data.get('formula', ''),
            concentration_ub=float(metabolite_data.get('concentration_ub', -3.912023005)),
            concentration_lb=float(metabolite_data.get('concentration_lb', -14.50865774)),
            annotation=metabolite_data.get('annotation',{})
        )
        metabolites[metabolite.id] = metabolite
    
    return metabolites

def read_json_model_genes(filename):
    """Read a JSON file and return a list of genes."""
    with open(filename, 'r') as file:
        data = json.load(file)
    genes = {}
    for gene_data in data['genes']:
        gene = Gene(
            id=gene_data.get('id'),
            name=gene_data.get('name', ''),
            uniprot = gene_data.get('uniprot', ''),
            annotation = gene_data.get("annotation",{})
        )
        genes[gene.id] = gene
    
    return genes

def read_json_model_environment(filename):
    """Read a JSON file and return an environment."""
    with open(filename, 'r') as file:
        data = json.load(file)

    environment = Environment(
        p_h=data.get('p_h', 7.5),
        p_mg=data.get('p_mg', 3.0),
        ionic_strength=data.get('ionic_strength', '0.1M'),
        temperature=data.get('temperature', '310.15K'),
    )
    
    return environment

def read_json_model(file_path):
    with open(file_path, 'r') as file:
        model_data = json.load(file)
        annotation_data = model_data["additional_properties"]
        # 解析模型数据
        genes = read_json_model_genes(file_path)
        reactions = read_json_model_reactions(file_path)
        metabolites = read_json_model_metabolites(file_path)
        environment = read_json_model_environment(file_path)
        obj_name = annotation_data.get('obj_name', '')
        medium_data = annotation_data.get('medium', {})
        K_value = annotation_data.get('K_value', 0.0)
        E_total = annotation_data.get('E_total', 0.0)
        solver = annotation_data.get('solver', 'cplex_direct')
        # 创建空的 Model 对象
        model = Model(
            genes=genes,
            reactions=reactions,
            metabolites=metabolites,
            environment=environment,
            obj_name=obj_name,
            medium=medium_data,
            K_value=K_value,
            E_total=E_total,
            solver=solver
        )
    return model
        

def write_json_model(model: Model, filename: str):
    """将模型写入 JSON 文件。"""
    data = {
        "genes": _serialize_genes(model.genes),
        "reactions": _serialize_reactions(model.reactions),
        "metabolites": _serialize_metabolites(model.metabolites),
        "environment": _serialize_environment(model.environment),
        "additional_properties": {
            "obj_name": model.obj_name,
            "medium": model.medium,
            "K_value": model.K_value,
            "E_total": model.E_total,
            "solver": model.solver
        }
    }
    
    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(data, file, indent=4, ensure_ascii=False)

def _serialize_genes(genes: Dict[str, Gene]) -> List[Dict[str, Any]]:
    """基因对象转化为字典"""
    genelist = []
    for gene in genes:
        if "uniprot" in gene.annotation.keys():
            genelist.append({"id": gene.id, "name": gene.name, "uniprot": gene.annotation["uniprot"], "annotation": gene.annotation})
        else:
            genelist.append({"id": gene.id, "name": gene.name, "uniprot": "", "annotation": gene.annotation})
        
    # return [{"id": gene.id, "name": gene.name, "uniprot": gene.uniprot, "annotation": gene.annotation} for gene in genes]
    return genelist

def _serialize_reactions(reactions: Dict[str, Reaction]) -> List[Dict[str, Any]]:
    """反应对象转化为字典"""
    serialized_reactions = []
    # for reaction in reactions.values():
    for reaction in reactions:
        serialized_reaction = {
            "id": reaction.id,
            "name": reaction.name,
            "subsystem": reaction.subsystem,
            "metabolites": {metabolit.id: value for metabolit, value in reaction.metabolites.items()},
            "gene_reaction_rule": reaction.gene_reaction_rule,
            "lower_bound": reaction.lower_bound,
            "upper_bound": reaction.upper_bound,
            "kcat": reaction.kcat,
            "mw": reaction.mw,
            "g0": reaction.g0,
            "annotation": reaction.annotation
        }
        serialized_reactions.append(serialized_reaction)
    return serialized_reactions

def _serialize_metabolites(metabolites: Dict[str, Metabolite]) -> List[Dict[str, Any]]:
    """代谢物对象转化为字典"""
    metabolitelist = []
    try:
        for metabolite in metabolites.values:
            serialize_metabolite = {"id": metabolite.id,
            "name": metabolite.name,
            "compartment": metabolite.compartment,
            "charge": metabolite.charge,
            "formula": metabolite.formula,
            "concentration_ub": metabolite.concentration_ub,
            "concentration_lb": metabolite.concentration_lb,
            "annotation": metabolite.annotation
            }
            metabolitelist.append(serialize_metabolite)
    except:
        for metabolite in metabolites:
            serialize_metabolite = {"id": metabolite.id,
            "name": metabolite.name,
            "compartment": metabolite.compartment,
            "charge": metabolite.charge,
            "formula": metabolite.formula,
            "concentration_ub": metabolite.concentration_ub,
            "concentration_lb": metabolite.concentration_lb,
            "annotation": metabolite.annotation
            }
            metabolitelist.append(serialize_metabolite)
    return metabolitelist

def _serialize_environment(environment: Environment) -> Dict[str, Any]:
    """环境对象转化为字典"""
    return {
        "p_h": environment.p_h,
        "p_mg": environment.p_mg,
        "ionic_strength": environment.ionic_strength,
        "temperature": environment.temperature
    }