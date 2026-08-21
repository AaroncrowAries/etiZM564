from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
import cobra
from cobra.util.solver import set_objective
import pyomo.environ as pyo
from pyomo.environ import *
from pyomo.opt import SolverFactory
from .flux_balance_analysis import *
import pandas as pd

import pandas as pd
def Get_Fusion_Protein(reaction_list,model):
    """Get the fusion protein of the reaction list.
    Arguments:
    ----------
    * reaction_list: The reaction list of the model.
    * model: a dictionary of model.
    Returns:
    ----------
    * result: The fusion protein of the reaction list. The GPR relationship is the same for every pair of reactions
    , with one having a positive ΔG0 and the other having a negative ΔG0.
    """
    # get the reaction with positive flux
    reaction_gpr = {}
    result = []
    # gpr of the reactions
    for reaction in model.reactions.values():
        if reaction.get_g0 < 0 and reaction.g0 != -99999:
            reaction_gpr[reaction.id] = reaction.gene_reaction_rule.split(" and ")
    reaction_gpr_filtered = {key : value for key,value in reaction_gpr.items() if value != ['']}
    gene_reaction_dict = {}
    for key, values in reaction_gpr_filtered.items():
        for value in values:
            gene_reaction_dict.setdefault(value, []).append(key)
    for reaction_id in reaction_list:
        reaction = model.get_reaction_by_id(reaction_id)
        gene_list = reaction.gene_reaction_rule.split(" and ")
        for gene in gene_list:
            if gene in gene_reaction_dict.keys():
                for reaction_result in gene_reaction_dict[gene]:
                    if(reaction_result.split("_")[0] != reaction_id.split("_")[0]):
                        result.append([reaction_id,reaction_result])
    columns = ['Bottleneck_reaction', 'Bottleneck_reaction_g0', 'Bottleneck_reaction_gpr', 'Bottleneck_reaction_equation',
           'Coupling_reaction', 'Coupling_reaction_g0', 'Coupling_reaction_gpr', 'Coupling_reaction_equation',
           'final_g0']
    result_df = pd.DataFrame(columns=columns)
    for item in result:
        Bottleneck_reaction = model.get_reaction_by_id(item[0])
        Coupling_reaction = model.get_reaction_by_id(item[1])
        # result_row = [item[0], Bottleneck_reaction.g0, Bottleneck_reaction.gene_reaction_rule, Bottleneck_reaction.build_reaction_string,]
        result_row = {'Bottleneck_reaction' : Bottleneck_reaction.id, 'Bottleneck_reaction_g0' : Bottleneck_reaction.g0,'Bottleneck_reaction_gpr' : Bottleneck_reaction.gene_reaction_rule,'Bottleneck_reaction_equation' : Bottleneck_reaction.build_reaction_string()
                                    ,'Coupling_reaction' : Coupling_reaction.id, 'Coupling_reaction_g0' : Coupling_reaction.g0, 'Coupling_reaction_gpr' : Coupling_reaction.gene_reaction_rule,'Coupling_reaction_equation' : Coupling_reaction.build_reaction_string(), 'final_g0': Coupling_reaction.g0 + Bottleneck_reaction.g0}
        result_df = result_df.append(result_row,ignore_index=True)
    # result_df['index'] = result  
    # result_df = result_df.set_index('index')
    return result_df

# def Find_Energy_Coupling_Strategy(model, model_type, target, B_value):
#     if model_type == "TCM":
#         opt = calculate_TCM(model, target, B_value)
#     if model_type == "ETM":
#         opt = calculate_ETM(model, target, B_value)        
#     energy_list = ['atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J',
#                 'nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0',
#                 'nadh_c','NADH','C00004','C00004[c]','MNXM10','HMDB01487','cpd00004','BOPGDPNILDQYTO-NNYOXOHSSA-L','nadh','META:NADH','cpd00004_c0']
#     common_metabolites = ["h2o_c","h_c","co2_c","nadph_c","nadp_c","atp_c","adp_c","amp_c","nadh_c","o2_c","nad_c","nh4_c","cu2_c","mg2_p","cu2_p","mg2_c",
#                         "ni2_c","pi_p","zn2_p","k_p","fe2_p","ni2_p","so4_p","zn2_c","fe3_p","ppi_c","fe3_c","h_p","pi_c","k_c","so4_c","fe2_c"]
#     reactant_dict = {}
#     product_dict = {}
#     energy_reaction_dict = {}
#     energy_reaction_dict_2 = {}
#     for reaction_id in model.reactions.keys():
#         if opt.reaction[reaction_id].value != 0:
#             reaction = model.get_reaction_by_id(reaction_id)
#             for metabolite_id in reaction.reactants:
#                 if metabolite_id not in reactant_dict.keys():
#                     reactant_dict[metabolite_id] = [reaction_id]
#                 else:
#                     reactant_dict[metabolite_id].append(reaction_id)
#             for metabolite_id in reaction.products:
#                 if metabolite_id not in product_dict.keys():
#                     product_dict[metabolite_id] = [reaction_id]
#                 else:
#                     product_dict[metabolite_id].append(reaction_id)
#     for reaction_id in model.reactions.keys():
#         has_energy = False
#         reaction = model.get_reaction_by_id(reaction_id)
#         for metabolite_id in reaction.reactants:
#             if metabolite_id in energy_list:
#                 has_energy = True
#         if opt.reaction[reaction_id].value == 0 and reaction.g0 < 0 and has_energy:
#             if all(item in reactant_dict.keys() for item in reaction.reactants) and all(item in product_dict.keys() for item in reaction.products):
#                 reaction_unique_id = reaction_id.rsplit('_',1)[0]
#                 result_dict = {}
#                 for reactant_id in reaction.reactants:
#                     if reactant_id not in common_metabolites:
#                         result_dict[reactant_id] = reactant_dict[reactant_id]
#                 energy_reaction_dict[reaction_id] = result_dict
#                 for key_1,value_1 in energy_reaction_dict[reaction_id].items():
#                     new_list = [item for item in value_1 if item.rsplit('_',1)[0] != reaction_unique_id]
#                     energy_reaction_dict[reaction_id][key_1] = new_list
#                 energy_reaction_dict[reaction_id] = {key: value for key, value in energy_reaction_dict[reaction_id].items() if value}
#                 energy_reaction_dict = {key: value for key, value in energy_reaction_dict.items() if value}
#                 result_dict = {}
#                 for product_id in reaction.products:
#                     if product_id not in common_metabolites:
#                         result_dict[product_id] = product_dict[product_id]
#                 energy_reaction_dict_2[reaction_id] = result_dict
#                 for key_1,value_1 in energy_reaction_dict_2[reaction_id].items():
#                     new_list = [item for item in value_1 if item.rsplit('_',1)[0] != reaction_unique_id]
#                     energy_reaction_dict_2[reaction_id][key_1] = new_list
#                 energy_reaction_dict_2[reaction_id] = {key: value for key, value in energy_reaction_dict_2[reaction_id].items() if value}
#                 energy_reaction_dict_2 = {key: value for key, value in energy_reaction_dict_2.items() if value}
#     result = {key: {**energy_reaction_dict[key], **energy_reaction_dict_2[key]} for key in energy_reaction_dict if key in energy_reaction_dict_2}
#     return result