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

def Get_Continuous_Reaction(reaction_list,model):
    """Get the continuous reaction of the reaction list.
    Arguments:
    ----------
    * reaction_list: The reaction list of the model.
    * model: a dictionary of model.
    Returns:
    ----------
    * result: The continuous reaction of the reaction list. The product of the previous reaction is consumed by the next reaction.
    These two reactions have opposite signs for ΔG0, with one being positive and the other being negative.
    """
    reaction_negative_g0_list = []
    result = []
    metabolites_standard = ["h_c","h2o_c","co2_c","atp_c","adp_c","pi_c","nadph_c","nadp_c","o2_c","nad_c","nadh_c","so4_c","coa_c","accoa_c"]
    # reaction_g0 = pd.read_csv(reaction_g0_file,index_col=0,sep='\t')
    # get the reaction with positive flux and positive g0
    for reaction in model.reactions.values():
        if reaction.g0 < 0 and reaction.g0 != -99999:
            reaction_negative_g0_list.append(reaction.id)
    # get the coefficient of main metabolites
    for reaction_id in reaction_negative_g0_list:
        reaction = model.get_reaction_by_id(reaction_id)
        reaction_metabolites_negative_g0 = {}
        for key,value in reaction.metabolites.items():
            if key not in metabolites_standard:
                reaction_metabolites_negative_g0[key] = value
        for reaction_obj_id in reaction_list:
            reaction_obj = model.get_reaction_by_id(reaction_obj_id)
            for key, value in reaction_obj.metabolites.items():
                if key in reaction_metabolites_negative_g0.keys() and reaction_metabolites_negative_g0[key]*value < 0:
                    if reaction.id.split("_")[0] !=  reaction_obj_id.split("_")[0]:
                        result.append([reaction_obj_id, reaction.id, key])
    columns = ['Bottleneck_reaction', 'Bottleneck_reaction_g0', 'Bottleneck_reaction_gpr', 'Bottleneck_reaction_equation',
           'Coupling_reaction', 'Coupling_reaction_g0', 'Coupling_reaction_gpr', 'Coupling_reaction_equation',
           'Coupling_metabolite']
    result_df = pd.DataFrame(columns=columns)
    for item in result:
        Bottleneck_reaction = model.get_reaction_by_id(item[0])
        Coupling_reaction = model.get_reaction_by_id(item[1])
        # result_row = [item[0], Bottleneck_reaction.g0, Bottleneck_reaction.gene_reaction_rule, Bottleneck_reaction.build_reaction_string,]
        result_row = {'Bottleneck_reaction' : Bottleneck_reaction.id, 'Bottleneck_reaction_g0' : Bottleneck_reaction.g0,'Bottleneck_reaction_gpr' : Bottleneck_reaction.gene_reaction_rule,'Bottleneck_reaction_equation' : Bottleneck_reaction.build_reaction_string()
                                    ,'Coupling_reaction' : Coupling_reaction.id, 'Coupling_reaction_g0' : Coupling_reaction.g0, 'Coupling_reaction_gpr' : Coupling_reaction.gene_reaction_rule,'Coupling_reaction_equation' : Coupling_reaction.build_reaction_string(), 'Coupling_metabolite': item[2]}
        result_df = result_df.append(result_row,ignore_index=True)
    # result_df['index'] = result  
    # result_df = result_df.set_index('index')
    return result_df

def Find_Energy_Coupling_Strategy(reaction_list,model):
    energy_list = ['atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J',
                'nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0',
                'nadh_c','NADH','C00004','C00004[c]','MNXM10','HMDB01487','cpd00004','BOPGDPNILDQYTO-NNYOXOHSSA-L','nadh','META:NADH','cpd00004_c0']
    init_result = Get_Continuous_Reaction(reaction_list,model)
    continue_result = []
    final_result = []
    for index,row in init_result.iterrows():
        continue_result.append([row["Bottleneck_reaction"],row["Coupling_reaction"],row["Coupling_metabolite"]])
    for continue_info in continue_result:
        Bottleneck_reaction = model.get_reaction_by_id(continue_info[0])
        Coupling_reaction = model.get_reaction_by_id(continue_info[1])
        for metabolite in Bottleneck_reaction.reactants:
            if metabolite in energy_list:
                if metabolite in Coupling_reaction.products:
                    continue_info.append(metabolite)
                    final_result.append(continue_info)
        for metabolite in Bottleneck_reaction.products:
            if metabolite in energy_list:
                if metabolite in Coupling_reaction.reactants:
                    continue_info.append(metabolite)
                    final_result.append(continue_info)
    columns = ['Bottleneck_reaction', 'Bottleneck_reaction_g0', 'Bottleneck_reaction_gpr', 'Bottleneck_reaction_equation',
           'Coupling_reaction', 'Coupling_reaction_g0', 'Coupling_reaction_gpr', 'Coupling_reaction_equation',
           'Coupling_metabolite']
    result_df = pd.DataFrame(columns=columns)
    for item in final_result:
        Bottleneck_reaction = model.get_reaction_by_id(item[0])
        Coupling_reaction = model.get_reaction_by_id(item[1])
        # result_row = [item[0], Bottleneck_reaction.g0, Bottleneck_reaction.gene_reaction_rule, Bottleneck_reaction.build_reaction_string,]
        result_row = {'Bottleneck_reaction' : Bottleneck_reaction.id, 'Bottleneck_reaction_g0' : Bottleneck_reaction.g0,'Bottleneck_reaction_gpr' : Bottleneck_reaction.gene_reaction_rule,'Bottleneck_reaction_equation' : Bottleneck_reaction.build_reaction_string()
                                    ,'Coupling_reaction' : Coupling_reaction.id, 'Coupling_reaction_g0' : Coupling_reaction.g0, 'Coupling_reaction_gpr' : Coupling_reaction.gene_reaction_rule,'Coupling_reaction_equation' : Coupling_reaction.build_reaction_string(), 'Coupling_metabolite': item[2], 'Coupling_energy':item[3]}
        result_df = result_df.append(result_row,ignore_index=True)
    return result_df
