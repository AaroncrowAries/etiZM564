from ..core.model import Model 
from ..core.environment import Environment 
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.gene import Gene
from .g0 import *
import pandas as pd
import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from ..analysis.flux_balance_analysis import *

def get_K_value(model):
    max_min_Df_list = pd.DataFrame()
    starttime = datetime.datetime.now()
    filtered_reactions = []
    obj = model.obj_name
    # print(obj)
    for reaction in model.reactions.values():
    # for reaction in model.reactions:
        if abs(reaction.g0) != 0 and abs(reaction.g0) != 99999:
            filtered_reactions.append(reaction.id)
    # print(filtered_reactions)
    # Maximize (DFi), without metabolites ratio constraint
    with ProcessPoolExecutor() as executor:
        print('maximize')
        futures = {executor.submit(Get_Max_Min_Df, model, each_reaction, 'maximize', 'cplex_direct'): each_reaction for each_reaction in filtered_reactions}
        for future in as_completed(futures):
            tmp = future.result()
            # print(tmp)
            for each_index in tmp.index:
                max_min_Df_list.loc[each_index, 'max_value'] = tmp.loc[each_index, 'max_value']

    # Minimize (DFi), without metabolites ratio constraint
    with ProcessPoolExecutor() as executor:
        print('minimize')
        futures = {executor.submit(Get_Max_Min_Df, model, each_reaction, 'minimize', 'cplex_direct'): each_reaction for each_reaction in filtered_reactions}
        for future in as_completed(futures):
            tmp = future.result()
            # print(tmp)
            for each_index in tmp.index:
                max_min_Df_list.loc[each_index, 'min_value'] = tmp.loc[each_index, 'min_value']
    print(max_min_Df_list)
    max_min_Df_list['diff'] = max_min_Df_list['max_value'] - max_min_Df_list['min_value']
    K_value = np.max(max_min_Df_list['max_value']) - np.min(max_min_Df_list['min_value'])
    endtime = datetime.datetime.now()
    print(endtime - starttime)
    print("Max Df : " + str(np.max(max_min_Df_list['max_value'])))
    print("Min Df : " + str(np.min(max_min_Df_list['min_value'])))
    print("K : " + str(K_value))
    model.set_obj_name(obj)
    model.K_value = K_value
    return K_value

def Get_Max_Min_Df(model, reaction, obj_target,solver):
    model.set_obj_name(reaction)
    Concretemodel = Template_Concretemodel(model,obj_target=obj_target,\
        set_obj_TM_value=True,set_metabolite=True,set_Df=True)
    """Calculation both of maximum and minimum thermodynamic driving force for reactions,without metabolite ratio constraints.
    
    Notes:
    ----------
    * reaction_list: List of reaction IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolite_list: List of metabolite IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.

    * coef_matrix: The model coefficient matrix.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolites_lnC: Upper and lower bound information of metabolite concentration (Natural logarithm).
    The format is as follows (in .txt file):
    id	lnClb	lnCub
    2pg_c	-14.50865774	-3.912023005
    13dpg_c	-14.50865774	-3.912023005
    ...
    
    * reaction_g0: Thermodynamic parameter of reactions (drG'°) .
    The format is as follows (in .txt file):
    reaction	g0
    13PPDH2	-21.3
    13PPDH2_reverse	21.3
    ...
    
    * obj_name: Object reaction ID.
    * obj_target: Type of object function (maximize or minimize).    
    * set_obj_TM_value: set_obj_TM_value: The object function is the thermodynamic driving force of a reaction (True or False).
    * set_metabolite: Set the upper and lower bounds of metabolite concentration (True or False).
    * set_Df: Adding thermodynamic driving force expression for reactions (True or False).
    """ 
    max_min_Df_list=pd.DataFrame()  
    opt=Model_Solve(Concretemodel,solver)
    if obj_target=='maximize':
        max_min_Df_list.loc[model.obj_name,'max_value']=opt.obj() 
    elif obj_target=='minimize':
        max_min_Df_list.loc[model.obj_name,'min_value']=opt.obj() 

    return max_min_Df_list

def tcmodel_construction(model):
    """Constructs a thermodynamically consistent model.
    Arguments
    ----------
    * model:   The model to be added
    """
    get_reaction_g0(model)
    print("Successfully get reaction g0!")
    print("Calculating K value...")
    model.K_value = get_K_value(model)
    print("Successfully get K value!")
    return model

def find_bottleneck_reaction_by_id(model, reaction, obj_target, B_value):
    """To find the reactions that result in a thermodynamic model solution of zero
    Arguments:
    ----------
    * model: a dictionary of model.
    * reaction: The reaction to be tested.
    * obj_target: The target of the objective function.
    * B_value: The value of MDF(B).
    Returns:
    ----------
    * result: The reaction that results in a thermodynamic model solution of zero.
    """
    result = {}     
    if reaction.get_g0 > 0:
        tem_g0 = reaction.get_g0
        reaction.set_g0(0)
        opt_TCM = calculate_TCM(model,obj_target,B_value)
        reaction.set_g0(tem_g0)
        if opt_TCM.obj() > 0:
            result[reaction.id] = opt_TCM.obj()
    return result

def find_bottleneck_reaction(model, obj_target, B_value, max_workers):
    """To find the reactions that result in a thermodynamic model solution of zero
    Arguments:
    ----------
    * model: a dictionary of model.
    * obj_target: The target of the objective function.
    * B_value: The value of MDF(B).
    Returns:
    ----------
    * result_dir: The reaction that results in a thermodynamic model solution of zero.
    """
    result_dir = {}
    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(find_bottleneck_reaction_by_id, model, eachreaction, obj_target, B_value): 
                   eachreaction for eachreaction in model.reactions.values()}
        # futures = {executor.submit(find_bottleneck_reaction_by_id, model, eachreaction, obj_target, B_value): 
        #            eachreaction for eachreaction in model.reactions}
        for future in as_completed(futures):
            tmp=future.result()
            if tmp is not None:  
                result_dir.update(tmp)
    return result_dir

def find_double_bottleneck_reaction(model, obj_target, B_value, max_workers):
    """To find the reactions that result in a thermodynamic model solution of zero
    Arguments:
    ----------
    * model: a dictionary of model.
    * obj_target: The target of the objective function.
    * B_value: The value of MDF(B).
    Returns:
    ----------
    * result_dir: The reaction that results in a thermodynamic model solution of zero.
    """
    reaction_list = []
    double_reaction_list = []
    result_dir = {}
    # for reaction in model.reactions.values():
    for reaction in model.reactions:
        if reaction.get_g0 < 0:
            reaction_list.append(reaction)
    for i in range(len(reaction_list)):
        for j in range(i+1,len(reaction_list)):
            double_reaction_list.append([reaction_list[i],reaction_list[j]])
    def find_double_bottleneck_reaction_by_id(model, reaction, obj_target, B_value):
        result = {}     
        tem_g0_0 = reaction[0].get_g0
        tem_g0_1 = reaction[1].get_g0
        reaction[0].set_g0(0)
        reaction[1].set_g0(0)
        opt_TCM = calculate_TCM(model, obj_target, B_value)
        reaction[0].set_g0(tem_g0_0)
        reaction[1].set_g0(tem_g0_1)
        if opt_TCM.obj() > 0:
            result[reaction["id"]] = opt_TCM.obj()
        return result
    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(find_double_bottleneck_reaction_by_id, model, eachreaction, obj_target
                               , B_value): eachreaction for eachreaction in double_reaction_list}
        for future in as_completed(futures):
            tmp=future.result()
            if tmp is not None:  
                result_dir.update(tmp)
    return result_dir