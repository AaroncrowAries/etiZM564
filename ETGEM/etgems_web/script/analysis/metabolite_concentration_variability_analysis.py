import pandas as pd
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
from .flux_balance_analysis import *
from concurrent.futures import ProcessPoolExecutor, as_completed

def Get_Max_Min_Met_Concentration(model,met_id,obj_target,B_value,max_biomass_under_mdf,biomass_id,Bottleneck_reaction_list):
    tem_obj = model.obj_name
    model.set_obj_name(met_id)
    Concretemodel = Template_Concretemodel(model=model,obj_target=obj_target,set_obj_Met_value=True,set_metabolite=True,set_Df=True,set_metabolite_ratio=True,set_bound=True,\
        set_stoi_matrix=True,set_substrate_ini=True,set_thermodynamics=True,B_value=B_value,\
        set_biomass_ini=True,biomass_id=biomass_id,biomass_value=max_biomass_under_mdf,set_integer=True,set_enzyme_constraint=True,\
        Bottleneck_reaction_list=Bottleneck_reaction_list,set_Bottleneck_reaction=True)
    """Calculation of the maximum and minimum concentrations for metabolites in a specific list.

    Notes：
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
    
    * reaction_kcat_MW: The enzymatic data of kcat (turnover number, /h) divided by MW (molecular weight, kDa).
    The format is as follows (in .csv file):
    kcat,MW,kcat_MW
    AADDGT,49389.2889,40.6396,1215.299582180927
    GALCTND,75632.85923999999,42.5229,1778.63831582512
    ...
    
    * lb_list: Lower bound of reaction flux. It is obtained by analyzing SBML model with get_data(model) function.
    * ub_list: Upper bound of reaction flux. It is obtained by analyzing SBML model with get_data(model) function.
    * obj_name: Object reaction ID.
    * obj_target: Type of object function (maximize or minimize).   
    * set_obj_Met_value: The object function is the concentration of a metabolite (True or False).
    * set_metabolite: Set the upper and lower bounds of metabolite concentration (True or False).
    * set_Df: Adding thermodynamic driving force expression for reactions (True or False).
    * set_metabolite_ratio: Adding concentration ratio constraints for metabolites (True or False).
    * set_bound: Set the upper and lower bounds of reaction flux (True or False).
    * set_stoi_matrix: Adding flux balance constraints (True or False).
    * set_substrate_ini: Adding substrate amount constraints (True or False).
    * K_value: the maximum value minus the minimum value of reaction thermodynamic driving force (1249, maxDFi-minDFi).
    * substrate_name: Substrate input reaction ID in the model (such as EX_glc__D_e_reverse).
    * substrate_value: Set the upper bound for substrate input reaction flux (10 mM).
    * set_thermodynamics: Adding thermodynamic constraints (True or False).
    * B_value: The value of maximizing the minimum thermodynamic driving force (MDF).
    * set_biomass_ini: Set the lower bound for biomass synthesis reaction flux (True or False).
    * biomass_id: Biomass synthesis reaction ID in the model (BIOMASS_Ec_iML1515_core_75p37M).
    * biomass_value: The lower bound of biomass synthesis reaction flux.
    * set_enzyme_constraint: Adding enzymamic constraints (True or False).
    * set_integer: Adding binary variables constraints (True or False)
    * E_total: Total amount constraint of enzymes (0.114).
    * Bottleneck_reaction_list: A list extracted from the result file automatically.
    * set_Bottleneck_reaction: Adding integer variable constraints for specific reaction (True or False).
    """  
    max_min_Df_list=pd.DataFrame()
    try:
        opt=Model_Solve(Concretemodel,model.solver)   
    except:
        print('Model can not solve!')
    else: 
        if obj_target=='maximize':
            max_min_Df_list.loc[model.obj_name,'max_value']=opt.obj() 
        elif obj_target=='minimize':
            max_min_Df_list.loc[model.obj_name,'min_value']=opt.obj() 
    model.set_obj_name(tem_obj)
    return max_min_Df_list

def Get_Max_Min_Met_Concentration_for_list(model,B_value,max_biomass_under_mdf,biomass_id,Bottleneck_reaction_list,Bottleneck_reaction_met,max_workers):
    max_min_concentration_list_fixed = pd.DataFrame()

    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(Get_Max_Min_Met_Concentration,model,eachmet,'maximize',B_value,\
            max_biomass_under_mdf,biomass_id,Bottleneck_reaction_list): eachmet for eachmet in Bottleneck_reaction_met}
        for future in as_completed(futures):
            tmp = future.result()
            for eachindex in tmp.index:
                #print(eachindex,tmp.loc[eachindex,'max_value'])
                max_min_concentration_list_fixed.loc[eachindex,'max_concentration'] = tmp.loc[eachindex,'max_value']

    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(Get_Max_Min_Met_Concentration,model,eachmet,'minimize',B_value,\
            max_biomass_under_mdf,biomass_id,Bottleneck_reaction_list): eachmet for eachmet in Bottleneck_reaction_met}
        for future in as_completed(futures):
            tmp = future.result()
            for eachindex in tmp.index:
                #print(eachindex,tmp.loc[eachindex,'max_value'])
                max_min_concentration_list_fixed.loc[eachindex,'min_concentration'] = tmp.loc[eachindex,'min_value']
    return max_min_concentration_list_fixed