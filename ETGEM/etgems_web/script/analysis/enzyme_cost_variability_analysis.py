from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
from .flux_balance_analysis import *
from concurrent.futures import ProcessPoolExecutor, as_completed

#Determination of key enzymes by analysing the enzyme cost variability
def Get_Max_Min_E(model,obj_name,obj_target,B_value,biomass_id,max_biomass_under_mdf):
    tem_obj = model.obj_name
    model.set_obj_name(obj_name)
    Concretemodel = Template_Concretemodel(model=model,obj_target=obj_target,set_obj_single_E_value=True,set_metabolite=True,set_Df=True,set_metabolite_ratio=True,set_bound=True,\
        set_stoi_matrix=True,set_substrate_ini=True,set_thermodynamics=True,B_value=B_value,\
        set_biomass_ini=True,biomass_id=biomass_id,biomass_value=max_biomass_under_mdf,set_integer=True,set_enzyme_constraint=True)
    """Calculation of the maximum and minimum enzyme cost for reactions in a specific list.

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
    * set_obj_single_E_value: The object function is the enzyme cost of a reaction (True or False).
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
    * set_integer: Adding binary variables constraints (True or False).
    * set_enzyme_constraint: Adding enzymamic constraints (True or False).
    * E_total: Total amount constraint of enzymes (0.114).
    """  

    max_min_Df_list=pd.DataFrame()
    opt=Model_Solve(Concretemodel,model.solver)
    if obj_target=='maximize':
        max_min_Df_list.loc[model.obj_name,'max_value']=opt.obj() 
    elif obj_target=='minimize':
        max_min_Df_list.loc[model.obj_name,'min_value']=opt.obj() 
    model.set_obj_name(tem_obj)
    return max_min_Df_list

def max_Max_Min_E_for_list(model,B_value,biomass_id,max_biomass_under_mdf,enz_use_reaction_list,max_workers):

    max_min_E_list_fixed = pd.DataFrame()

    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(Get_Max_Min_E,model,eachreaction,'maximize',B_value,\
                        biomass_id,max_biomass_under_mdf): eachreaction for eachreaction in enz_use_reaction_list}
        for future in as_completed(futures):
            tmp = future.result()
            for eachindex in tmp.index:
                #print(eachindex,tmp.loc[eachindex,'max_value'])
                max_min_E_list_fixed.loc[eachindex,'max_E']=tmp.loc[eachindex,'max_value']

    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(Get_Max_Min_E,model,eachreaction,'minimize',B_value,\
                        biomass_id,max_biomass_under_mdf): eachreaction for eachreaction in enz_use_reaction_list}
        for future in as_completed(futures):
            tmp = future.result()
            for eachindex in tmp.index:
                #print(eachindex,tmp.loc[eachindex,'max_value'])
                max_min_E_list_fixed.loc[eachindex,'min_E']=tmp.loc[eachindex,'min_value']
    return max_min_E_list_fixed

#Minimum enzyme cost calculation
def Min_Enzyme_Cost_Calculation(model,biomass_value,biomass_id,B_value):
    """Calculation of minimum enzyme cost for a model.
    Arguments:
    ----------
    * model: The model to be calculated.
    * biomass_value: Biomass value.
    * biomass_id: Biomass reaction ID.
    * substrate_name: Substrate name.
    * substrate_value: Substrate value.
    * K_value: The value of MDF.
    * E_total: The total enzyme mass.
    * B_value: The value of MDF(B).
    * solver: The solver used to solve the model.
    Return:
    ----------
    * min_E: The value of minimum enzyme cost.
    """
    Concretemodel=Template_Concretemodel(model=model,\
        set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,\
        set_metabolite=True,set_Df=True,set_stoi_matrix=True,\
        set_bound=True,set_enzyme_constraint=True,set_integer=True,set_metabolite_ratio=True,\
        set_thermodynamics=True,B_value=B_value,set_obj_E_value=True)
    opt=Model_Solve(Concretemodel,model.solver)
    min_E=opt.obj()
    return min_E