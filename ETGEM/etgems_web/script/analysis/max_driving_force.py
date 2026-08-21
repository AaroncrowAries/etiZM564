import pandas as pd
import datetime
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
import cobra
from .flux_balance_analysis import *
from concurrent.futures import ProcessPoolExecutor, as_completed

def Max_MDF_By_model(model,biomass_value,biomass_id,obj_enz_constraint,obj_no_enz_constraint):
    """Calculation of the maximum MDF value under preset growth rate.
    Arguments:
    * model: The model to be analyzed.
    * substrate_name: Substrate input reaction ID in the model (such as EX_glc__D_e_reverse).
    * substrate_value: Set the upper bound for substrate input reaction flux (10 mM).
    * biomass_value: The lower bound of biomass synthesis reaction flux.
    * biomass_id: Biomass synthesis reaction ID in the model (BIOMASS_Ec_iML1515_core_75p37M).
    * K_value: the maximum value minus the minimum value of reaction thermodynamic driving force (1249, maxDFi-minDFi).
    * E_total: Total amount constraint of enzymes (0.114).
    * obj_enz_constraint: The maximum enzyme cost of a reaction.
    * obj_no_enz_constraint: The minimum enzyme cost of a reaction.
    * solver: The solver used to solve the model.
    Returns:
    * MDF_list: The maximum MDF value under preset growth rate.
    """
    MDF_list=pd.DataFrame()

    if biomass_value<=obj_no_enz_constraint:
        TCM=Template_Concretemodel(model,\
            set_substrate_ini=True,set_biomass_ini=True,\
            biomass_value=biomass_value,biomass_id=biomass_id,set_metabolite=True,set_Df=True,set_obj_B_value=True,set_stoi_matrix=True,\
            set_bound=True,set_integer=True,set_metabolite_ratio=True)

        opt=Model_Solve(TCM,model.solver)
        MDF_list.loc[biomass_value,'TCM']=opt.obj()
    else:
        MDF_list.loc[biomass_value,'TCM']=None

    if biomass_value<=obj_enz_constraint:
        ETM=Template_Concretemodel(model=model,\
            set_substrate_ini=True,set_biomass_ini=True,\
            biomass_value=biomass_value,biomass_id=biomass_id,set_metabolite=True,set_Df=True,set_obj_B_value=True,set_stoi_matrix=True,\
            set_bound=True,set_enzyme_constraint=True,set_integer=True,set_metabolite_ratio=True)

        opt=Model_Solve(ETM,model.solver)
        MDF_list.loc[biomass_value,'ETM']=opt.obj()
    else:
        MDF_list.loc[biomass_value,'ETM']=None
        
    return MDF_list

def calculate_MDF_list(norm_model, model_type, output_path):
    if model_type == "ETM":
        opt_ecm = calculate_ECM(norm_model,"maximize")
        obj_max = opt_ecm.obj()
    elif model_type == "TCM":
        opt_gem = calculate_GEM(norm_model,"maximize")
        obj_max = opt_gem.obj()
    obj_min = obj_max/10#0.1
    obj_value = obj_min
    MDF_list = []
    starttime = datetime.datetime.now()
    function_name = f"calculate_{model_type}"
    B_value = MDF_Calculation(norm_model,obj_value,norm_model.obj_name)
    if B_value !='none':
        MDF_list.append([obj_value,B_value])
    while(obj_value <= obj_max):
        obj_value += 0.01
        if obj_value <= obj_max:
            B_value = MDF_Calculation(norm_model,obj_value,norm_model.obj_name)
            if B_value !='none':
                MDF_list.append([obj_value,B_value])
                obj_value = globals()[function_name](norm_model,"maximize",B_value).obj()
                MDF_list.append([obj_value,B_value])
    endtime = datetime.datetime.now()
    print (endtime - starttime)
    MDF_list = [[round(num, 3) for num in inner_list] for inner_list in MDF_list]
    MDF_df = pd.DataFrame()
    for item in MDF_list:
        MDF_df.loc[item[0],model_type] = item[1]
    MDF_df.to_csv(output_path)
    return [MDF_list, MDF_df]

def Get_Results_Thermodynamics(model,Concretemodel):
    """The formatting of the calculated results, includes the metabolic flux, binary variable values, thermodynamic driving force of reactions, the enzyme amount and the metabolite concentrations. The value of "-9999" means that the missing of kinetic (kcat) or thermodynamickcat (drG'°) parameters.
    
    Notes:
    ----------
    * model: is in SBML format (.xml).
    * Concretemodel: Pyomo Concretemodel.
    
    * reaction_kcat_MW: The enzymatic data of kcat (turnover number, /h) divided by MW (molecular weight, kDa).
    The format is as follows (in .csv file):
    kcat,MW,kcat_MW
    AADDGT,49389.2889,40.6396,1215.299582180927
    GALCTND,75632.85923999999,42.5229,1778.63831582512
    ...
    
    * reaction_g0: Thermodynamic parameter of reactions (drG'°) .
    The format is as follows (in .txt file):
    reaction	g0
    13PPDH2	-21.3
    13PPDH2_reverse	21.3
    ...
    
    * coef_matrix: The model coefficient matrix.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolite_list: List of metabolite IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.
    """
    reaction_list,metabolite_list,lb_list,ub_list,coef_matrix,reaction_g0_list,reaction_kcat_MW_list,metabolites_lnC_list = Get_Model_Data(model)
    metabolites_lnC = pd.DataFrame.from_dict(metabolites_lnC_list, orient='index')
    metabolites_lnC.columns = ["lnClb","lnCub"]
    metabolites_lnC = metabolites_lnC.rename_axis('id')
    reaction_g0 = pd.DataFrame.from_dict(reaction_g0_list, orient='index')
    reaction_g0.columns = ["g0"]
    for index,row in reaction_g0.iterrows():
        if row['g0'] == -99999 or abs(row['g0']) == 0.0:
            reaction_g0.drop(index, axis=0, inplace=True)
    reaction_kcat_MW = pd.DataFrame.from_dict(reaction_kcat_MW_list, orient='index')
    reaction_kcat_MW.columns = ["kcat_MW"]
    result_dataframe = pd.DataFrame()
    for eachreaction in Concretemodel.reaction:
        flux=Concretemodel.reaction[eachreaction].value
        z=Concretemodel.z[eachreaction].value
        result_dataframe.loc[eachreaction,'flux']=flux
        result_dataframe.loc[eachreaction,'z']=z  
        if eachreaction in reaction_g0.index:
            result_dataframe.loc[eachreaction,'f']=-reaction_g0.loc[eachreaction,'g0']-2.579*sum(coef_matrix[i,eachreaction]*Concretemodel.metabolite[i].value  for i in metabolite_list if (i,eachreaction) in coef_matrix.keys())
        else:
            result_dataframe.loc[eachreaction,'f']=-9999
        if eachreaction in reaction_kcat_MW.index:
            result_dataframe.loc[eachreaction,'enz']= flux/(reaction_kcat_MW.loc[eachreaction,'kcat_MW'])
        else:
            result_dataframe.loc[eachreaction,'enz']= -9999 
            
        tmp=model.get_reaction_by_id(eachreaction)
        met_list=''
        for met in tmp.metabolites.items():    
            met_list=met_list+';'+str(met[0])+' : '+str(np.exp(Concretemodel.metabolite[met[0]].value))
        result_dataframe.loc[eachreaction,'met_concentration']= met_list  
        
    return(result_dataframe)

#Determination of bottleneck reactions by analysing the variability of thermodynamic driving force
def Get_Max_Min_Df_Complete(model,obj_target,path_reac_list,B_value,max_biomass_under_mdf,biomass_id):
    """Determination of bottleneck reactions by analysing the variability of thermodynamic driving force.
    Arguments:
    ----------
    * model: The model to be calculated.
    * path_reac_list: The list of reactions in the pathway.
    * K_value: The value of MDF.
    * E_total: The total enzyme mass.
    * B_value: The value of MDF(B).
    * max_biomass_under_mdf: The maximum biomass under MDF.
    * biomass_id: Biomass reaction ID.
    * substrate_name: Substrate name.
    * substrate_value: Substrate value.
    * solver: The solver used to solve the model.
    * max_workers: The maximum number of workers.
    Return:
    ----------
    * max_min_Df_list_fixed: The list of reactions in the pathway and the corresponding maximum and minimum values of MDF.
    """
    max_min_Df_list_fixed=pd.DataFrame()
    tem_obj_name = model.obj_name
    def Get_Max_Min_Df_single_Complete(model,obj_id,obj_target,B_value,max_biomass_under_mdf,biomass_id):
        model.set_obj_name(obj_id)
        Concretemodel = Template_Concretemodel(model=model,\
            obj_target=obj_target,set_obj_TM_value=True,set_metabolite=True,set_Df=True,set_metabolite_ratio=True,set_bound=True,\
            set_stoi_matrix=True,set_substrate_ini=True,set_thermodynamics=True,B_value=B_value,\
            set_biomass_ini=True,biomass_id=biomass_id,biomass_value=max_biomass_under_mdf,set_integer=True,set_enzyme_constraint=True)
        """Calculation both of maximum and minimum thermodynamic driving force for reactions in a special list.
        
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
        * ub_list[substrate_name]: substrate_value (the upper bound for substrate input reaction flux)
        * obj_name: Object reaction ID.
        * obj_target: Type of object function (maximize or minimize). 
        * set_obj_TM_value: The object function is the thermodynamic driving force of a reaction (True or False).
        * set_metabolite: Set the upper and lower bounds of metabolite concentration (True or False).
        * set_Df: Adding thermodynamic driving force expression for reactions (True or False).
        * set_metabolite_ratio: Adding concentration ratio constraints for metabolites (True or False).
        * set_bound: Set the upper and lower bounds of reaction flux (True or False).
        * set_stoi_matrix: Adding flux balance constraints (True or False).
        * set_biomass_ini: Set the lower bound for biomass synthesis reaction flux (True or False).
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
            max_min_Df_list.loc[obj_id,obj_target]=opt.obj() 
        elif obj_target=='minimize':
            max_min_Df_list.loc[obj_id,obj_target]=opt.obj() 

        return max_min_Df_list
    for reaction in path_reac_list:
        single_result = Get_Max_Min_Df_single_Complete(model, reaction, obj_target, B_value, max_biomass_under_mdf, biomass_id)
        max_min_Df_list_fixed = pd.concat([max_min_Df_list_fixed, single_result])
    model.set_obj_name(tem_obj_name)
    return max_min_Df_list_fixed

def MDF_Calculation(model,biomass_value,biomass_id):
    """Calculation of maximum thermodynamic driving force (MDF) for a model.
    Arguments:
    ----------
    * model: The model to be calculated.
    * biomass_value: Biomass value.
    * biomass_id: Biomass reaction ID.
    * substrate_name: Substrate name.
    * substrate_value: Substrate value.
    * K_value: The value of MDF.
    * E_total: The total enzyme mass.
    * solver: The solver used to solve the model.
    return:
    ----------
    * B_value: The value of MDF.
    """

    Concretemodel=Template_Concretemodel(model,\
        set_substrate_ini=True,set_biomass_ini=True,\
        biomass_value=biomass_value,biomass_id=biomass_id,set_metabolite=True,set_Df=True,set_obj_B_value=True,set_stoi_matrix=True,\
        set_bound=True,set_enzyme_constraint=True,set_integer=True,set_metabolite_ratio=True)
    opt=Model_Solve(Concretemodel,model.solver)
    try:
        opt.obj()-0.000001
    except:
        B_value='none'
    else:
        B_value=opt.obj()-0.000001
    return B_value

def Get_Max_Min_Df_Ratio(model,obj_target,solver):

    Concretemodel = Template_Concretemodel(model,obj_target=obj_target,\
        set_obj_TM_value=True,set_metabolite=True,set_Df=True,set_metabolite_ratio=True)
    """Calculation both of maximum and minimum thermodynamic driving force for reactions,with metabolite ratio constraints.
    
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
    * set_metabolite_ratio: Adding concentration ratio constraints for metabolites (True or False).   
    """
    max_min_Df_list=pd.DataFrame()
    opt = Model_Solve(Concretemodel,solver)   
    if obj_target=='maximize':
        max_min_Df_list.loc[model.obj_name,'max_value']=opt.obj() 
    elif obj_target=='minimize':
        max_min_Df_list.loc[model.obj_name,'min_value']=opt.obj() 

    return max_min_Df_list