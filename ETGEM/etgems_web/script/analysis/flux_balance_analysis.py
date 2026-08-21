import json
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
import cobra
from typing import Dict,Any
from cobra.util.solver import set_objective
import re
import pandas as pd
import matplotlib.pyplot as plt
import pyomo.environ as pyo
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np

def Get_Model_Data(model):
    """Returns reaction_list,metabolite_list,lb_list,ub_list,coef_matrix from model.
    
    Notes: 
    ----------
    *model： is in SBML format (.xml).
    """
    reaction_list=[]
    metabolite_list=[]
    reaction_g0_list = {}
    reaction_kcat_MW_list = {}
    metabolites_lnc = {}
    lb_list={}
    ub_list={}
    coef_matrix={}
    try:
        for rea in model.reactions.values():
            reaction_list.append(rea.id)
            if rea.lower_bound + rea.upper_bound < 0:
                lb_list[rea.id] = rea.lower_bound + 1000
                ub_list[rea.id] = rea.upper_bound + 1000
                try:
                    for met in model.metabolites.values():
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                            # print(rea.get_coefficient(met.id))
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=-rea.metabolites[met.id]
                            # print(-rea.get_coefficient(met.id))
                except:
                    for met in model.metabolites:
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                            # print(rea.get_coefficient(met.id))
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=-rea.metabolites[met.id]
                            # print(-rea.get_coefficient(met.id))
            else:
                lb_list[rea.id]=rea.lower_bound
                ub_list[rea.id]=rea.upper_bound
                try:
                    for met in model.metabolites.values():
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=rea.metabolites[met.id]
                except:
                    for met in model.metabolites:
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=rea.metabolites[met.id]
    except:
        for rea in model.reactions:
            reaction_list.append(rea.id)

            if rea.lower_bound + rea.upper_bound < 0:
                lb_list[rea.id] = rea.lower_bound + 1000
                ub_list[rea.id] = rea.upper_bound + 1000
                try:
                    for met in model.metabolites.values():
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                            # print(rea.get_coefficient(met.id))
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=-rea.metabolites[met.id]
                            # print(-rea.get_coefficient(met.id))
                except:
                    for met in model.metabolites:
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                            # print(rea.get_coefficient(met.id))
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=-rea.metabolites[met.id]
                            # print(-rea.get_coefficient(met.id))
            else:
                lb_list[rea.id]=rea.lower_bound
                ub_list[rea.id]=rea.upper_bound
                try:
                    for met in model.metabolites.values():
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=rea.metabolites[met.id]
                except:
                    for met in model.metabolites:
                        metabolite_list.append(met.id)
                        try:
                            rea.metabolites[met.id]
                        except:
                            pass
                        else:
                            coef_matrix[met.id,rea.id]=rea.metabolites[met.id]
    reaction_list=list(set(reaction_list))
    metabolite_list=list(set(metabolite_list))
    for rea in reaction_list:
        # print(rea)
        try:
            reaction = model.get_reaction_by_id(rea)
        except KeyError:
            # 如果rea不是有效的ID，则尝试另一种方法（假设rea可能是一个对象）
            try:
                # 检查rea是否有id属性，并尝试通过该属性获取反应对象
                reaction = model.reactions.get_by_id(rea)
            except:
                reaction = model.reactions[rea.id]
        except TypeError:
            # 如果rea不是预期的字符串或对象类型，则捕获TypeError（此异常类型可能需要根据实际情况调整）
            try:
                # 检查rea是否有id属性，并尝试通过该属性获取反应对象
                reaction = model.reactions.get_by_id(rea)
            except:
                reaction = model.reactions[rea.id]
        reaction_g0_list[rea] = reaction.g0
        if reaction.kcat != "" and reaction.mw != "":
            reaction_kcat_MW_list[rea] = reaction.kcat/reaction.mw*3600000
        else:
            reaction_kcat_MW_list[rea] = 0
    for met in metabolite_list:
        try:
            metabolite = model.get_metabolite_by_id(met)
            # print(met)
        except KeyError:
            # 如果met不是有效的ID，则尝试另一种方法（假设met可能是一个对象）
            try:
                # 检查met是否有id属性，并尝试通过该属性获取反应对象
                metabolite = model.metabolites.get_by_id(met)
                # print(met.id)
            except:
                metabolite = model.metabolites[met.id]  
                # print(met.id)
        except TypeError:
            try:
                # 检查met是否有id属性，并尝试通过该属性获取反应对象
                metabolite = model.metabolites.get_by_id(met)
            except:
                metabolite = model.metabolites[met.id]  
            
        metabolites_lnc[met] = [metabolite.concentration_lb,metabolite.concentration_ub] 
    return(reaction_list,metabolite_list,lb_list,ub_list,coef_matrix,reaction_g0_list,reaction_kcat_MW_list,metabolites_lnc)

def Template_Concretemodel(model,obj_target=None,set_obj_value=False,set_substrate_ini=False,\
    set_biomass_ini=False,biomass_value=None,biomass_id=None,set_metabolite=False,set_Df=False,set_obj_B_value=False,set_stoi_matrix=False,\
    set_bound=False,set_enzyme_constraint=False,set_integer=False,set_metabolite_ratio=False,set_thermodynamics=False,B_value=None,\
    set_obj_E_value=False,set_obj_V_value=False,set_obj_TM_value=False,set_obj_Met_value=False,set_obj_single_E_value=False,E_total=None,\
    Bottleneck_reaction_list=None,set_Bottleneck_reaction=False,custom_constraints=None):
    
    """According to the parameter conditions provided by the user, the specific pyomo model is returned.

    Notes
    ----------
    * reaction_list: List of reaction IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolite_list: List of metabolite IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * coef_matrix: The model coefficient matrix.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolites_lnC:Upper and lower bound information of metabolite concentration (Natural logarithm).
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
    * obj_name: Name of object, such as set_obj_value, set_obj_single_E_value, set_obj_TM_value and set_obj_Met_value.    
    * K_value: the maximum value minus the minimum value of reaction thermodynamic driving force (1249, maxDFi-minDFi).
    * obj_target: Type of object function (maximize or minimize).
    * set_obj_value: Set the flux as the object function (True or False).
    * set_substrate_ini: Adding substrate amount constraints (True or False).
    * substrate_name: Substrate input reaction ID in the model (such as EX_glc__D_e_reverse).
    * substrate_value: Set the upper bound for substrate input reaction flux (10 mmol/h/gDW).
    * set_biomass_ini: Set the lower bound for biomass synthesis reaction flux (True or False).
    * biomass_value: The lower bound of biomass synthesis reaction flux.
    * biomass_id: Biomass synthesis reaction ID in the model (BIOMASS_Ec_iML1515_core_75p37M).
    * set_metabolite: Set the upper and lower bounds of metabolite concentration (True or False).
    * set_Df: Adding thermodynamic driving force expression for reactions (True or False).
    * set_obj_B_value: The object function is the maximizing thermodynamic driving force of a pathway (True or False)
    * set_stoi_matrix: Adding flux balance constraints (True or False).
    * set_bound: Set the upper and lower bounds of reaction flux (True or False).
    * set_enzyme_constraint: Adding enzymamic constraints (True or False).
    * set_integer: Adding binary variables constraints (True or False).
    * set_metabolite_ratio: Adding concentration ratio constraints for metabolites (True or False).
    * set_thermodynamics: Adding thermodynamic constraints (True or False).
    * B_value: The value of maximizing the minimum thermodynamic driving force (MDF).
    * set_obj_E_value: The object function is the minimum enzyme cost of a pathway (True or False).
    * set_obj_V_value: The object function is the pFBA of a pathway (True or False)
    * set_obj_TM_value: The object function is the thermodynamic driving force of a reaction (True or False).
    * set_obj_Met_value: The object function is the concentration of a metabolite (True or False).
    * set_obj_single_E_value: The object function is the enzyme cost of a reaction (True or False).     
    * E_total: Total amount constraint of enzymes (0.114).
    * Bottleneck_reaction_list: A list extracted from the result file automatically.
    * set_Bottleneck_reaction: Adding integer variable constraints for specific reaction (True or False).
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
    for index,row in reaction_kcat_MW.iterrows():
        if row["kcat_MW"] < 0.00002:
            reaction_kcat_MW.drop(index, axis=0, inplace=True)
    Concretemodel = ConcreteModel()
    Concretemodel.metabolite = pyo.Var(metabolite_list,  within=Reals)
    Concretemodel.Df = pyo.Var(reaction_list,  within=Reals)
    Concretemodel.B = pyo.Var()
    Concretemodel.reaction = pyo.Var(reaction_list,  within=NonNegativeReals)
    Concretemodel.z = pyo.Var(reaction_list,  within=pyo.Binary)
    obj_name = model.obj_name
    K_value = model.K_value
    substrate = model.medium
    E_total = model.E_total
    # non_exchange_reaction_list = []
    # for reaction in model["reactions"]:
    #     has_positive = False
    #     has_negative = False
    #     for value in reaction["metabolites"].values():
    #         if value > 0:
    #             has_positive = True
    #         elif value < 0:
    #             has_negative = True
    #     if has_positive and has_negative:
    #         non_exchange_reaction_list.append(reaction["id"])

    #Set upper and lower bounds of metabolite concentration
    if set_metabolite:
        def set_metabolite(m,i):
            return  inequality(metabolites_lnC.loc[i,'lnClb'], m.metabolite[i], metabolites_lnC.loc[i,'lnCub'])
        Concretemodel.set_metabolite= Constraint(metabolite_list,rule=set_metabolite)           

    #thermodynamic driving force expression for reactions
    if set_Df:
        def set_Df(m,j):
            return  m.Df[j]==-reaction_g0.loc[j,'g0']-2.579*sum(coef_matrix[i,j]*m.metabolite[i]  for i in metabolite_list if (i,j) in coef_matrix.keys())
        rg0=list(set(list(reaction_g0.index)).intersection(set(reaction_list))) 
        Concretemodel.set_Df = Constraint(rg0,rule=set_Df)
    
    #Set the maximum flux as the object function
    if set_obj_value:   
        def set_obj_value(m):
            return m.reaction[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_value, sense=minimize)
            
    #Set the value of maximizing the minimum thermodynamic driving force as the object function
    if set_obj_B_value:             
        def set_obj_B_value(m):
            return m.B  
        Concretemodel.obj = Objective(rule=set_obj_B_value, sense=maximize)  

    #Set the minimum enzyme cost of a pathway as the object function
    if set_obj_E_value:             
        def set_obj_E_value(m):
            return sum(m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index)
        Concretemodel.obj = Objective(rule=set_obj_E_value, sense=minimize)  

    #Minimizing the flux sum of pathway (pFBA)
    if set_obj_V_value:             
        def set_obj_V_value(m):
            return sum(m.reaction[j] for j in reaction_list)
        Concretemodel.obj = Objective(rule=set_obj_V_value, sense=minimize)  

    #To calculate the variability of enzyme usage of single reaction.
    if set_obj_single_E_value:             
        def set_obj_single_E_value(m):
            return m.reaction[obj_name]/(reaction_kcat_MW.loc[obj_name,'kcat_MW'])
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_single_E_value, sense=maximize) 
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_single_E_value, sense=minimize)

    #To calculate the variability of thermodynamic driving force (Dfi) of single reaction.
    if set_obj_TM_value:   
        def set_obj_TM_value(m):
            return m.Df[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_TM_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_TM_value, sense=minimize)

    #To calculate the concentration variability of metabolites.
    if set_obj_Met_value:   
        def set_obj_Met_value(m):
            return m.metabolite[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_Met_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_Met_value, sense=minimize)

    #Adding flux balance constraints （FBA）
    if set_stoi_matrix:
        def set_stoi_matrix(m,i):
            return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
        # print(metabolite_list)
        Concretemodel.set_stoi_matrix = Constraint( metabolite_list,rule=set_stoi_matrix)

    #Adding the upper and lower bound constraints of reaction flux
    if set_bound:
        def set_bound(m,j):
            return inequality(lb_list[j],m.reaction[j],ub_list[j])
        Concretemodel.set_bound = Constraint(reaction_list,rule=set_bound) 

    #Set the upper bound for substrate input reaction flux
    if set_substrate_ini:
        def set_substrate_ini(m): 
            for substrate_name,substrate_value in substrate.items():
                return m.reaction[substrate_name] <= substrate_value 
        Concretemodel.set_substrate_ini = Constraint(rule=set_substrate_ini)   
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_biomass_ini:
        def set_biomass_ini(m): 
            return m.reaction[biomass_id] >=biomass_value
        Concretemodel.set_biomass_ini = Constraint(rule=set_biomass_ini)  

    #Adding enzymamic constraints
    if set_enzyme_constraint:
        def set_enzyme_constraint(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index)<= E_total
        Concretemodel.set_enzyme_constraint = Constraint(rule=set_enzyme_constraint)

    #Adding thermodynamic MDF(B) object function
    if set_obj_B_value:
        def set_obj_B_value(m,j):
            return m.B<=(m.Df[j]+(1-m.z[j])*K_value)
        Concretemodel.set_obj_B_value = Constraint(reaction_list, rule=set_obj_B_value)

    #Adding thermodynamic constraints
    if set_thermodynamics:
        def set_thermodynamics(m,j):
            return (m.Df[j]+(1-m.z[j])*K_value)>= B_value
        Concretemodel.set_thermodynamics = Constraint(reaction_list, rule=set_thermodynamics)
        
    #Adding binary variables constraints
    if set_integer:
        def set_integer(m,j):
            return m.reaction[j]<=m.z[j]*ub_list[j] 
        Concretemodel.set_integer = Constraint(reaction_list,rule=set_integer)    

    #Adding concentration ratio constraints for metabolites
    if set_metabolite_ratio:
        def set_atp_adp(m):
            # print(m)
            # print(m.metabolite)
            # return m.metabolite['atp_c']-m.metabolite['adp_c']==np.log(10)
            try:
                # print(m.metabolite['cpd00002_c0'])
                return m.metabolite['cpd00002_c0']-m.metabolite['cpd00008_c0']==np.log(10)
            except:
                # print(m.metabolite.get_by_id('cpd00002_c0'))
                return m.metabolite.get_by_id('cpd00002_c0')-m.metabolite.get_by_id('cpd00008_c0')==np.log(10)
        def set_adp_amp(m):
            # return m.metabolite['adp_c']-m.metabolite['amp_c']==np.log(1)
            try:
                return m.metabolite['cpd00008_c0']-m.metabolite['cpd00018_c0']==np.log(1)
            except:
                return m.metabolite.get_by_id('cpd00008_c0')-m.metabolite.get_by_id('cpd00018_c0')==np.log(10)
        def set_nad_nadh(m):
            # return m.metabolite['nad_c']-m.metabolite['nadh_c']==np.log(10)
            try:
                return m.metabolite['cpd00003_c0']-m.metabolite['cpd00004_c0']==np.log(10)
            except:
                return m.metabolite.get_by_id('cpd00003_c0')-m.metabolite.get_by_id('cpd00004_c0')==np.log(10)
        def set_nadph_nadp(m):
            # return m.metabolite['nadph_c']-m.metabolite['nadp_c']==np.log(10)
            try:
                return m.metabolite['cpd00005_c0']-m.metabolite['cpd00006_c0']==np.log(10)
            except:
                return m.metabolite.get_by_id('cpd00005_c0')-m.metabolite.get_by_id('cpd00006_c0')==np.log(10)
        def set_hco3_co2(m):
            # return m.metabolite['hco3_c']-m.metabolite['co2_c']==np.log(2)
            try:
                return m.metabolite['cpd00242_c0']-m.metabolite['cpd00011_c0']==np.log(2)
            except:
                return m.metabolite.get_by_id('cpd00242_c0')-m.metabolite.get_by_id('cpd00011_c0')==np.log(10)

        Concretemodel.set_atp_adp = Constraint(rule=set_atp_adp) 
        Concretemodel.set_adp_amp = Constraint(rule=set_adp_amp) 
        Concretemodel.set_nad_nadh = Constraint(rule=set_nad_nadh) 
        Concretemodel.set_nadph_nadp = Constraint(rule=set_nadph_nadp) 
        Concretemodel.set_hco3_co2 = Constraint(rule=set_hco3_co2)

    #Adding Bottleneck reaction constraints
    if set_Bottleneck_reaction:
        def set_Bottleneck_reaction(m,j):
            return m.z[j]==1 
        Concretemodel.set_Bottleneck_reaction = Constraint(Bottleneck_reaction_list,rule=set_Bottleneck_reaction) 

    #Adding custom constraints (list of Constraint objects)
    if custom_constraints:
        if isinstance(custom_constraints, list):
            for con in custom_constraints:
                if isinstance(con, pyo.Constraint):
                    new_con = pyo.Constraint(con.index_set(), con.rule, con.bounds, con.domain)
                    Concretemodel.add_component(new_con.name, new_con)
                else:
                    Concretemodel.add_component('custom_con_%d' % i, pyo.Constraint(rule=con))
        elif isinstance(custom_constraints, dict):
            for i, (con_name, con) in enumerate(custom_constraints.items()):
                if callable(con):
                    Concretemodel.add_component(con_name, pyo.Constraint(rule=con))
                elif isinstance(con, pyo.Constraint):
                    new_con = pyo.Constraint(con.index_set(), con.rule, con.bounds, con.domain)
                    Concretemodel.add_component(new_con.name, new_con)
                else:
                    Concretemodel.add_component(con_name, pyo.Constraint(rule=lambda m: con))

    return Concretemodel

#Constructing a GEM (iML1515 model) using Pyomo Concretemodel framework
def calculate_GEM(model,obj_target,custom_constraints=None):
    """Constructing a GEM using Pyomo Concretemodel framework
    
    Arguments:
    ------------
    * model: a dictionary of model
    * obj_name: the name of objective reaction
    * obj_target: the target of objective reaction
    * substrate_name: the name of substrate input reaction
    * substrate_value: the value of substrate input reaction
    * solver: the name of solver
    Return:
    * opt: the optimal solution of the model
    """
    for substrate_name,substrate_value in model.medium.items():
        try:
            model.get_reaction_by_id(substrate_name).set_upper_bound(substrate_value)
        except:
            model.reactions.get_by_id(substrate_name).set_upper_bound(substrate_value)

    GEM=Template_Concretemodel(model,obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
            set_stoi_matrix=True,set_bound=True,custom_constraints=custom_constraints)
    opt = Model_Solve(GEM,model.solver)
    return opt

def Model_Solve(model,solver):
    """Solving programming problems
    Arguments:
    ------------
    * model: a dictionary of model
    * solver: the name of solver
    """
    opt = pyo.SolverFactory(solver)
    opt.solve(model)
    return model

def calculate_ECM(model,obj_target,custom_constraints=None):
    """Constructing a enzyme-constraint GEM using Pyomo Concretemodel framework
    
    Arguments:
    ------------
    * model: a dictionary of model
    * obj_name: the name of objective reaction
    * obj_target: the target of objective reaction
    * substrate_name: the name of substrate input reaction
    * substrate_value: the value of substrate input reaction
    * E_total: the total enzyme mass
    * solver: the name of solver
    * custom_constraints: custom constraints dict
    Return:
    * opt: the optimal solution of the model
    """
    for substrate_name,substrate_value in model.medium.items():
        try:
            model.get_reaction_by_id(substrate_name).set_upper_bound(substrate_value)
        except:
            model.reactions.get_by_id(substrate_name).set_upper_bound(substrate_value)

    ECM=Template_Concretemodel(model,obj_target=obj_target,set_obj_value=True,\
        set_substrate_ini=True,set_stoi_matrix=True,set_bound=True,\
        set_enzyme_constraint=True,custom_constraints=custom_constraints)
    opt = Model_Solve(ECM,model.solver)
    return opt


#Constructing a thermodynamic constraints model (EcoTCM) using Pyomo Concretemodel framework
def calculate_TCM(model,obj_target,B_value,custom_constraints=None):
    """Constructing a thermodynamics-constraint GEM using Pyomo Concretemodel framework
    
    Arguments:
    ------------
    * model: a dictionary of model
    * obj_name: the name of objective reaction
    * obj_target: the target of objective reaction
    * substrate_name: the name of substrate input reaction
    * substrate_value: the value of substrate input reaction
    * K_value: the value of MDF
    * B_value: the value of thermodynamic MDF(B)
    * solver: the name of solver
    Return:
    * opt: the optimal solution of the model
    """
    for substrate_name,substrate_value in model.medium.items():
        try:
            model.get_reaction_by_id(substrate_name).set_upper_bound(substrate_value)
        except:
            model.reactions.get_by_id(substrate_name).set_upper_bound(substrate_value)

    TCM=Template_Concretemodel(model,\
        obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
        set_stoi_matrix=True,set_bound=True,set_metabolite=True,\
        set_Df=True,set_integer=True,set_metabolite_ratio=True,set_thermodynamics=True,B_value=B_value,custom_constraints=custom_constraints)
    opt = Model_Solve(TCM,model.solver)
    return opt

#Constructing a enzymatic and thermodynamic constraints model (EcoETM) using Pyomo Concretemodel framework
def calculate_ETM(model,obj_target,B_value,custom_constraints=None):
    """Constructing a thermodynamics-constraint and enzyme-constraint GEM using Pyomo Concretemodel framework
    
    Arguments:
    ------------
    * model: a dictionary of model
    * obj_name: the name of objective reaction
    * obj_target: the target of objective reaction
    * substrate_name: the name of substrate input reaction
    * substrate_value: the value of substrate input reaction
    * E_total: the total enzyme mass
    * K_value: the value of MDF
    * B_value: the value of thermodynamic MDF(B)
    * solver: the name of solver
    * custom_constraints: custom constraints dict
    Return:
    * opt: the optimal solution of the model
    """
    for substrate_name,substrate_value in model.medium.items():
        try:
            model.get_reaction_by_id(substrate_name).set_upper_bound(substrate_value)
        except:
            model.reactions.get_by_id(substrate_name).set_upper_bound(substrate_value)
    
    ETM=Template_Concretemodel(model,\
        obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
        set_stoi_matrix=True,set_bound=True,set_metabolite=True,\
        set_Df=True,set_integer=True,set_metabolite_ratio=True,set_thermodynamics=True,B_value=B_value,\
        set_enzyme_constraint=True,custom_constraints=custom_constraints)
    opt = Model_Solve(ETM,model.solver)
    return opt

#Minimum flux sum calculation（pFBA）
def Min_Flux_Sum_Calculation(model,model_type,biomass_value,biomass_id,B_value):
    """Calculation of minimum flux sum for a model.
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
    * min_V: The value of minimum flux sum.
    """
    if model_type == "GEM":
        Concretemodel = Template_Concretemodel(model,
            set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,
            set_stoi_matrix=True,set_bound=True,
            set_integer=True,set_obj_V_value=True)
    elif model_type == "ECM":
        Concretemodel = Template_Concretemodel(model,
            set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,
            set_stoi_matrix=True,set_bound=True,
            set_enzyme_constraint=True,set_integer=True,set_obj_V_value=True)
    elif model_type == "TCM":
        Concretemodel = Template_Concretemodel(model,
            set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,
            set_metabolite=True,set_Df=True,set_metabolite_ratio=True,set_stoi_matrix=True,
            set_bound=True,set_thermodynamics=True,B_value=B_value,set_integer=True,set_obj_V_value=True)
    elif model_type == "ETM":
        Concretemodel = Template_Concretemodel(model=model,\
            set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,\
            set_metabolite=True,set_Df=True,set_stoi_matrix=True,\
            set_bound=True,set_enzyme_constraint=True,set_integer=True,set_metabolite_ratio=True,\
            set_thermodynamics=True,B_value=B_value,set_obj_V_value=True)
    else:
        raise ValueError(
                f"{model_type} is not a model type"
        )
    opt=Model_Solve(Concretemodel,model.solver)
    min_V=opt.obj()
    return [min_V,Concretemodel]

def Max_Growth_Rate_Calculation(model,obj_target,B_value):
    """Calculation of maximum growth rate for a model.
    Arguments:
    ----------
    * model: The model to be calculated.
    * obj_name: Object reaction ID.
    * obj_target: Type of object function (maximize or minimize).
    * substrate_name: Substrate name.
    * substrate_value: Substrate value.
    * K_value: The value of MDF.
    * E_total: The total enzyme mass.
    * B_value: The value of MDF(B).
    * solver: The solver used to solve the model.
    Return:
    ----------
    * opt.obj(): The value of maximum growth rate.
    """
    Concretemodel=Template_Concretemodel(model=model,\
        obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
        set_metabolite=True,set_Df=True,set_stoi_matrix=True,\
        set_bound=True,set_enzyme_constraint=True,set_integer=True,set_metabolite_ratio=True,\
        set_thermodynamics=True,B_value=B_value)
    opt=Model_Solve(Concretemodel,model.solver)
    return opt.obj()