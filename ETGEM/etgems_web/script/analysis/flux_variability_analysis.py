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
from .flux_balance_analysis import *
from concurrent.futures import ProcessPoolExecutor, as_completed

def flux_variability_analysis_single_reaction(model, model_type, reaction, obj_target, obj_value, fraction_of_optimum, B_value,biomass_id):
    """flux variablity analysis of single reaction, the reaction is set to be the objective function
    , and the biomass is set to be the constraint.
    Arguments:
    ----------
    * model: a dictionary of model.
    * reaction: The reaction that is set to be the objective function.
    * obj_target: The target of the objective function.
    * substrate_name: The name of the substrate.
    * substrate_value: The value of the substrate.
    * biomass_id: The id of the biomass.
    * biomass_value: The value of the biomass.
    * fraction_of_optimum: The fraction of the optimum.
    * K_value: The value of MDF.
    * B_value: The value of MDF(B).
    * E_total: The total enzyme.
    * solver: The solver used to solve the model.
    Returns:
    ----------
    * result: The flux variablity analysis of single reaction.
    """
    tem_obj_name = model.obj_name
    model.set_obj_name(reaction)
    result = pd.DataFrame()
    if model_type == "GEM":
        Concretemodel = Template_Concretemodel(model,obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
            set_stoi_matrix=True,set_bound=True)
    elif model_type == "ECM":
        Concretemodel = Template_Concretemodel(model,obj_target=obj_target,set_obj_value=True,\
            set_substrate_ini=True,set_stoi_matrix=True,set_bound=True,\
            set_enzyme_constraint=True)
    elif model_type == "TCM":
        Concretemodel = Template_Concretemodel(model,\
            obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
            set_stoi_matrix=True,set_bound=True,set_metabolite=True,\
            set_Df=True,set_integer=True,set_metabolite_ratio=True,set_thermodynamics=True,B_value=B_value)
    elif model_type == "ETM":
        Concretemodel = Template_Concretemodel(model,\
            obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,\
            set_stoi_matrix=True,set_bound=True,set_metabolite=True,\
            set_Df=True,set_integer=True,set_metabolite_ratio=True,set_thermodynamics=True,B_value=B_value,\
            set_enzyme_constraint=True)
    else:
        raise ValueError(
                f"{model_type} is not a model type"
        )

    def set_biomass_value(m, obj_id, obj_value, fraction_of_optimum):
        return inequality(fraction_of_optimum*obj_value, m.reaction[obj_id], obj_value)
    
    # Concretemodel.set_biomass_value = Constraint(rule=set_biomass_value)
    Concretemodel.set_biomass_value = Constraint(expr=set_biomass_value(Concretemodel, biomass_id, obj_value, fraction_of_optimum))

    opt=Model_Solve(Concretemodel,model.solver)
    result.loc[reaction,obj_target]=opt.obj()
    model.set_obj_name(tem_obj_name)
    return result

def flux_variability_analysis(model, model_type, obj_value, fraction_of_optimum, B_value, biomass_id, max_workers):
    """flux variablity analysis of total reactions in model, each reaction is set to be the objective function
    , and the biomass is set to be the constraint.
    Arguments:
    ----------
    * model: a dictionary of model.
    * obj_target: The target of the objective function.
    * substrate_name: The name of the substrate.
    * substrate_value: The value of the substrate.
    * biomass_id: The id of the biomass.
    * biomass_value: The value of the biomass.
    * fraction_of_optimum: The fraction of the optimum.
    * K_value: The value of MDF.
    * B_value: The value of MDF(B).
    * E_total: The total enzyme.
    * solver: The solver used to solve the model.
    * max_workers: The maximum number of workers.
    Returns:
    ----------
    * result: The flux variablity analysis of total reactions in model.
    """
    reaction_list = []
    flux_variability_analysis = pd.DataFrame()
    for reaction in model.reactions.values():
        reaction_list.append(reaction.id)
    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(flux_variability_analysis_single_reaction,model, model_type, eachreaction, "maximize", 
                              obj_value, fraction_of_optimum, B_value, biomass_id): eachreaction for eachreaction in reaction_list}
        for future in as_completed(futures):
            tmp = future.result()
            for eachindex in tmp.index:
                flux_variability_analysis.loc[eachindex,'maximize']=tmp.loc[eachindex,'maximize']
                
    with ProcessPoolExecutor(max_workers) as executor:
        futures = {executor.submit(flux_variability_analysis_single_reaction,model, model_type, eachreaction, "minimize", 
                              obj_value, fraction_of_optimum, B_value, biomass_id): eachreaction for eachreaction in reaction_list}
        for future in as_completed(futures):
            tmp = future.result()
            for eachindex in tmp.index:
                flux_variability_analysis.loc[eachindex,'minimize']=tmp.loc[eachindex,'minimize']
    return flux_variability_analysis