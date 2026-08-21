import cobra
from cobra.util.solver import set_objective
from typing import Dict,Any
import json
from cobra.core import Reaction,Model
import re
import pandas as pd

def json_load(path: str) -> Dict[Any, Any]:
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder
def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)

def trans_model2standard_json_etgem(model_file, output_path,new_model_file):
    """Converts the given model file to a standard ETGEMs JSON file.
    Arguments
    ----------
    * model_file:   The path of the JSON file of GEMs that shall be converted   
    """
    if model_file.split(".")[-1] == 'json':
        model = cobra.io.load_json_model(model_file)
    elif model_file.split(".")[-1] == 'xml':
        model = cobra.io.read_sbml_model(model_file)
    remove_unused_metabolites(model)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = output_path + "%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    cobra_model = cobra.io.load_json_model(json_path)
    dictionary_model = json_load(json_path)
    dictionary_model["additional_properties"] = {}
    for eachreaction in dictionary_model['reactions']:
        eachreaction["g0"] = -99999
        eachreaction["kcat"] = 99999
        eachreaction["MW"] = 0.001
        if "objective_coefficient" in eachreaction.keys():
            dictionary_model["additional_properties"]["obj_name"] = eachreaction["id"]
        for key,value in eachreaction["metabolites"].items():
            formatted_value = "{:.5g}".format(value)
            eachreaction["metabolites"][key] = float(formatted_value)
    for eachmetabolite in dictionary_model['metabolites']:
        # if eachmetabolite["id"][-1] != 'c':
        #     eachmetabolite["concentration_ub"] = 0.0
        #     eachmetabolite["concentration_lb"] = 0.0
        # elif eachmetabolite["id"][-1] == 'c':
        eachmetabolite["concentration_ub"] = -3.912023005
        eachmetabolite["concentration_lb"] = -14.50865774
    for eachgene in dictionary_model['genes']:
        cobra_gene = cobra_model.genes.get_by_id(eachgene["id"])
        if "uniprot" in cobra_gene.annotation.keys():
            eachgene["uniprot"] = cobra_gene.annotation["uniprot"]
        else:
            eachgene["uniprot"] = ""
    dictionary_model["additional_properties"]["medium"] = {}
    dictionary_model["additional_properties"]["K_value"] = 0.0
    dictionary_model["additional_properties"]["E_total"] = 1.0
    dictionary_model["additional_properties"]["solver"] = "cplex_direct"
    dictionary_model["environment"] = {}
    dictionary_model["environment"]["p_h"] = 7.5
    dictionary_model["environment"]["p_mg"] = 3.0
    dictionary_model["environment"]["ionic_strength"] = "0.1M"
    dictionary_model["environment"]["temperature"] = "310.15K"

    model_output_path = output_path+new_model_file
    json_write(model_output_path,dictionary_model)
    return dictionary_model

def remove_unused_metabolites(model):
    """Remove unused metabolites from model
    Arguments
    ----------
    * model: cobra.Model.
    """
    metabolite_used_list = []
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            metabolite_used_list.append(metabolite.id)
            metabolite_used_list = list(set(metabolite_used_list))
    metabolite_unused_num = 1
    while(metabolite_unused_num > 0):
        metabolite_unused_num = 0
        for metabolite in model.metabolites:
            if metabolite.id not in metabolite_used_list:
                metabolite_unused_num += 1
                model.remove_metabolites([metabolite])

def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value.replace("(","").replace(")","")
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value.replace("(","").replace(")","")
                    model.add_reactions([r_add])
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model

def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    """
    #warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction.metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            reverse_reaction.name = str(reaction.name)+' reverse'
            reverse_reaction.annotation = reaction.annotation
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reaction.gene_reaction_rule = reaction.gene_reaction_rule.replace("(","").replace(")","")
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)

def result2flux_solution_df(model,opt_model):
    reactions_list = []
    fluxes = []
    for reaction in model.reactions.keys():
        reactions_list.append(reaction)
        fluxes.append(opt_model.reaction[reaction].value)
    final_data = {'reaction' : reactions_list,
                'fluxes' : fluxes}
    flux_solution_df = pd.DataFrame(final_data)
    flux_solution_df = flux_solution_df.set_index('reaction')
    flux_solution_df_select=flux_solution_df[abs(flux_solution_df['fluxes'])>0.0000001]#只储存大于0的值e-6
    flux_solution_df_sorted = flux_solution_df_select.sort_values(by='fluxes', ascending=False)

    flux_solution_return=pd.DataFrame() 
    for index, row in flux_solution_df_sorted.iterrows():
        flux_solution_return.loc[index,'reaction_id']=index
        flux_solution_return.loc[index,'reaction_name']=model.get_reaction_by_id(index).name
        # flux_solution_return.loc[index,'fluxes']=abs(round(row['fluxes'],3))
        flux_solution_return.loc[index,'fluxes']=round(row['fluxes'],7)
        flux_solution_return.loc[index,'gene id']=model.get_reaction_by_id(index).gene_reaction_rule
        #flux_solution_return.loc[index,'ec-code']=model.reactions.get_by_id(index).annotation['ec-code'] 
        flux_solution_return.loc[index,'equ']=model.get_reaction_by_id(index).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
        flux_solution_return.loc[index,'abs_fluxes']=abs(round(row['fluxes'],7))
    
    return flux_solution_return
