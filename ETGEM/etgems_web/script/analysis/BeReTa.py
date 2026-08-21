import numpy as np
from scipy.stats import pearsonr,linregress
import cobra
import  re
import json
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.model import Model
from ..core.environment import Environment
from ..core.gene import Gene
import cobra
from typing import Dict,Any
from cobra.util.solver import set_objective
import pandas as pd
import matplotlib.pyplot as plt
import pyomo.environ as pyo
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np
from .flux_balance_analysis import *
from concurrent.futures import ProcessPoolExecutor, as_completed

def BeReTa_use(model, model_type, reg_net, exp_comp, biomass_id, product_rxn, pval_cut, n_target_cut, f_target_cut):
    # INPUTS
    # model: COBRA model structure
    # reg_net: Transcriptional regulatory network structure
    #   regulator: Cell containing regulator names
    #   target: Cell containing target names
    #   sign: Vector containing sign of regulation (activator, 1; repressor, -1; unknown, 0) (optional)
    #   regulator, target, and sign should have same length
    # exp_comp: Gene expression compendium data structure
    #   genes: Cell containing gene names
    #   expression: Gene expression levels
    #   genes and expression should have same length
    # product_rxn: Exchange reaction for target product
    #
    # OUTPUTS
    # BeReTaSolution: BeReTa solution structure (target regulator, beneficial score, p-value)

    # Generate regulatory strength matrix (RegStr) and flux slope vector (q_slope)
    tem_obj = model.obj_name
    #product_rxn = "DM_" + product
    #if product_rxn in model.reactions.keys():
    #    model.remove_reactions(product_rxn)
    #model.add_demand_reaction(product)
    model.obj_name = biomass_id
    RegStr, RegStr_gene = generateRegStrMatrix(reg_net, exp_comp, model)
    #print(len(RegStr['Reg_TF']))
    q_slope = calculateFluxSlope(model, model_type, product_rxn)##这步求解有问题
    # Trimming
    reg_list = RegStr['Reg_TF'] 
    Reg_RS = RegStr['Reg_RS']
    
    #print(RegStr,q_slope,RegStr.shape,q_slope.shape)
    # Calculate beneficial scores for all transcriptional regulators (TRs)
    beneficial_score = np.dot(Reg_RS, q_slope)
    #print(type(beneficial_score),beneficial_score, q_slope)
    # Perform permutation test
    beneficial_score_pval = permutationTestBeReTa(model, Reg_RS, q_slope, beneficial_score, reg_list)
    # Apply target criteria to select BeReTa targets
    BeReTaSolution = selectBeReTaTargets(model, RegStr, RegStr_gene, q_slope, beneficial_score, beneficial_score_pval, reg_list, pval_cut, n_target_cut, f_target_cut)
    model.set_obj_name(tem_obj)
    return BeReTaSolution

def BeReTa_v1(model, model_type, reg_net, exp_comp, biomass_id, product, pval_cut, n_target_cut, f_target_cut):
    # INPUTS
    # model: COBRA model structure
    # reg_net: Transcriptional regulatory network structure
    #   regulator: Cell containing regulator names
    #   target: Cell containing target names
    #   sign: Vector containing sign of regulation (activator, 1; repressor, -1; unknown, 0) (optional)
    #   regulator, target, and sign should have same length
    # exp_comp: Gene expression compendium data structure
    #   genes: Cell containing gene names
    #   expression: Gene expression levels
    #   genes and expression should have same length
    # product_rxn: Exchange reaction for target product
    #
    # OUTPUTS
    # BeReTaSolution: BeReTa solution structure (target regulator, beneficial score, p-value)

    # Generate regulatory strength matrix (RegStr) and flux slope vector (q_slope)
    tem_obj = model.obj_name
    product_rxn = "DM_" + product
    if product_rxn in model.reactions.keys():
        model.remove_reactions(product_rxn)
    model.add_demand_reaction(product)
    model.obj_name = biomass_id
    RegStr, RegStr_gene = generateRegStrMatrix(reg_net, exp_comp, model)
    #print(len(RegStr['Reg_TF']))
    q_slope = calculateFluxSlope(model, model_type, product_rxn)##这步求解有问题
    # Trimming
    reg_list = RegStr['Reg_TF'] 
    Reg_RS = RegStr['Reg_RS']
    
    #print(RegStr,q_slope,RegStr.shape,q_slope.shape)
    # Calculate beneficial scores for all transcriptional regulators (TRs)
    beneficial_score = np.dot(Reg_RS, q_slope)
    #print(type(beneficial_score),beneficial_score, q_slope)
    # Perform permutation test
    beneficial_score_pval = permutationTestBeReTa(model, Reg_RS, q_slope, beneficial_score, reg_list)
    # Apply target criteria to select BeReTa targets
    BeReTaSolution = selectBeReTaTargets(model, RegStr, RegStr_gene, q_slope, beneficial_score, beneficial_score_pval, reg_list, pval_cut, n_target_cut, f_target_cut)
    model.set_obj_name(tem_obj)
    return BeReTaSolution

def generateRegStrMatrix(reg_net, exp_comp, model):
    # Generate regulatory strength matrix (RS)
    regulator = reg_net['regulator']
    target = reg_net['target']
    expression = exp_comp['expression']
    expressionid = exp_comp['genes']
    model_genes = [gene.id for gene in model.genes.values()]
    model_rxns = [reaction.id for reaction in model.reactions.values()]
    corr_mat1_col = list(model_genes)

    if reg_net.shape[1] == 3:
        sign = reg_net['sign']
    else:
        sign = np.zeros((reg_net.shape[0], reg_net.shape[1]))
    # Discard regulatory interactions for which correlation cannot be computed.    
    regulator_temp = []
    target_temp = []
    sign_temp = []
    for i in range(len(regulator)):
        regulator_id = np.where(expressionid == regulator[i])[0]
        target_id = np.where(expressionid == target[i])[0]
        if (len(regulator_id) == 1) and (len(target_id) == 1) and target[i] in corr_mat1_col:#有表达值，并且在模型中
            regulator_temp.append(regulator[i])
            target_temp.append(target[i])
            sign_temp.append(sign[i])
    regulator = regulator_temp
    target = target_temp
    sign = sign_temp

    reg_list = np.unique(regulator)
    corr_mat1 = np.zeros((len(reg_list), len(model_genes)))
    corr_mat1_row = list(reg_list)

    # Generate regulatory strength (RS) matrix for TR-gene interactions
    for i in range(len(regulator)):
        row_id = corr_mat1_row.index(regulator[i])
        col_id = corr_mat1_col.index(target[i])
        #print(row_id,col_id)
        x1 = expression[expressionid==regulator[i]].values[0]
        x2 = expression[expressionid==target[i]].values[0]
        ## Mr.Mao
        corr_val = pearsonr(np.array(x1).flatten(), np.array(x2).flatten())[0]
        #print(corr_val)
        if sign[i] == 0:
            corr_mat1[row_id, col_id] = corr_val
        elif sign[i] == 1:
            #print(row_id,col_id,corr_val)
            corr_mat1[row_id, col_id] = abs(corr_val)
        elif sign[i] == -1:
            corr_mat1[row_id, col_id] = -abs(corr_val)
    #print(corr_mat1.shape)
    #Remove zero columns in the matrix
    row_index = np.where(np.sum(np.abs(corr_mat1), axis=1) != 0)[0]
    col_index = np.where(np.sum(np.abs(corr_mat1), axis=0) != 0)[0]
    #print(row_index,col_index)
    corr_mat2 = corr_mat1[row_index, :][:, col_index]
    corr_mat2_row = [corr_mat1_row[i] for i in row_index]
    corr_mat2_col = [corr_mat1_col[i] for i in col_index]
    RegStr_gene = {}
    RegStr_gene['RSmat_gene'] = corr_mat2_col
    RegStr_gene['RSmat_TF'] = corr_mat2_row
    RegStr_gene['RSmat_mat'] = corr_mat2
    #RegStr_gene = [['RSmat_gene'] + list(corr_mat2_col)] + list(zip(corr_mat2_row, corr_mat2.tolist()))

    # Generate regulatory strength (RS) matrix for TR-rxn interactions
    # GPR mapping using average function
    corr_mat3 = np.zeros((len(corr_mat2_row), len(model_rxns)))
    corr_mat3_row = corr_mat2_row
    corr_mat3_col = np.array(model_rxns)

    for i in range(corr_mat2.shape[0]):
        gpr_genes = np.array(corr_mat2_col)
        gpr_levels = np.array(corr_mat2[i, :]).reshape(-1, 1)
        gpr_genes2 = gpr_genes[np.nonzero(gpr_levels)[0]]
        gpr_levels2 = gpr_levels[np.nonzero(gpr_levels)[0]]
        levels = gene_to_reaction_levels(model, list(gpr_genes2), gpr_levels2, lambda x, y: (x + y) / 2, lambda x, y: (x + y) / 2)
        corr_mat3[i, np.nonzero(~np.isnan(levels))] = levels[np.nonzero(~np.isnan(levels))]
    
    RegStr = {}
    RegStr['Reg_rxn'] = corr_mat3_col
    RegStr['Reg_TF'] = corr_mat3_row
    RegStr['Reg_RS'] = corr_mat3

    return [RegStr, RegStr_gene]

def gene_to_reaction_levels(model, genes, levels, f_and, f_or):
    # Original authors: Daniel Machado and Markus Herrgard, PLoS Computational Biology, 2014
    # Code obtained from https://github.com/cdanielmachado/transcript2flux
    #
    # Convert gene expression levels to reaction levels using GPR associations.
    # Level is NaN if there is no GPR for the reaction or no measured genes.
    #
    # INPUTS
    #       model - cobra model
    #       genes - gene names
    #       levels - gene expression levels
    #       f_and - function to replace AND
    #       f_or - function to replace OR
    #
    # OUTPUTS
    #       reaction_levels - reaction expression levels
    #
    # Author: Daniel Machado, 2013
    reaction_levels = np.zeros(len(model.reactions.values()))
    i = 0 
    for reaction in model.reactions.values():
        #print(reaction.gene_reaction_rule)
        [result, level] = eval_gpr(reaction.gene_reaction_rule, genes, levels, f_and, f_or)
        #print(result)
        reaction_levels[i] = result
        i=i+1
    return reaction_levels

def eval_gpr(rule, genes, levels, f_and, f_or):
    # Original authors: Daniel Machado and Markus Herrgard, PLoS Computational Biology, 2014
    # Code obtained from https://github.com/cdanielmachado/transcript2flux
    #
    # Evaluate the expression level for a single reaction using the GPRs.
    # Note: Computes the expression level even if there are missing measured
    # values for the given rule. This implementation is a modified version of
    # an implementation provided in [Lee et al, BMC Sys Biol, 2012]
    EVAL_OK = 1
    PARTIAL_MEASUREMENTS = 0
    NO_GPR_ERROR = -1
    NO_MEASUREMENTS = -2
    MAX_EVALS_EXCEEDED = -3
    MAX_EVALS = 1000
    NONETYPE = 'NaN'
    NUMBER = r'[0-9\.\-e]+'
    MAYBE_NUMBER = f"{NUMBER}|{NONETYPE}"
    expression = rule
    result = np.nan
    status = EVAL_OK
    if not expression:
        status = NO_GPR_ERROR
    else:
        rule_genes = list(set(re.findall(r'\b(\w+|-)+\b', expression)) - {'and', 'or'})
        total_measured = 0
        for gene in rule_genes:
            #print(gene,rule_genes,genes)
            if gene in genes:
                j = genes.index(gene)
                level = str(levels[j])
                total_measured += 1
            else:
                level = NONETYPE
            expression = re.sub(r'\b{}\b'.format(gene), level, expression)
        if total_measured == 0:
            status = NO_MEASUREMENTS
        else:
            if total_measured < len(rule_genes):
                status = PARTIAL_MEASUREMENTS
            maybe_and = lambda a, b: maybe_functor(f_and, a, b)
            maybe_or = lambda a, b: maybe_functor(f_or, a, b)
            str_wrapper = lambda f, a, b: str(f(float(a), float(b)))
            counter = 0
            while np.isnan(result):
                counter += 1
                if counter > MAX_EVALS:
                    status = MAX_EVALS_EXCEEDED
                    break
                try:
                    result = eval(expression)[0]
                except Exception as e:
                    paren_expr = r'\(\s*({})\s*\)'.format(MAYBE_NUMBER)
                    and_expr = r'({})\s+and\s+({})'.format(MAYBE_NUMBER, MAYBE_NUMBER)
                    or_expr = r'({})\s+or\s+({})'.format(MAYBE_NUMBER, MAYBE_NUMBER)
                    expression = re.sub(paren_expr, r'\1', expression)
                    expression = re.sub(and_expr, lambda m: str_wrapper(maybe_and, m.group(1), m.group(2)), expression)
                    expression = re.sub(or_expr, lambda m: str_wrapper(maybe_or, m.group(1), m.group(2)), expression)
    return result, status

def maybe_functor(f, a, b):
    # Original authors: Daniel Machado and Markus Herrgard, PLoS Computational Biology, 2014
    # Code obtained from https://github.com/cdanielmachado/transcript2flux
    if np.isnan(a) and np.isnan(b):
        c = np.nan
    elif not np.isnan(a) and np.isnan(b):
        c = a
    elif np.isnan(a) and not np.isnan(b):
        c = b
    else:
        c = f(a, b)
    return c

def calculate(model, model_type, obj_target, B_value):
    if model_type == "GEM":
        initial_sol = calculate_GEM(model,obj_target)
        final_sol = Min_Flux_Sum_Calculation_for_BeReTa(model,model_type,initial_sol.reaction[model.obj_name].value,model.obj_name,B_value)
    elif model_type == "ecGEM":
        initial_sol = calculate_ECM(model,obj_target)
        final_sol = Min_Flux_Sum_Calculation_for_BeReTa(model,model_type,initial_sol.reaction[model.obj_name].value,model.obj_name,B_value)
    elif model_type == "tcGEM":
        initial_sol = calculate_TCM(model,obj_target,B_value)
        final_sol = Min_Flux_Sum_Calculation_for_BeReTa(model,model_type,initial_sol.reaction[model.obj_name].value,model.obj_name,B_value)
    elif model_type == "etGEM":
        initial_sol = calculate_ETM(model,obj_target,B_value)
        final_sol = Min_Flux_Sum_Calculation_for_BeReTa(model,model_type,initial_sol.reaction[model.obj_name].value,model.obj_name,B_value)
    return final_sol

#Minimum flux sum calculation（pFBA）
def Min_Flux_Sum_Calculation_for_BeReTa(model,model_type,biomass_value,biomass_id,B_value):
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
    elif model_type == "ecGEM":
        Concretemodel = Template_Concretemodel(model,
            set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,
            set_stoi_matrix=True,set_bound=True,
            set_enzyme_constraint=True,set_integer=True,set_obj_V_value=True)
    elif model_type == "tcGEM":
        Concretemodel = Template_Concretemodel(model,
            set_biomass_ini=True,biomass_value=biomass_value,biomass_id=biomass_id,set_substrate_ini=True,
            set_metabolite=True,set_Df=True,set_metabolite_ratio=True,set_stoi_matrix=True,
            set_bound=True,set_thermodynamics=True,B_value=B_value,set_integer=True,set_obj_V_value=True)
    elif model_type == "etGEM":
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
    return opt

def calculateFluxSlope(model, model_type, product_rxn):
    ## Calculate flux slope vector (q_slope)
    # Calculate initial and maximum flux of product Rxn.
    reaction_ids = [reaction.id for reaction in model.reactions.values()]
    product_rxn_id =  reaction_ids.index(product_rxn)
    initial_sol = calculate(model, model_type, "maximize", 0.0)
    initial_product_flux = initial_sol.reaction[product_rxn].value
    model.set_obj_name(product_rxn)
    maximum_sol = calculate(model, model_type,"maximize", 0.0)
    maximum_product_flux = maximum_sol.reaction[product_rxn].value
    # Flux scanning while increasing the product flux
    num_steps = 20
    flux_values = np.zeros((len(model.reactions.values()), num_steps))
    for i in range(num_steps):
        product_flux = initial_product_flux + ((maximum_product_flux - initial_product_flux) * ((i) / num_steps))
        model.get_reaction_by_id(product_rxn).set_upper_bound(product_flux)
        model.get_reaction_by_id(product_rxn).set_lower_bound(product_flux)
        enforced_sol = calculate(model, model_type,"maximize", 0.0)
        #enforced_sol = optimizeCbModel(model_enforced, 'max', 'one')
        for j in range(len(model.reactions.values())):
            flux_values[j, i] = np.abs(enforced_sol.reaction[list(model.reactions.keys())[j]].value)
    # Use linear regression to estimate flux slopes
    flux_values_product = flux_values[product_rxn_id, :]
    q_slope = np.zeros((len(model.reactions.values()), 1))
    for i in range(len(model.reactions.values())):
        linear_regression = linregress(flux_values_product, flux_values[i, :])
        q_slope[i, 0] = linear_regression.slope
    q_slope[q_slope < 0] = 0
    return q_slope

def exchanges(model):
    """Exchange reactions in model.
    Reactions that exchange mass with the exterior. Uses annotations
    and heuristics to exclude non-exchanges such as sink reactions.
    """
    exchange_list=[]
    for reaction in model.reactions.values():
        metabolites = reaction.metabolites
        if all(value > 0 for value in metabolites.values()) or all(value > 0 for value in metabolites.values()):
            exchange_list.append(reaction.id)
    return exchange_list

def findorphanreaction(model):
    orphan_reactions = []
    for reaction in model.reactions.values():
        if not reaction.gene_reaction_rule:
            orphan_reactions.append(reaction.id)
    return orphan_reactions

def permutationTestBeReTa(model, RegStr, q_slope, beneficial_score, reg_list):
    ## Calculate p-values for beneficial scores
    # Permute only for gene-associated rxns
    selExc = exchanges(model)
    orphans = findorphanreaction(model)
    reaction_ids = [reaction.id for reaction in model.reactions.values()]
    gene_associated = np.setdiff1d(reaction_ids, selExc)
    gene_associated = np.setdiff1d(gene_associated, orphans)
    rand_rxn_ids = np.where(np.isin(reaction_ids, gene_associated))[0]
    # Calculate permuted beneficial scores
    rand_num = 10000
    beneficial_score_rand = np.zeros((len(reg_list), rand_num))
    for i in range(rand_num):
        rand_order = np.argsort(np.random.rand(len(rand_rxn_ids)))
        q_slope_rand = q_slope.copy()
        q_slope_rand[rand_rxn_ids[rand_order]] = q_slope[rand_rxn_ids]
        #print(RegStr.shape,q_slope_rand.shape,np.dot(RegStr, q_slope_rand))
        tmp=np.dot(RegStr, q_slope_rand)
        j=0
        for value in np.nditer(tmp):
            beneficial_score_rand[j, i] = value
            j=j+1

    # Calculate p-values
    beneficial_score_pval = np.ones(len(reg_list))
    for i in range(len(reg_list)):
        #print(beneficial_score[i])
        if beneficial_score[i][0] > 0:
            beneficial_score_pval[i] = len(np.where(beneficial_score_rand[i, :] > beneficial_score[i][0])[0]) / rand_num
        elif beneficial_score[i][0] < 0:
            beneficial_score_pval[i] = len(np.where(beneficial_score_rand[i, :] < beneficial_score[i][0])[0]) / rand_num
    return beneficial_score_pval

def findGenesFromRxns(model, effective_rxns):
    gene_ids_list=[]
    for reaction in model.reactions.values():
        if reaction.id in effective_rxns:
            gene_list = reaction.gene_reaction_rule.split(" and ")
            gene_ids_list.append(gene_list)
    return gene_ids_list

def selectBeReTaTargets(model, RegStr, RegStr_gene, q_slope, beneficial_score, beneficial_score_pval, reg_list, pval_cut, n_target_cut, f_target_cut):
    ## Apply target criteria to select BeReTa targets
    # Define target criteria   (1) TR should have non-zero beneficial score.
    #pval_cut = 0.05 #(2) The p-value of the beneficial score should be less than 0.05.
    #n_target_cut = 1 #(3) TR should have two or more effective gene/reaction targets.
    #f_target_cut = 0.1 #(4) At least 10% of target metabolic genes of TR should be beneficial, i.e. have positive flux slopes.
    # Calculate the metrics
    n_effective_genes = np.zeros(len(reg_list))
    n_effective_rxns = np.zeros(len(reg_list))
    f_effective_genes = np.zeros(len(reg_list))
    Reg_rxn = RegStr['Reg_rxn']
    #print(Reg_rxn)
    Reg_RS = RegStr['Reg_RS']
    #print(type(q_slope),type(Reg_RS[1,:]))
    for i in range(len(reg_list)):
        regulated_genes = RegStr_gene['RSmat_gene']
        regulated_genes_num = np.count_nonzero(RegStr['Reg_RS'][i])
        j=0
        Reg_RStmp=np.zeros((Reg_RS.shape[1],1))
        #print(Reg_RS[i,:].shape)
        for value in np.nditer(Reg_RS[i,:]):
            Reg_RStmp[j, 0] = value
            j=j+1
        tmp = q_slope * Reg_RStmp
        #print(q_slope.shape,Reg_RStmp.shape,tmp.shape)
        # 创建布尔类型的掩码数组，判断元素是否为0
        mask = (tmp != 0) 
        #print(mask)
        # 使用 np.where() 获取满足条件的位置索引
        indices = np.where(mask)[0]
        #print(indices)
        #print(q_slope.shape,Reg_RS.shape)
        effective_rxns = Reg_rxn[indices]
        n_effective_rxns[i] = len(effective_rxns)
        if len(effective_rxns) > 0:
            effective_genes_temp = findGenesFromRxns(model, effective_rxns)
            effective_genes = []
            for j in range(len(effective_genes_temp)):
                effective_genes.extend(effective_genes_temp[j])
            effective_genes = np.unique(effective_genes)
            effective_genes = np.intersect1d(effective_genes, regulated_genes)
            n_effective_genes[i] = len(effective_genes)
            f_effective_genes[i] = len(effective_genes) / regulated_genes_num
        else:
            n_effective_genes[i] = 0
            f_effective_genes[i] = 0
    #Select BeReTa targets
    BeReTa_metrics = np.column_stack((beneficial_score, beneficial_score_pval, n_effective_genes, n_effective_rxns, f_effective_genes))
    
    beneficial_scoretmp=np.zeros((beneficial_score.shape[0]))
    #print(beneficial_scoretmp.shape)
    j=0
    for value in np.nditer(beneficial_score):
        #print(value)
        beneficial_scoretmp[j] = abs(value)
        j=j+1
    #print(beneficial_scoretmp)
    score_order = np.argsort(-beneficial_scoretmp)
    #print(score_order)
    # 将索引数组转换为整数类型
    score_order = score_order.astype(int)
    # 将转换后的整数类型索引数组转换为普通的Python列表
    score_order = score_order.tolist()
    reg_list2 = [reg_list[i] for i in score_order]
    BeReTa_metrics = BeReTa_metrics[score_order, :]
    
    BeReTaSolution_reg_list = []
    BeReTaSolution_metrics = []
    for i in range(len(reg_list2)):
        if (BeReTa_metrics[i, 1] < pval_cut) and (BeReTa_metrics[i, 2] >= n_target_cut) and (BeReTa_metrics[i, 3] >= n_target_cut) and (BeReTa_metrics[i, 4] >= f_target_cut):
            BeReTaSolution_reg_list.append(reg_list2[i])
            BeReTaSolution_metrics.append(BeReTa_metrics[i, :])
    #print(BeReTaSolution_reg_list)
    #print(BeReTaSolution_metrics)
    #BeReTaSolution={}
    #BeReTaSolution['reglist']=BeReTaSolution_reg_list
    #BeReTaSolution['metrics']=BeReTaSolution_metrics
    #BeReTaSolution = np.column_stack((BeReTaSolution_reg_list, BeReTaSolution_metrics[:, :2])).tolist()
    # 定义列名
    columns = ['beneficial_score', 'beneficial_score_pval', 'n_effective_genes', 'n_effective_rxns', 'f_effective_genes']
    # 定义索引
    index_values = BeReTaSolution_reg_list
    # 使用pandas的DataFrame构造函数
    BeReTaSolution = pd.DataFrame(BeReTaSolution_metrics, columns=columns)
    # 设置DataFrame的索引
    BeReTaSolution.index = index_values
    return BeReTaSolution
