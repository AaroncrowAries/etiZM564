from etgems_web.script.io.json import *
from etgems_web.script.analysis.flux_balance_analysis import *
import pandas as pd
import numpy as np
import datetime
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pyomo.environ as pyo
from concurrent.futures import ProcessPoolExecutor, as_completed
from .max_driving_force import *
from .metabolite_concentration_variability_analysis import *
from .enzyme_cost_variability_analysis import *
from ..plotting.plotting import *
from ..cave.spav_function import *
from ..core.reaction import *
from .coupling_strategy import *
from .continuous_reaction import *
from .BeReTa import *
from ..core.model import *
import re
import warnings
import copy

# 字典append进dataframe老是报futurewarning
warnings.simplefilter(action='ignore', category=FutureWarning)

O2_list=['EX_o2_e_reverse','EX_o2_e','EX_cpd00007_b','EX_o2(e)','EX_MNXM735438@BOUNDARY','mnxr102090c2b','Ex_o2','Ex_o2_reverse','R1287']
ATPM_list=['ATPM','rxn00062','rxn00062_c0','DM_atp_c_','mnxr153054c3','R00086','META:ATPASE-RXN','ATPASE-RXN','MNXR96131','rxn09694', 'rxn11300','rxn09694_c0', 'rxn11300_c0','ATPM_c', 'ATPM_LPAREN_NGAM_RPAREN_', 'DM_atp_c', 'DM_atp[c]']
CO2_list=['EX_co2', 'EX_co2_e', 'EX_co2_e_', 'EX_co2(e)', 'EX_co2[e]', 'EX_co2_LPAREN_e_RPAREN_','EX_co2_reverse', 'EX_co2_e_reverse','META:TRANS-RXN0-545','MNXR96810','rxn05467', 'rxn08237', 'rxn08238', 'rxn09706', 'rxn09775', 'rxn09821', 'rxn09860', 'rxn09876','R-ATH-1237042', 'R-ATH-1247649', 'R-BTA-1237042', 'R-BTA-1237069', 'R-BTA-1247645', 'R-BTA-1247649', 'R-CEL-1237069', 'R-CEL-1247645', 'R-CFA-1237042', 'R-CFA-1237069', 'R-CFA-1247645', 'R-CFA-1247649', 'R-DDI-1237042', 'R-DDI-1237069', 'R-DDI-1247645', 'R-DDI-1247649', 'R-DME-1237042', 'R-DME-1237069', 'R-DME-1247645', 'R-DME-1247649', 'R-DRE-1237042', 'R-DRE-1237069', 'R-DRE-1247645', 'R-DRE-1247649', 'R-GGA-1237042', 'R-GGA-1237069', 'R-GGA-1247645', 'R-GGA-1247649', 'R-HSA-1237042', 'R-HSA-1237069', 'R-HSA-1247645', 'R-HSA-1247649', 'R-MMU-1237042', 'R-MMU-1237069', 'R-MMU-1247645', 'R-MMU-1247649', 'R-OSA-1237042', 'R-OSA-1247649', 'R-RNO-1237042', 'R-RNO-1237069', 'R-RNO-1247645', 'R-RNO-1247649', 'R-SCE-1237042', 'R-SCE-1247649', 'R-SSC-1237042', 'R-SSC-1237069', 'R-SSC-1247645', 'R-SSC-1247649', 'R-TGU-1237042', 'R-TGU-1237069', 'R-TGU-1247645', 'R-TGU-1247649', 'R-XTR-1237069', 'R-XTR-1247645']
ATP_list=['atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J']
nadph_list=['nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0']
nadh_list=['nadh_c','NADH','C00004','C00004[c]','MNXM10','HMDB01487','cpd00004','BOPGDPNILDQYTO-NNYOXOHSSA-L','nadh','META:NADH','cpd00004_c0']
ATPM_id = 'ATPM'
o2_model_id = 'EX_o2_e'

def set_model_param(norm_model,model_substrate_file,model_product_file,substrate_id,product_id,substrate_uptake_rate,o2_value,ATPM,allow_co2):
    # 清除原先模型的medium
    medium_reaction=[]
    for key,value in norm_model.get_medium().items():
        medium_reaction.append(key)

    for each_r in medium_reaction:
        norm_model.remove_medium(each_r)  
      
    model_substrate = pd.read_csv(model_substrate_file,sep='\t',index_col='substrate_model_correspondind')
    model_product = pd.read_csv(model_product_file, sep='\t', index_col="obj_model_correspondind")  
    
    #substrate_name = model_substrate.loc[substrate_id,'substrate_name']
    #product_name = model_product.loc[product_id,'obj_name']
    product_m_id=model_product.loc[product_id,'obj_id']
    substrate_m_id=model_substrate.loc[substrate_id,'substrate_id']
    #print(product_id)
    if re.search('user_add_DM_',product_id):
        try:
            #demand = norm_model.add_demand_reaction(product_m_id)#这样操作后demand.id没有继承reaction的内容
            norm_model.add_demand_reaction(product_m_id)
            demand = norm_model.get_reaction_by_id("DM_"+product_m_id)
            norm_model.obj_name = demand.id #add demand reaction as the objective
            #print(demand.id)
        except:
            print('Unable to add demand reaction for the objective metabolite, please check the objective product or the model!')
            data =  {
                "statusCode":500,
                "Error":"Unable to add demand reaction for the objective metabolite, please check the objective product or the model!"
            }
        else:
            #norm_model.obj_name = demand.id
            product_id=demand.id
            #print(product_id)
    else:
        norm_model.obj_name = product_id
    #print(substrate_id,product_id)
    #代谢物浓度
    if substrate_id !='substrate_use_in_model':
        if substrate_id.startswith('SK_'):
            norm_model.add_sink_reaction(substrate_m_id)
            sink_reaction = norm_model.get_reaction_by_id("SK_"+substrate_m_id)
            substrate_id=sink_reaction.id
    norm_model.add_medium(substrate_id,substrate_uptake_rate)

    #calculate_GEM(norm_model,"maximize")
    print(substrate_id,product_id,norm_model.medium)
    #有氧厌氧
    if o2_model_id in O2_list:
        try:
            norm_model.get_reaction_by_id(o2_model_id)
        except:
            print('No this O2 reaction!')                   
        else:
            if re.search('reverse',o2_model_id):
                norm_model.get_reaction_by_id(o2_model_id).lower_bound = 0
                norm_model.get_reaction_by_id(o2_model_id).upper_bound = o2_value
            else:
                norm_model.get_reaction_by_id(o2_model_id).lower_bound = -o2_value
                norm_model.get_reaction_by_id(o2_model_id).upper_bound = o2_value
    else:
        Warning_string="The program failed to find the ID for oxygen uptake reaction, and the parameter you set did not take effect and the default value in the model was used. It is recommended to modify the corresponding reaction ID to 'EX_o2_e' in the 'Model Modification' page."

    #allow CO2 assimilation from environment
    if allow_co2=='YES':
        kk=0
        for each_co2 in CO2_list:
            try:
                norm_model.get_reaction_by_id(each_co2)
            except:
                pass                   
            else:
                kk=1
                norm_model.get_reaction_by_id(each_co2).lower_bound = -1000
                norm_model.get_reaction_by_id(each_co2).upper_bound = 1000
        if kk==0:
            Warning_string="The program failed to find the ID for CO2 assimilation reaction, and the parameter you set did not take effect and the default value in the model was used. It is recommended to modify the corresponding reaction ID to 'EX_co2' in the 'Model Modification' page."
      
    # ATPMs
    if ATPM_id in ATPM_list:
        try:
            norm_model.get_reaction_by_id(ATPM_id)
        except:
            print("This is model does not contain ATP Maintenance Requirement reaction!")
        else:
            norm_model.get_reaction_by_id(ATPM_id).lower_bound = ATPM
            norm_model.get_reaction_by_id(ATPM_id).upper_bound = ATPM
    else:
        Warning_string="The program failed to find the ID for non-growth ATP maintenance (NGAM) reaction, and the parameter you set did not take effect and the default value in the model was used. It is recommended to modify the corresponding reaction ID to 'ATPM' in the 'Model Modification' page."

    return [norm_model,product_id,substrate_id]

def pathway_compute(model_name,local_model_file,model_substrate_file,model_product_file,substrate_id,substrate_uptake_rate,o2_value,ATPM,allow_co2,product_id,model_type,workdir):
    """_summary_

    Args:
        event (dict): _description_
        workdir (str): _description_
    """
    model = read_json_model(local_model_file)
    norm_model = copy.deepcopy(model)
    #模型参数设置
    [norm_model,product_id_new] = set_model_param(norm_model,model_substrate_file,model_product_file,
                substrate_id,product_id,substrate_uptake_rate,o2_value,ATPM,allow_co2)

    # save
    model_change_file = os.path.join(workdir,'%s_modified_model_change.json' % model_name)
    write_json_model(norm_model,model_change_file)

    # 后续计算都用这个product id 
    #event['transform_product_id'] = product_id_new
    if model_type in ["GEM","ecGEM"]:
        # 直接返回结果
        flux_outfile,result_static = pathway_compute_result(norm_model,model_type,model_name,substrate_id,product_id_new,workdir)
        return flux_outfile,result_static
    else:
        # tcGEM etGEM  meETGEM 返回计算拐点
        mdf_file, mdf_plot,inflection_points = pathway_compute_mdf(norm_model,workdir)
        #对每个拐点进行计算
        pathway_compute_set_mdf(norm_model,model_name,inflection_points,product_id_new,workdir)
        return flux_outfile,result_static,mdf_file,mdf_plot
    
def pathway_compute_set_mdf(norm_model: object,model_name:str,inflection_points:pd.DataFrame,product_id:str,substrate_id:str,model_type:str,met_C_file:str,workdir:str):
    """_summary_

    Args:
        event (dict): _description_
        workdir (str): _description_
        
        event = {
            "jobid":"",
            "set_MDF_substrate":""
        }
    """
    #inflection_points_list =list(inflection_points['ETM']+0.05)
    inflection_points_list = [eachp + 0.05 if eachp > 0 else eachp - 0.05 for eachp in list(inflection_points.index)]

    pathway_compute_dict ={}
    #norm_model = read_json_model(local_compute_model_file)
    if len(inflection_points_list)>0:
        with ProcessPoolExecutor(len(inflection_points_list)) as executor:
            futures = {executor.submit(MDF_Calculation,norm_model, eachpoint, product_id): eachpoint for eachpoint in inflection_points_list}
            for future in as_completed(futures):
                B_value = future.result()
                #B_value = MDF_Calculation(norm_model,set_MDF_substrate,product_id)
                print("B_value:",B_value)
                # compute
                flux_outfile,result_static = pathway_compute_result(norm_model,model_type,model_name,substrate_id,product_id,met_C_file,workdir,B_value)
                #print("flux_outfile:",flux_outfile)
                #print("result_static:",result_static)
                pathway_compute_dict[B_value] = {}
                pathway_compute_dict[B_value]["flux_outfile"] = flux_outfile
                pathway_compute_dict[B_value]["result_static"] = result_static

    return pathway_compute_dict
                

def get_me_strategies(norm_model,B_value,biomass_id:str,local_exp_comp_genes_file:str,local_reg_net_file:str,workdir:str):
    exp_comp_genes = pd.read_csv(local_exp_comp_genes_file,index_col=0)
    exp_comp={}
    exp_comp['genes']=exp_comp_genes.index
    exp_comp['expression']=exp_comp_genes
    
    reg_net = pd.read_csv(local_reg_net_file)
    
    [Bottleneck_reaction_result_file,Limiting_metabolite_result_file,Bottleneck_enzyme_result_file,ETM_solution_file,Bottleneck_reaction_list] = me_strategies(norm_model,B_value,reg_net, exp_comp, biomass_id,workdir)
    [Fusion_Protein_file,Energy_Coupling_file,Continuous_Reaction_file] = meetgem_fusion_protein(norm_model,Bottleneck_reaction_list,workdir)
   
    return Bottleneck_reaction_result_file,Limiting_metabolite_result_file,Bottleneck_enzyme_result_file,ETM_solution_file,Fusion_Protein_file,Energy_Coupling_file,Continuous_Reaction_file


def meetgem_fusion_protein(norm_model:object,Bottleneck_reaction_list:list,workdir:str):
    tasks = ["Fusion_Protein","Energy_Coupling","Continuous_Reaction",]
    
    output_files =  []
    for job_type in tasks:
        filename = f'{job_type}.tsv'
        if job_type == "Fusion_Protein":
            result_df = Get_Fusion_Protein(Bottleneck_reaction_list,norm_model)
        elif job_type == "Energy_Coupling":
            result_df = Find_Energy_Coupling_Strategy(Bottleneck_reaction_list,norm_model)
        elif job_type == "Continuous_Reaction":
            result_df = Get_Continuous_Reaction(Bottleneck_reaction_list,norm_model)
        result_file = os.path.join(workdir,filename)
        result_df.to_csv(result_file,sep="\t") 
        output_files.append(result_file)
        
    return output_files

def pathway_compute_mdf(norm_model: object,workdir:str):
    MDF_result_file = os.path.join(workdir,"MDF_clc_result.csv")
    [MDF_list_ETM,MDF_df_ETM] = calculate_MDF_list(norm_model, "ETM", MDF_result_file)
    df = pd.read_csv(MDF_result_file, index_col=0)
    inflection_points = df[df['ETM'].diff() != 0]
    inflection_points = inflection_points.dropna()
    inflection_points.to_csv(MDF_result_file)
    
    MDF_dict = {"ETM" : MDF_list_ETM}  
    if MDF_df_ETM['ETM'].min()<0:
        ytick = list(np.arange(MDF_df_ETM['ETM'].min(), MDF_df_ETM['ETM'].max() + 2, 2))
    else:
        ytick = list(np.arange(0, MDF_df_ETM['ETM'].max() + 2, 2))
        
    png_file = os.path.join(workdir,"max_MDF_by_four_model.png")
    #xlab = "Growth rate ($\mathregular{h^-1}$)"
    xlab = "Product rate (mmol/gDW/h)"
    Draw_MDF_By_Product_rate(MDF_dict,0,xlab, MDF_df_ETM.index.max()+1,1,ytick,png_file)  
    
    return MDF_result_file, png_file,inflection_points
        
def pathway_compute_result(norm_model: object,model_type:str,model_name:str,substrate_id:str,product_id:str,met_C_file:str,workdir:str,B_value=""):
    if model_type == 'GEM':
        opt_etm = calculate_GEM(norm_model,"maximize")
    elif model_type == 'ecGEM':    
        opt_etm = calculate_ECM(norm_model,"maximize")
    elif model_type == 'tcGEM':
        opt_etm = calculate_TCM(norm_model,"maximize",B_value)
    elif model_type == 'etGEM':    
        opt_etm = calculate_ETM(norm_model,"maximize",B_value)
    elif model_type == "meETGEM":
        opt_etm = calculate_ETM(norm_model,"maximize",B_value)

    met_C = pd.read_csv(met_C_file,sep='\t', index_col=0)

    reactions_list = []
    fluxes = []
    for reaction in norm_model.reactions.keys():
        reactions_list.append(reaction)
        fluxes.append(opt_etm.reaction[reaction].value)
    final_data = {'reaction' : reactions_list,
                'fluxes' : fluxes}
    flux_solution_df = pd.DataFrame(final_data)
    flux_solution_df = flux_solution_df.set_index('reaction')
    #flux_solution_df.head()
    
    model_pfba_solution=pd.DataFrame() 
    result_static={}
    # flux_solution_df=flux_solution.to_frame()
    #print(flux_solution_df)
    Optimal_rate=flux_solution_df.loc[product_id,'fluxes']
    flux_solution_df_select=flux_solution_df[abs(flux_solution_df['fluxes'])>0.0000001]#只储存大于0的值e-6
    result_static['Reaction number']=flux_solution_df_select.shape[0]#反应个数
    result_static['Objective flux']=round(Optimal_rate,7)

    if flux_solution_df_select.shape[0] > 0 and Optimal_rate>0:
        result_static['Exchange reaction']={}
        result_static['ATP consumption reaction']={}
        result_static['ATP generation reaction']={}
        result_static['NADH consumption reaction']={}
        result_static['NADH generation reaction']={}
        result_static['NADPH consumption reaction']={}
        result_static['NADPH generation reaction']={}

        for index, row in flux_solution_df_select.iterrows():
            model_pfba_solution.loc[index,'reaction_id']=index
            model_pfba_solution.loc[index,'reaction_name']=norm_model.get_reaction_by_id(index).name
            # model_pfba_solution.loc[index,'fluxes']=abs(round(row['fluxes'],3))
            model_pfba_solution.loc[index,'fluxes']=round(row['fluxes'],7)
            model_pfba_solution.loc[index,'gene id']=norm_model.get_reaction_by_id(index).gene_reaction_rule
            #model_pfba_solution.loc[index,'ec-code']=norm_model.reactions.get_by_id(index).annotation['ec-code'] 
            model_pfba_solution.loc[index,'equ']=norm_model.get_reaction_by_id(index).build_reaction_string()
            model_pfba_solution.loc[index,'equ']=norm_model.get_reaction_by_id(index).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
            model_pfba_solution.loc[index,'abs_fluxes']=abs(round(row['fluxes'],7))
            
            met_generation_reaction(norm_model,index,ATP_list,'ATP',flux_solution_df_select,result_static)    
            met_generation_reaction(norm_model,index,nadh_list,'NADH',flux_solution_df_select,result_static)  
            met_generation_reaction(norm_model,index,nadph_list,'NADPH',flux_solution_df_select,result_static)  
            met_consumption_reaction(norm_model,index,ATP_list,'ATP',flux_solution_df_select,result_static)    
            met_consumption_reaction(norm_model,index,nadh_list,'NADH',flux_solution_df_select,result_static)  
            met_consumption_reaction(norm_model,index,nadph_list,'NADPH',flux_solution_df_select,result_static)

        for eachr in norm_model.reactions.keys():
            if "EX_" in eachr:
                if eachr in list(flux_solution_df_select.index):
                    result_static['Exchange reaction'][eachr]={}
                    result_static['Exchange reaction'][eachr]['reaction_id']=eachr
                    result_static['Exchange reaction'][eachr]['reaction_name']=norm_model.get_reaction_by_id(eachr).name
                    # result_static['Exchange reaction'][eachr.id]['fluxes']=abs(round(flux_solution_df_select.loc[eachr.id,'fluxes'],3))#交换反应速率
                    result_static['Exchange reaction'][eachr]['fluxes']=round(flux_solution_df_select.loc[eachr,'fluxes'],7)#交换反应速率
                    result_static['Exchange reaction'][eachr]['equ']=norm_model.get_reaction_by_id(eachr).build_reaction_string()
                    result_static['Exchange reaction'][eachr]['equ_name']=norm_model.get_reaction_by_id(eachr).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
                if substrate_id in list(flux_solution_df_select.index):
                    result_static['Exchange reaction'][substrate_id]={}
                    result_static['Exchange reaction'][substrate_id]['reaction_id']=substrate_id
                    result_static['Exchange reaction'][substrate_id]['reaction_name']=norm_model.get_reaction_by_id(substrate_id).name
                    # result_static['Exchange reaction'][eachr.id]['fluxes']=abs(round(flux_solution_df_select.loc[eachr.id,'fluxes'],3))#交换反应速率
                    result_static['Exchange reaction'][substrate_id]['fluxes']=round(flux_solution_df_select.loc[substrate_id,'fluxes'],3)#交换反应速率
                    result_static['Exchange reaction'][substrate_id]['equ']=norm_model.get_reaction_by_id(substrate_id).build_reaction_string()
                    result_static['Exchange reaction'][substrate_id]['equ_name']=norm_model.get_reaction_by_id(substrate_id).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
        #产物碳原子数---其实底物和产物这个可以写一个函数放spav_function.py的                          
        product_c_num_r=0
        product_c_num_l=0
        product_equ=norm_model.get_reaction_by_id(product_id)
        # product_equ = norm_model.get_reaction_by_id("ATPM")
    for eachm in product_equ.metabolites.keys():
        if product_equ.metabolites[eachm]<0:
            met_coef=abs(product_equ.metabolites[eachm])
            try:
                met_formula=get_composition_from_formula(norm_model.get_metabolite_by_id(eachm).formula)
            except:
                try:
                    met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                except:
                    print('No formula!')
                else:
                    if 'C' in met_formula.keys():
                        #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                        product_c_num_r=product_c_num_r+met_coef*met_formula['C']
            else:
                if 'C' in met_formula.keys():
                    #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                    product_c_num_r=product_c_num_r+met_coef*met_formula['C']

    for eachm in product_equ.metabolites.keys():
        if product_equ.metabolites[eachm]>0:
            met_coef=abs(product_equ.metabolites[eachm])
            try:
                met_formula=get_composition_from_formula(norm_model.get_metabolite_by_id(eachm).formula)
            except:
                try:
                    met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                except:
                    print('No formula!')
                else:
                    if 'C' in met_formula.keys():
                        #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                        product_c_num_r=product_c_num_r+met_coef*met_formula['C']
            else:
                if 'C' in met_formula.keys():
                    #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                    product_c_num_r=product_c_num_r+met_coef*met_formula['C']

        product_c_num=product_c_num_r-product_c_num_l
        #只有产品碳原子有净生成时计算 
        if product_c_num>0:
            #底物碳原子数
            substrate_c_num_r=0
            substrate_c_num_l=0
            substrate_equ=norm_model.get_reaction_by_id(substrate_id)
            for eachm in substrate_equ.metabolites.keys():
                if substrate_equ.metabolites[eachm]<0:
                    met_coef=abs(substrate_equ.metabolites[eachm])
                    try:
                        met_formula=get_composition_from_formula(norm_model.get_metabolite_by_id(eachm).formula)
                    except:
                        try:
                            met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                        except:
                            print('Substrate right no formula!')
                        else:
                            if 'C' in met_formula.keys():
                                #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                                substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                    else:
                        if 'C' in met_formula.keys():
                            #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                            substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                        
            for eachm in substrate_equ.metabolites.keys():
                if substrate_equ.metabolites[eachm]>0:
                    met_coef=abs(substrate_equ.metabolites[eachm])
                    try:
                        met_formula=get_composition_from_formula(norm_model.get_metabolite_by_id(eachm).formula)
                    except:
                        try:
                            met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                        except:
                            print('Substrate right no formula!')
                        else:
                            if 'C' in met_formula.keys():
                                #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                                substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                    else:
                        if 'C' in met_formula.keys():
                            #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                            substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']                                     
            substrate_c_num=substrate_c_num_r-substrate_c_num_l

            if substrate_id in list(flux_solution_df_select.index):
                result_static['C Yield']=abs(round((product_c_num*flux_solution_df_select.loc[product_id,'fluxes'])/(substrate_c_num*flux_solution_df_select.loc[substrate_id,'fluxes']),7))
            elif substrate_id+"_reverse" in list(flux_solution_df_select.index):
                result_static['C Yield']=abs(round((product_c_num*flux_solution_df_select.loc[product_id,'fluxes'])/(substrate_c_num*flux_solution_df_select.loc[substrate_id+"_reverse",'fluxes']),7))
                
        #计算得率
        try:
            flux_solution_df_select.loc[product_id,'fluxes']
        except:
            result_static['Yield']=flux_solution_df.loc[product_id,'fluxes']
        else:
            if substrate_id in list(flux_solution_df_select.index):
                result_static['Yield']=abs(round(flux_solution_df_select.loc[product_id,'fluxes']/flux_solution_df_select.loc[substrate_id,'fluxes'],7))
            elif substrate_id+"_reverse" in list(flux_solution_df_select.index):
                result_static['Yield']=abs(round(flux_solution_df_select.loc[product_id,'fluxes']/flux_solution_df_select.loc[substrate_id+"_reverse",'fluxes'],7))
    if model_type == "ecGEM":
        for index,row in model_pfba_solution.iterrows():
            reaction = norm_model.get_reaction_by_id(row["reaction_id"])
            model_pfba_solution.loc[index,"E"] = opt_etm.reaction[reaction.id].value*reaction.mw/reaction.kcat/3600000
    if model_type == "tcGEM":
        for index,row in model_pfba_solution.iterrows():
            reaction = norm_model.get_reaction_by_id(row["reaction_id"])
            metabolite_concentration = {}
            for metabolite in reaction.metabolites:
                metabolite_concentration[metabolite] = opt_etm.metabolite[metabolite].value
            model_pfba_solution.loc[index,"met_concentration"] = str(metabolite_concentration).replace("'","")
            if reaction.g0 != -99999:
                model_pfba_solution.loc[index,"MDF"] = -reaction.g0-2.579*sum(metabolite_coe*opt_etm.metabolite[metabolite].value  for metabolite,metabolite_coe in reaction.metabolites.items())
            else:
                model_pfba_solution.loc[index,"MDF"] = -9999
    if model_type == "etGEM" or model_type == "meETGEM":
        for index,row in model_pfba_solution.iterrows():
            reaction = norm_model.get_reaction_by_id(row["reaction_id"])
            metabolite_concentration = {}
            for metabolite in reaction.metabolites:
                metabolite_concentration[metabolite] = opt_etm.metabolite[metabolite].value
            model_pfba_solution.loc[index,"met_concentration"] = str(metabolite_concentration).replace("'","")
            if reaction.g0 != -99999:
                model_pfba_solution.loc[index,"MDF"] = -reaction.g0-2.579*sum(metabolite_coe*opt_etm.metabolite[metabolite].value  for metabolite,metabolite_coe in reaction.metabolites.items())
            else:
                model_pfba_solution.loc[index,"MDF"] = -9999
        model_pfba_solution.loc[index,"E"] = opt_etm.reaction[reaction.id].value*reaction.mw/reaction.kcat/3600000
    # output
    flux_outfile=os.path.join(
        workdir,
        f'{model_name}-d3flux.tsv'
    )
    model_pfba_solution=model_pfba_solution.sort_values(by=['abs_fluxes'],ascending=False)
    model_pfba_solution.to_csv(flux_outfile,index=False,sep='\t')

    return flux_outfile,result_static

def me_strategies(norm_model: object,B_value,reg_net:pd.DataFrame, exp_comp:pd.DataFrame, biomass_id:str,workdir:str):
    pval_cut = 0.05 #(2) The p-value of the beneficial score should be less than 0.05.
    n_target_cut = 1 #(3) TR should have two or more effective gene/reaction targets.
    f_target_cut = 0.1 #(4) At least 10% of target metabolic genes of TR should be beneficial, i.e. have positive flux slopes.

    norm_model_me = copy.deepcopy(norm_model)
    opt_etm = calculate_ETM(norm_model,"maximize",B_value)

    #pFBA
    [min_V,Concretemodel]=Min_Flux_Sum_Calculation(norm_model_me,"ETM",opt_etm.obj(),norm_model.obj_name,B_value)
    print("Min flux amount : " +str(min_V))
    use_result = Get_Results_Thermodynamics(norm_model_me,Concretemodel)
    use_result = use_result[use_result['flux'] > 1e-10] 
    use_result = use_result.sort_values(by = 'flux',axis = 0,ascending = False)
    use_result["reaction"] = use_result.apply(lambda row: norm_model_me.get_reaction_by_id(row.name).build_reaction_string(), axis = 1)
    use_result["gpr"] = use_result.apply(lambda row: norm_model_me.get_reaction_by_id(row.name).gene_reaction_rule, axis = 1)
    use_result_tmp=use_result[use_result['f']>-norm_model_me.K_value]
    use_result_select=use_result_tmp[abs(use_result_tmp['f']-B_value)<=1e-05]
    path_reac_list=list(use_result_select.index)

    #基于计算结果，挖关键酶
    use_result_sort = use_result.sort_values(by='enz',ascending = False)
    use_result_sort.head()
    e_threshold = norm_model.E_total*0.01
    enz_use_reaction_list = list(use_result_sort[use_result_sort['enz'] > e_threshold].index)
    enz_use_reaction_list
    # Calculate the minimum enzyme cost for reactions in above list, and the higher the value, the more critical (key) the enzyme.
    max_min_E_list_fixed = max_Max_Min_E_for_list(norm_model_me,B_value,norm_model.obj_name,obj_enz_constraint,enz_use_reaction_list,4)
    max_min_E_list_fixed = max_min_E_list_fixed.sort_values(by=['min_E'],ascending=False)
    Bottleneck_enzyme = max_min_E_list_fixed[(max_min_E_list_fixed['max_E'] - max_min_E_list_fixed['min_E']) <= 0.001]
    column_names = ["name", "GPR", "min_E", "max_E", "kcat", "mw", "g0"]
    Bottleneck_enzyme_result = pd.DataFrame(columns=column_names)
    for index in Bottleneck_enzyme.index:
        limit_enzyme = norm_model_me.get_reaction_by_id(index)
        Bottleneck_enzyme_result.loc[index,"name"] = limit_enzyme.name
        Bottleneck_enzyme_result.loc[index,"GPR"] = limit_enzyme.gene_reaction_rule
        Bottleneck_enzyme_result.loc[index,"min_E"] = Bottleneck_enzyme.loc[index,"min_E"]
        Bottleneck_enzyme_result.loc[index,"max_E"] = Bottleneck_enzyme.loc[index,"max_E"] 
        Bottleneck_enzyme_result.loc[index,"kcat"] = limit_enzyme.kcat
        Bottleneck_enzyme_result.loc[index,"mw"] = limit_enzyme.mw
        Bottleneck_enzyme_result.loc[index,"g0"] = limit_enzyme.g0
    Bottleneck_enzyme_result_file = os.path.join(workdir,"Bottleneck_enzyme_result.csv")

    # 挖调控因子靶点
    norm_model_me = copy.deepcopy(norm_model)
    ETM_solution = BeReTa_use(norm_model_me,'etGEM', reg_net, exp_comp, biomass_id,norm_model.obj_name, pval_cut, n_target_cut, f_target_cut)
    for index, row in ETM_solution.iterrows():
        if row["beneficial_score"] < 0:
            ETM_solution.loc[index,"Strategy"] = "Weaken"
        elif row["beneficial_score"] > 0:
            ETM_solution.loc[index,"Strategy"] = "Enhance"
    ETM_solution_file = os.path.join(workdir,"ETM_solution.csv")  
    
    # Calculate the maximum thermodynamic driving force for reactions in above list, if the value is still equal to B, it is the bottleneck reaction.
    reaction_Df_list_fixed=Get_Max_Min_Df_Complete(norm_model_me,"maximize",path_reac_list,B_value,opt_etm.obj(),norm_model.obj_name)
    print(reaction_Df_list_fixed)
    Bottleneck_reaction=reaction_Df_list_fixed[(reaction_Df_list_fixed['maximize']-B_value)<=0.01]
    print(Bottleneck_reaction)
    column_names = ['id', 'name', 'GPR', 'lower_bound', 'upper_bound', 'kcat', 'mw', 'g0']
    Bottleneck_reaction_result = pd.DataFrame(columns=column_names)
    column_names = ["name", "charge", "formula", "concentration_lb", "concentration_ub","concentration", "Bottleneck_reaction"]
    Limiting_metabolite_result = pd.DataFrame(columns=column_names)

    if len(list(Bottleneck_reaction.index))>0:#如果有瓶颈反应
        Bottleneck_reaction_list = list(Bottleneck_reaction.index)
        for bottleneck in Bottleneck_reaction_list:
            bottleneck_reaction = norm_model_me.get_reaction_by_id(bottleneck)
            Bottleneck_reaction_result.loc[bottleneck,"id"] = bottleneck_reaction.id
            Bottleneck_reaction_result.loc[bottleneck,"name"] = bottleneck_reaction.name
            Bottleneck_reaction_result.loc[bottleneck,"GPR"] = bottleneck_reaction.gene_reaction_rule
            Bottleneck_reaction_result.loc[bottleneck,"lower_bound"] = bottleneck_reaction.lower_bound
            Bottleneck_reaction_result.loc[bottleneck,"upper_bound"] = bottleneck_reaction.upper_bound
            Bottleneck_reaction_result.loc[bottleneck,"kcat"] = bottleneck_reaction.kcat
            Bottleneck_reaction_result.loc[bottleneck,"mw"] = bottleneck_reaction.mw
            Bottleneck_reaction_result.loc[bottleneck,"g0"] = bottleneck_reaction.g0

        #有了瓶颈反应，计算限制性代谢物
        norm_model_me = copy.deepcopy(norm_model)
        if len(list(Bottleneck_reaction.index))>0:
            use_result_select.loc[Bottleneck_reaction.index[0],'met_concentration']
        Bottleneck_reaction_list=list(Bottleneck_reaction.index)
        Bottleneck_reaction_met=[]
        for rea in norm_model_me.reactions.values():
            if rea.id in Bottleneck_reaction_list:
                #print(rea)
                for met in rea.metabolites.keys():
                    if met !='h_c' and met !='h2o_c':
                        Bottleneck_reaction_met.append(met)         
        Bottleneck_reaction_met=list(set(Bottleneck_reaction_met))
        print(Bottleneck_reaction_met)
        max_min_concentration_list_fixed = Get_Max_Min_Met_Concentration_for_list(norm_model_me,B_value,obj_enz_constraint,norm_model.obj_name,Bottleneck_reaction_list,Bottleneck_reaction_met,4)
        print(max_min_concentration_list_fixed)
        Limiting_metabolite = max_min_concentration_list_fixed[(max_min_concentration_list_fixed['max_concentration'] - max_min_concentration_list_fixed['min_concentration']) <= 0.001]
        if len(list(Limiting_metabolite.index))>0:
            for index in Limiting_metabolite.index:
                limit_metabolite = norm_model_me.get_metabolite_by_id(index)
                Limiting_metabolite_result.loc[index,"name"] = limit_metabolite.name
                Limiting_metabolite_result.loc[index,"charge"] = limit_metabolite.charge
                Limiting_metabolite_result.loc[index,"formula"] = limit_metabolite.formula
                Limiting_metabolite_result.loc[index,"concentration_lb"] = limit_metabolite.concentration_lb 
                Limiting_metabolite_result.loc[index,"concentration_ub"] = limit_metabolite.concentration_ub
                Limiting_metabolite_result.loc[index,"concentration"] = Limiting_metabolite.loc[index, 'max_concentration']
                corr_reaction = []
                for reaction in Bottleneck_reaction_list:
                    reaction = norm_model_me.get_reaction_by_id(reaction)
                    if index in reaction.metabolites.keys():
                        corr_reaction.append(reaction.id)
                Limiting_metabolite_result.loc[index,"Bottleneck_reaction"] = corr_reaction

    Bottleneck_reaction_result_file = os.path.join(workdir,"Bottleneck_reaction_result.csv")
    Limiting_metabolite_result_file = os.path.join(workdir,"Limiting_metabolite_result.csv")
    return Bottleneck_reaction_result_file,Limiting_metabolite_result_file,Bottleneck_enzyme_result_file,ETM_solution_file,Bottleneck_reaction_list
