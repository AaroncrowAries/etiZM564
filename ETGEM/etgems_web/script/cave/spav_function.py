import json
import pandas as pd
import cobra
import itertools
import re
import os
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import set_objective
from cobra.manipulation.delete import knock_out_model_genes
from cobra.medium import find_boundary_types, find_external_compartment, sbo_terms
import matplotlib.pyplot as plt
import numpy as np
from cobra.core.model import Model as Model2
from cobra.core.reaction import Reaction as Reaction2
import d3flux

def download_reference_from_s3(s3,bucket,obj):
    """
        download reference data /tmp/
    """
    object_name = obj.split('/')[-1]
    local_reference = os.path.join('/tmp',object_name)
    s3.Object(bucket, obj).download_file(local_reference)
    return local_reference

def download_file_from_s3(s3,bucket,obj,workdir):
    """
        download file data workdir
    """
    object_name = obj.split('/')[-1]
    local_reference = os.path.join(workdir,object_name)
    s3.Object(bucket, obj).download_file(local_reference)
    return local_reference

def change_file_acl(s3,bucket,object_key):
    """
    修改文件的acl，便于客户下载
    """
    try:
        s3.ObjectAcl(bucket,object_key).put(ACL="public-read")
        url = "https://"+bucket+".s3.amazonaws.com/%s" % object_key
        return url
    except Exception as e:
        print(str(e))
        print('%s is not existed in s3' %object_key )
        return False

def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def get_fluxes_detail_in_model(model, call_method, fluxes_outfile):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """
    try:
        if call_method=='pFBA':
            model_pfba_solution = cobra.flux_analysis.pfba(model)
        elif call_method=='FBA':   
            model_pfba_solution = model.optimize()

    except:
        print('Can not solve the model!')
    else:
        model_pfba_solution = model_pfba_solution.to_frame()
        model_pfba_solution_detail = {}
        for index, row in model_pfba_solution.iterrows():
            reaction_detail = model.reactions.get_by_id(index)
            model_pfba_solution_detail[index]={}
            if abs(row['fluxes'])>1e-9:
                useflux=round(row['fluxes'],6)
            else:
                useflux=0
                
            model_pfba_solution_detail[index]['fluxes'] = useflux    
            #for eachkey in reaction_detail.annotation.keys():
            #    model_pfba_solution_detail[index][eachkey] = reaction_detail.annotation[eachkey]   
            if 'ec-code' in reaction_detail.annotation.keys():
                model_pfba_solution_detail[index]['ec-code'] = reaction_detail.annotation['ec-code']  
            else:
                model_pfba_solution_detail[index]['ec-code'] = 'none'
            if reaction_detail.gene_reaction_rule:
                model_pfba_solution_detail[index]['gene id'] =reaction_detail.gene_reaction_rule
            else:
                model_pfba_solution_detail[index]['gene id'] = 'none'
            model_pfba_solution_detail[index]['equ'] = reaction_detail.reaction
            
        model_pfba_solution_detail_df = pd.DataFrame(model_pfba_solution_detail).T
        model_pfba_solution_detail_df=model_pfba_solution_detail_df[abs(model_pfba_solution_detail_df['fluxes'])>0]
        model_pfba_solution_detail_df.to_csv(fluxes_outfile,sep='\t')
        return model_pfba_solution_detail_df
    
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
        
def get_metpair(rea, pi_pairs1, h_pairs1, pi_pairs2, h_pairs2, nh4_pairs, other_pairs, currency_mets):
    #print(rea.id)
    # get the metabolite links for a reaction, excluding links through currency metabolites
    # processing in the order of P & H transfer, N transfer and other transfers
    sub_pro = []
    mark = 0  # mark if there is currency metabolite pairs for P/H transfer in the reaction
    c2mark = 0  # mark if there is currency metabolite pairs for P/H transfer in the reaction, second batch
    nmark = 0  # mark if there is nh4 transfer currency metabolite pairs in the reaction, deal seperately
    omark = 0  # mark if there is other group transfer currency metabolite pairs in the reaction, deal seperately
    cmet = []  # temorary currency metabolite list
    c2met = []  # temorary currency metabolite list, second batch
    ncmet = []  # temorary currency metabolite list for N transfer
    ocmet = []  # temorary currency metabolite list for other pair transfer
    phmet = []  # metabolite for P/H transfer, to be excluded again in the last step if >1 pairs still remaining, especially for UDP ADPsugars
    ex_pairs1 = pi_pairs1+h_pairs1  # excluded currency metabolite pairs, first batch
    ex_pairs2 = pi_pairs2+h_pairs2  # excluded currency metabolite pairs, second batch
    for sp in ex_pairs1:
        phmet.append(sp[0])
        phmet.append(sp[1])
    for sp in ex_pairs2:
        phmet.append(sp[0])
        phmet.append(sp[1])
    phmet = list(set(phmet))  # remove repeats
    # also check if CoA is a currency metabolite in the last step
    phmet += ['coa_c', 'coa_p', 'coa_e']
    subs = [
        m.id for m in rea.reactants if 'C' in m.elements if m.id not in currency_mets]
    pros = [m.id for m in rea.products if 'C' in m.elements if m.id not in currency_mets]
    for s, p in itertools.product(subs, pros):
        if (s, p) in ex_pairs1 or (p, s) in ex_pairs1:  # need to consider direction
            mark = 1
            cmet.append(s)
            cmet.append(p)
        if (s, p) in ex_pairs2 or (p, s) in ex_pairs2:  # need to consider direction
            c2mark = 1
            c2met.append(s)
            c2met.append(p)
        if (s, p) in nh4_pairs or (p, s) in nh4_pairs:  # for nh4 transfer
            nmark = 1
            ncmet.append(s)
            ncmet.append(p)
        if (s, p) in other_pairs or (p, s) in other_pairs:  # for other pair transfer
            omark = 1
            ocmet.append(s)
            ocmet.append(p)
    # if len(sub_pro)>1 and mark==1:
    if mark == 1:  # process in order
        subs = [m for m in subs if m not in cmet]
        pros = [m for m in pros if m not in cmet]
    if c2mark == 1:  # process in order
        subsn = [m for m in subs if m not in c2met]
        prosn = [m for m in pros if m not in c2met]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if nmark == 1:
        subsn = [m for m in subs if m not in ncmet]
        prosn = [m for m in pros if m not in ncmet]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if omark == 1:
        subsn = [m for m in subs if m not in ocmet]
        prosn = [m for m in pros if m not in ocmet]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if len(subs) > 1:  # to remove UTP, UDP etc.
        subsn = subs
        for m in subsn:
            if m in phmet:
                subs.remove(m)
                if len(subs) == 1:
                    break
    if len(pros) > 1:  # to remove UTP, UDP etc.
        prosn = pros
        for m in prosn:
            if m in phmet:
                pros.remove(m)
                if len(pros) == 1:
                    break
    for s, p in itertools.product(subs, pros):
        sub_pro.append((s, p))
    return sub_pro

def remove_currency_metabolite(pairs_of_currency_metabolites_file,special_currency_metabolites_file):
    pairs_of_currency_metabolites = pd.read_csv(pairs_of_currency_metabolites_file)
    special_currency_metabolites = pd.read_csv(special_currency_metabolites_file)

    pi_pairs1 = [tuple(x.split(",")) for x in pairs_of_currency_metabolites['pi_pairs1'] if str(x) != 'nan']
    pi_pairs2 = [tuple(x.split(",")) for x in pairs_of_currency_metabolites['pi_pairs2'] if str(x) != 'nan']
    h_pairs1 = [tuple(x.split(",")) for x in pairs_of_currency_metabolites['h_pairs1'] if str(x) != 'nan']
    h_pairs2 = [tuple(x.split(",")) for x in pairs_of_currency_metabolites['h_pairs2'] if str(x) != 'nan']
    nh4_pairs = [tuple(x.split(",")) for x in pairs_of_currency_metabolites['nh4_pairs'] if str(x) != 'nan']
    other_pairs = [tuple(x.split(",")) for x in pairs_of_currency_metabolites['other_pairs'] if str(x) != 'nan']

    excluded=[x for x in special_currency_metabolites['excluded'] if str(x) != 'nan']
    nocarbon=[x for x in special_currency_metabolites['no carbon metabolites'] if str(x) != 'nan'] 
    ions=[x for x in special_currency_metabolites['ions'] if str(x) != 'nan'] 
    currency_mets = excluded+ nocarbon +ions
    
    return(pi_pairs1,h_pairs1,pi_pairs2,h_pairs2,nh4_pairs,other_pairs,currency_mets)

def get_metlist_and_metlinklist_cb(model,model_type,select_reaction_list_df,pairs_of_currency_metabolites_file,special_currency_metabolites_file):
    [pi_pairs1,h_pairs1,pi_pairs2,h_pairs2,nh4_pairs,other_pairs,currency_mets]=remove_currency_metabolite(pairs_of_currency_metabolites_file,special_currency_metabolites_file)
    select_reaction_list=list(select_reaction_list_df.index)
    metabolite_list=[]
    metlinks={} 
    
    for rea in model.reactions:        
        if rea.id in select_reaction_list:
            #去除流通代谢物
            sub_pro = get_metpair(rea,pi_pairs1,h_pairs1,pi_pairs2,h_pairs2,nh4_pairs,other_pairs,currency_mets)
            for sp in sub_pro:
                #if re.search('_e',sp[0]) or re.search('_e',sp[1]):
                #    pass
                #elif re.search('_p',sp[0]) or re.search('_p',sp[1]):
                #    pass
                #else:
                    metname_s=model.metabolites.get_by_id(sp[0]).name
                    metabolite_list.append(sp[0])
                    metname_p=model.metabolites.get_by_id(sp[1]).name
                    metabolite_list.append(sp[1])
                    #flux: mmol gDWh−1 E: mg gDW−1                 
                    key_id=rea.id+'_'+sp[0]+'_'+sp[1]
                    
                    metlinks[key_id]={}
                    metlinks[key_id]['source_node']=sp[0]
                    metlinks[key_id]['target_node']=sp[1]
                    metlinks[key_id]['source_node_hovering']=metname_s
                    metlinks[key_id]['target_node_hovering']=metname_p
                    metlinks[key_id]['reaction_id']=rea.id
                    metlinks[key_id]['reaction_equation']=select_reaction_list_df.loc[rea.id ,'equ']                        
                    metlinks[key_id]['EC_number'] = select_reaction_list_df.loc[rea.id ,'ec-code']                     
                    if model_type=='single_model': 
                        metlinks[key_id]['gene']=select_reaction_list_df.loc[rea.id ,'gene id']
                        if isinstance(select_reaction_list_df.loc[rea.id ,'fluxes'],str):
                            vis_flux=select_reaction_list_df.loc[rea.id ,'fluxes']
                        else:
                            vis_flux=str(abs(round(select_reaction_list_df.loc[rea.id ,'fluxes'],7)))
                        metlinks[key_id]['flux']=vis_flux
                    elif model_type=='model_comparison': 
                        metlinks[key_id]['gene_in_model1']=select_reaction_list_df.loc[rea.id ,'gene in model1']
                        metlinks[key_id]['gene_in_model2']=select_reaction_list_df.loc[rea.id ,'gene in model2']
                        if isinstance(select_reaction_list_df.loc[rea.id ,'fluxes in model1'],str):
                            vis_flux=select_reaction_list_df.loc[rea.id ,'fluxes in model1']
                        else:
                            vis_flux=str(abs(round(select_reaction_list_df.loc[rea.id ,'fluxes in model1'],7)))
                        metlinks[key_id]['flux_in_model1']=vis_flux
                        if isinstance(select_reaction_list_df.loc[rea.id ,'fluxes in model2'],str):
                            vis_flux2=select_reaction_list_df.loc[rea.id ,'fluxes in model2']
                        else:
                            vis_flux2=str(abs(round(select_reaction_list_df.loc[rea.id ,'fluxes in model2'],7)))                        
                        metlinks[key_id]['flux_in_model2']=vis_flux2
                        metlinks[key_id]['reaction_type']=select_reaction_list_df.loc[rea.id ,'reaction type']                        

        #break
    uni_metabolite_list=list(set(metabolite_list))
    
    return[metlinks,uni_metabolite_list]

def model_calculate(model,json_model_file,use_obj,use_substrate,concentration,ATPM,call_method):    
    #json_model_file="./data/%s_irreversible.json"%model
    #print(model,use_substrate,use_obj)
    GEM_fluxes_outfile = './analysis/%s_%s2%s_solution_df_pfba.tsv'%(model,use_substrate,use_obj)
    norm_model = cobra.io.json.load_json_model(json_model_file)
    norm_model.objective=use_obj
    [ori_obj_id,ori_substrate_id_list,ori_sub_concentration,ori_ATPM]=get_model_substrate_obj(norm_model)
    for eachsubid in ori_substrate_id_list:
        norm_model.reactions.get_by_id(eachsubid).bounds = (0, 0)
    norm_model.reactions.get_by_id(use_substrate).bounds = (-concentration, 0)  
    if ATPM=='No':
        norm_model.reactions.get_by_id('ATPM').bounds = (0, 0)   
    get_fluxes_detail_in_model(norm_model, call_method, GEM_fluxes_outfile) 
    
def get_fluxes_detail_in_model_comparison_old(model1_pfba_solution,model2_pfba_solution, model_comparison_outfile):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """
    try:
        model1_select_reaction_list_df=model1_pfba_solution[abs(model1_pfba_solution['fluxes'])>0]#
        model2_select_reaction_list_df=model2_pfba_solution[abs(model2_pfba_solution['fluxes'])>0]#
    except:
        print('Can not compare solutions in these two models!')
    else:
        model_comparison_pfba_solution_detail = {}
        reaction_intersection=list(set(model1_select_reaction_list_df.index).intersection(set(model2_select_reaction_list_df.index)))
        reaction_only_in_model1=list(set(model1_select_reaction_list_df.index).difference(set(model2_select_reaction_list_df.index)))
        reaction_only_in_model2=list(set(model2_select_reaction_list_df.index).difference(set(model1_select_reaction_list_df.index)))
        reaction_total=(list(set(model1_select_reaction_list_df.index).union(set(model2_select_reaction_list_df.index))))

        for eachreaction in reaction_total:
            model_comparison_pfba_solution_detail[eachreaction]={}
            model_comparison_pfba_solution_detail[eachreaction]['gene in model1']='none'
            model_comparison_pfba_solution_detail[eachreaction]['gene in model2']='none'
            model_comparison_pfba_solution_detail[eachreaction]['equ'] ='none'
            model_comparison_pfba_solution_detail[eachreaction]['fluxes in model1']='none'
            model_comparison_pfba_solution_detail[eachreaction]['fluxes in model2']='none'
            model_comparison_pfba_solution_detail[eachreaction]['reaction type']='none'
            model_comparison_pfba_solution_detail[eachreaction]['ec-code']='none'
            #只在model2
            if eachreaction in reaction_only_in_model2:
                model_comparison_pfba_solution_detail[eachreaction]['gene in model2'] = model2_select_reaction_list_df.loc[eachreaction]['gene id']                    
                model_comparison_pfba_solution_detail[eachreaction]['equ'] = model2_select_reaction_list_df.loc[eachreaction]['equ']    
                model_comparison_pfba_solution_detail[eachreaction]['ec-code'] = model2_select_reaction_list_df.loc[eachreaction]['ec-code'] 
                model_comparison_pfba_solution_detail[eachreaction]['fluxes in model2'] = model2_select_reaction_list_df.loc[eachreaction]['fluxes']            
                model_comparison_pfba_solution_detail[eachreaction]['reaction type']='exsit in model2'
            else:    
                #交集
                if eachreaction in reaction_intersection:
                    model_comparison_pfba_solution_detail[eachreaction]['gene in model1'] = model1_select_reaction_list_df.loc[eachreaction]['gene id']
                    model_comparison_pfba_solution_detail[eachreaction]['gene in model2'] = model2_select_reaction_list_df.loc[eachreaction]['gene id']
                    model_comparison_pfba_solution_detail[eachreaction]['equ'] = model1_select_reaction_list_df.loc[eachreaction]['equ'] 
                    model_comparison_pfba_solution_detail[eachreaction]['fluxes in model1'] = model1_select_reaction_list_df.loc[eachreaction]['fluxes']
                    model_comparison_pfba_solution_detail[eachreaction]['fluxes in model2'] = model2_select_reaction_list_df.loc[eachreaction]['fluxes']            
                    model_comparison_pfba_solution_detail[eachreaction]['ec-code'] = 'mode1: '+model1_select_reaction_list_df.loc[eachreaction]['ec-code']+'; model2: '+model2_select_reaction_list_df.loc[eachreaction]['ec-code']#list(set([model1_select_reaction_list_df.loc[eachreaction]['ec-code'],model2_select_reaction_list_df.loc[eachreaction]['ec-code']]))
                    model_comparison_pfba_solution_detail[eachreaction]['reaction type']='exsit in two model'            
                #只在classic model
                elif eachreaction in reaction_only_in_model1:
                    model_comparison_pfba_solution_detail[eachreaction]['gene in model1'] = model1_select_reaction_list_df.loc[eachreaction]['gene id']                    
                    model_comparison_pfba_solution_detail[eachreaction]['equ'] = model1_select_reaction_list_df.loc[eachreaction]['equ']    
                    model_comparison_pfba_solution_detail[eachreaction]['ec-code'] = model1_select_reaction_list_df.loc[eachreaction]['ec-code'] 
                    model_comparison_pfba_solution_detail[eachreaction]['fluxes in model1'] = model1_select_reaction_list_df.loc[eachreaction]['fluxes']
                    model_comparison_pfba_solution_detail[eachreaction]['reaction type']='exsit in classic model1'
                
        model_comparison_pfba_solution_detail_df = pd.DataFrame(model_comparison_pfba_solution_detail).T
        model_comparison_pfba_solution_detail_df.to_csv(model_comparison_outfile,sep='\t')
        return model_comparison_pfba_solution_detail_df
    
def get_model_met_from_reaction_string_cb(model,reaction_str,met_contain_C,fwd_arrow,rev_arrow,reversible_arrow,term_split):
    # set the arrows
    forward_arrow_finder = (
            _forward_arrow_finder
            if fwd_arrow is None
            else re.compile(re.escape(fwd_arrow))
        )
    reverse_arrow_finder = (
            _reverse_arrow_finder
            if rev_arrow is None
            else re.compile(re.escape(rev_arrow))
        )
    reversible_arrow_finder = (
            _reversible_arrow_finder
            if reversible_arrow is None
            else re.compile(re.escape(reversible_arrow))
        )
    compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
    found_compartments = compartment_finder.findall(reaction_str)
    if len(found_compartments) == 1:
        compartment = found_compartments[0]
        reaction_str = compartment_finder.sub("", reaction_str)
    else:
        compartment = ""
    # reversible case
    arrow_match = reversible_arrow_finder.search(reaction_str)
    if arrow_match is not None:
        pass
    else:  # irreversible
        # try forward
        arrow_match = forward_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            pass
        else:
            # must be reverse
            arrow_match = reverse_arrow_finder.search(reaction_str)
            if arrow_match is None:
                raise ValueError("no suitable arrow found in '%s'" % reaction_str)
            else:
                pass
    reactant_str = reaction_str[: arrow_match.start()].strip()
    product_str = reaction_str[arrow_match.end() :].strip()
    for substr, factor in ((reactant_str, -1), (product_str, 1)):
        if len(substr) == 0:
            continue
        for term in substr.split(term_split):
            term = term.strip()
            if term.lower() == "nothing":
                continue
            if " " in term:
                num_str, met_id = term.split()
                num = float(num_str.lstrip("(").rstrip(")")) * factor
            else:
                met_id = term
                num = factor
            met_id += compartment
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                #必须要加上formula，就是要判断是否含有C
                if met_id in met_contain_C.index:
                    met = Metabolite(met_id,formula=met_contain_C.loc[met_id,'formula'],name=met_contain_C.loc[met_id,'name'])
                else:
                    met = Metabolite(met_id)
                model.add_metabolites(met)
                #reaction.add_metabolites({met: num})
    #model.add_reaction(reaction)
    return model  

def get_model_met_from_reaction_string(model,reaction_str,fwd_arrow,rev_arrow,reversible_arrow,term_split):
    # set the arrows
    forward_arrow_finder = (
            _forward_arrow_finder
            if fwd_arrow is None
            else re.compile(re.escape(fwd_arrow))
        )
    reverse_arrow_finder = (
            _reverse_arrow_finder
            if rev_arrow is None
            else re.compile(re.escape(rev_arrow))
        )
    reversible_arrow_finder = (
            _reversible_arrow_finder
            if reversible_arrow is None
            else re.compile(re.escape(reversible_arrow))
        )
    compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
    found_compartments = compartment_finder.findall(reaction_str)
    if len(found_compartments) == 1:
        compartment = found_compartments[0]
        reaction_str = compartment_finder.sub("", reaction_str)
    else:
        compartment = ""
    # reversible case
    arrow_match = reversible_arrow_finder.search(reaction_str)
    if arrow_match is not None:
        pass
    else:  # irreversible
        # try forward
        arrow_match = forward_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            pass
        else:
            # must be reverse
            arrow_match = reverse_arrow_finder.search(reaction_str)
            if arrow_match is None:
                raise ValueError("no suitable arrow found in '%s'" % reaction_str)
            else:
                pass
    reactant_str = reaction_str[: arrow_match.start()].strip()
    product_str = reaction_str[arrow_match.end() :].strip()
    for substr, factor in ((reactant_str, -1), (product_str, 1)):
        if len(substr) == 0:
            continue
        for term in substr.split(term_split):
            term = term.strip()
            if term.lower() == "nothing":
                continue
            if " " in term:
                num_str, met_id = term.split()
                num = float(num_str.lstrip("(").rstrip(")")) * factor
            else:
                met_id = term
                num = factor
            met_id += compartment
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                met = Metabolite(met_id)
                model.add_metabolites(met)
                #reaction.add_metabolites({met: num})
    #model.add_reaction(reaction)
    return model  

def build_model_from_FBA_result(model_pfba_solution):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    metlink_model = Model('metlink_model')
    for index, row in model_pfba_solution.iterrows():
        metlink_model=get_model_met_from_reaction_string(metlink_model,row['equ'],fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
        reaction = Reaction(index) 
        metlink_model.add_reactions([reaction])
        reaction.build_reaction_from_string(row['equ'])   
    
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
        
    return metlink_model

def build_model_from_FBA_result_cb(model_pfba_solution):
    """由FBA结果构建模型

    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    * met_contain_C: BIGG 模型获得的代谢物列表（分子量）.

    :return: cobra model.
    """
    metlink_model = Model('metlink_model')
    for index, row in model_pfba_solution.iterrows():
        metlink_model=get_model_met_from_reaction_string_cb(metlink_model,row['equ'],met_contain_C,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
        reaction = Reaction(index) 
        metlink_model.add_reactions([reaction])
        reaction.build_reaction_from_string(row['equ'])   
    return metlink_model

def get_metlist_and_metlinklist(model,select_reaction_list_df):
    select_reaction_list=list(select_reaction_list_df.index)
    metabolite_list=[]
    metlinks={} 
    for rea in model.reactions:        
        if rea.id in select_reaction_list:
            subs = [m.id for m in rea.reactants]
            pros = [m.id for m in rea.products]
            sub_pro = []
            for s, p in itertools.product(subs, pros):
                sub_pro.append((s, p))
            for sp in sub_pro:
                metabolite_list.append(sp[0])
                metabolite_list.append(sp[1])
                #flux: mmol gDWh−1 E: mg gDW−1                 
                key_id=rea.id+'_'+sp[0]+'_'+sp[1]
                metlinks[key_id]={}
                if select_reaction_list_df.loc[rea.id ,'fluxes']>0:
                    metlinks[key_id]['source_node']=sp[0]
                    metlinks[key_id]['target_node']=sp[1]
                else:
                    metlinks[key_id]['source_node']=sp[1]
                    metlinks[key_id]['target_node']=sp[0]                    
                metlinks[key_id]['reaction_id']=rea.id
                metlinks[key_id]['reaction_equation']=select_reaction_list_df.loc[rea.id ,'equ']                        
                if isinstance(select_reaction_list_df.loc[rea.id ,'fluxes'],str):
                    vis_flux=select_reaction_list_df.loc[rea.id ,'fluxes']
                else:
                    vis_flux=str(abs(round(select_reaction_list_df.loc[rea.id ,'fluxes'],7)))
                    #vis_flux=str(round(select_reaction_list_df.loc[rea.id ,'fluxes'],3))
                metlinks[key_id]['flux']=vis_flux
    uni_metabolite_list=list(set(metabolite_list))
    return[metlinks,uni_metabolite_list]

def get_model_substrate_obj(use_model):
    ATPM='No' 
    #ATPM_list=['ATPM','rxn00062_c0','DM_atp_c_','mnxr153054c3']
    ATPM_list=['ATPM','rxn00062','rxn00062_c0','DM_atp_c_','mnxr153054c3','R00086','META:ATPASE-RXN','ATPASE-RXN','MNXR96131','rxn09694', 'rxn11300','rxn09694_c0', 'rxn11300_c0','ATPM_c', 'ATPM_LPAREN_NGAM_RPAREN_', 'DM_atp_c', 'DM_atp[c]']
    substrate_list=[]
    concentration_list=[]
    EX_exclude_reaction_list=['EX_pi_e','EX_h_e','EX_fe3_e','EX_mn2_e','EX_co2_e','EX_fe2_e','EX_h2_e','EX_zn2_e',\
                             'EX_mg2_e','EX_ca2_e','EX_so3_e','EX_ni2_e','EX_no_e','EX_cu2_e','EX_hg2_e','EX_cd2_e',\
                             'EX_h2o2_e','EX_h2o_e','EX_no2_e','EX_nh4_e','EX_so4_e','EX_k_e','EX_na1_e','EX_o2_e',\
                             'EX_o2s_e','EX_ag_e','EX_cu_e','EX_so2_e','EX_cl_e','EX_n2o_e','EX_cs1_e','EX_cobalt2_e']
    EX_exclude_reaction_list=EX_exclude_reaction_list+[i+'_reverse' for i in EX_exclude_reaction_list]
    for r in use_model.reactions:
        if r.objective_coefficient == 1:
            obj=r.id #Product name
        #elif not r.lower_bound==0 and not r.lower_bound==-1000 and not r.lower_bound==-999999 and abs(r.lower_bound)>0.1:#排除很小的值
        elif not r.upper_bound==0 and not r.upper_bound==1000 and not r.upper_bound==999999 and abs(r.upper_bound)>0.1:#排除很小的值
            #print(r.id,r.upper_bound,r.lower_bound)
            if r.id in ATPM_list:
                if r.upper_bound>0:
                    ATPM='Yes' #ATP maintenance requirement
            elif r.id not in EX_exclude_reaction_list:
                #print(r.id,r.upper_bound,r.lower_bound)
                #substrate=r.id #Substrate name
                substrate_list.append(r.id)
                #concentration=r.upper_bound #Substrate uptake rate  
                concentration_list.append(r.upper_bound)
    return(obj,substrate_list,concentration_list,ATPM)

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
        if reaction.lower_bound < 0 and reaction.upper_bound == 0:
            for metabolite in reaction.metabolites:
                original_coefficient = reaction.get_coefficient(metabolite)
                reaction.add_metabolites({metabolite: -2*original_coefficient})
            reaction.id += "_reverse"
            reaction.upper_bound = -reaction.lower_bound
            reaction.lower_bound = 0
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
                             for k, v in reaction._metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)
    

def gene_deletion(norm_model, ids, method):
    all_reactions = []
    for g_id in ids:
        all_reactions.extend(
            knock_out_model_genes(norm_model, (norm_model.genes.get_by_id(g_id),))
        )
    with norm_model as model:
        for reaction in all_reactions:
            reaction.knock_out()   
            
        if method=='FBA':
            solution = model.optimize()   
        elif method=='pFBA':
            solution = cobra.flux_analysis.pfba(model)   
        elif method=='MOMA':
            solution = cobra.flux_analysis.moma(model)  
        elif method=='ROOM':
            solution = cobra.flux_analysis.room(model) 
            
        return ([r.id for r in all_reactions], solution.objective_value, solution, model.solver.status)     
    
def reaction_deletion(norm_model, ids, method):
    all_reactions =[norm_model.reactions.get_by_id(r_id) for r_id in ids]
    with norm_model as model:
        for reaction in all_reactions:
            reaction.knock_out()   
            
        if method=='FBA':
            solution = model.optimize()   
        elif method=='pFBA':
            solution = cobra.flux_analysis.pfba(model)   
        elif method=='MOMA':
            solution = cobra.flux_analysis.moma(model)  
        elif method=='ROOM':
            solution = cobra.flux_analysis.room(model) 
        #elif calc_method=='FVA':#出来几千个反应，显示一片有啥用。。
        #    flux_solution = cobra.flux_analysis.variability.flux_variability_analysis(norm_model) 
        #elif calc_method=='Geometric FBA':#慢
        #    flux_solution = cobra.flux_analysis.geometric_fba(norm_model) 
        elif method=='loopless_solution':#Convert an existing solution to a loopless one. CycleFreeFlux: efficient removal of thermodynamically infeasible loops from flux distributions.
            solution = cobra.flux_analysis.loopless.loopless_solution(model)   
        return ([r.id for r in all_reactions], solution.objective_value, solution, model.solver.status)    
    
def get_model_information(model):#Mr.Mao改
    All_rxn = []
    All_met = []
    for met in model.metabolites:
        All_met.append(met.id)
    for rxn in model.reactions:
        All_rxn.append(rxn.id)  # 所有 rxn
    uni_All_met=list(set(All_met))
    uni_All_rxn=list(set(All_rxn))
    return uni_All_met, uni_All_rxn

def get_reaction_uni_met(equation):
    uni_str_line=[]
    if re.search('-->',equation):
        str_line=equation.replace('_e','').replace('_p','').replace('_c','').split(' --> ')
        str_line = [i for i in str_line if i != '']
        uni_str_line=list(set(str_line))
    elif re.search('<--',equation):
        str_line=equation.replace('_e','').replace('_p','').replace('_c','').split(' <-- ')
        str_line = [i for i in str_line if i != '']
        uni_str_line=list(set(str_line))
    elif re.search('<=>',equation):
        str_line=equation.replace('_e','').replace('_p','').replace('_c','').split(' <=> ')
        str_line = [i for i in str_line if i != '']
        uni_str_line=list(set(str_line))

    return uni_str_line

def build_model_from_FBA_result_plus_ET(inimodel,model_pfba_solution,rxn_switch,coordinates,substrate_id,product_id):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    #print(substrate_id,product_id)
    metlink_model = Model2('metlink_model')
    #print(substrate_id,product_id)
    #Demand反应 
    if product_id.startswith('DM_'):
        product_id=product_id.replace('DM_','')
        product_add = inimodel.get_metabolite_by_id(product_id)
        if 'DM_'+product_id not in inimodel.reactions.keys():
            inimodel.add_demand_reaction(product_add.id)
    product_id =product_id.replace('__D_e','_e').replace('__L_e','_e')    
        
    #Sink反应
    if substrate_id.startswith('SK_'):
        substrate_id=substrate_id.replace('SK_','')
        substract_add=inimodel.get_metabolite_by_id(substrate_id)
        inimodel.add_sink_reaction(substract_add.id)
        #substrate_id=sink_reaction.id
    elif re.search('_reverse',substrate_id):
        substrate_id =substrate_id.replace('_reverse','')
    substrate_id =substrate_id.replace('__D_e','_e').replace('__L_e','_e')    
    #    
    substrate_sub=substrate_id.replace('_e','').replace('_p','').replace('_c','')
    product_sub=product_id.replace('_e','').replace('_p','').replace('_c','')
    #print(substrate_id,product_id)
    for index, row in model_pfba_solution.iterrows():
        new_equ=row['equ'].replace('__D_e','_e').replace('__L_e','_e')
        uni_str_line=get_reaction_uni_met(new_equ)    
        if len(uni_str_line)==1 and substrate_sub== uni_str_line[0]:  
            pass
        elif len(uni_str_line)==1 and product_sub== uni_str_line[0]: 
            pass         
        else:
            metlink_model=get_model_met_from_reaction_string(metlink_model,new_equ,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction2(index) 
            try:
                reaction_infor=inimodel.get_reaction_by_id(index)
            except:
                reaction.name=index
                reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(index)+'<hr/>'+\
                        'Reaction subsystem: <br/>'+'none'+'<hr/>'+'Reaction equation: <br/>'+\
                        str(new_equ)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+'none'+'<hr/>'
            else:
                reaction.name=reaction_infor.name
                reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(reaction_infor.name)+'<hr/>'+\
                            'Reaction subsystem: <br/>'+str(reaction_infor.subsystem)+'<hr/>'+'Reaction equation: <br/>'+\
                            str(reaction_infor.build_reaction_string())+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+reaction_infor.gene_reaction_rule+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                #reaction.notes['map_info']['display_name'] = reaction.id + " "+ str("%.2f" %(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f" %(abs(row['fluxes'])))
            elif rxn_switch == 'F':
                #reaction.notes['map_info']['display_name'] = reaction.name + " "+str("%.2f"%(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f"%(abs(row['fluxes'])))
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])

            reaction.notes['map_info']['group']=3
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(new_equ)
             

    #去除H
    H_list=['h_c','h_p','h_e','C00080','C00080[c]', 'C00080[p]','C00080[e]','PROTON', 'MNXM1', 'cpd00067', 'cpd00067_c0','cpd00067_p0','cpd00067_e0', 'HMDB59597', 'GPRLSGONYQIRFK-UHFFFAOYSA-N']
    for each_h in H_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #去除H2O
    H2O_list=['h2o_c','h2o_p','h2o_e','C00001', 'C01328','C00001[c]','C00001[p]','C00001[e]', 'C01328[c]','C01328[p]','C01328[e]', 'CPD-15815', 'HYDROXYL-GROUP', 'WATER', 'MNXM2', 'cpd00001', 'cpd00001_c0', 'cpd00001_p0','cpd00001_e0','cpd15275', 'cpd27222', 'HMDB01039', 'HMDB02111', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    for each_h2o in H2O_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h2o)
        except:
            pass
        else:
            hc_id.remove_from_model()
        #去除Pi
    Pi_list=['pi_c','pi_p','pi_e','C00009', 'C00009[c]','C00009[p]','C00009[e]','CPD-16459', 'PHOSPHATE-GROUP', 'Pi', 'MNXM9', 'cpd00009', 'cpd00009_c0', 'cpd00009_p0','cpd00009_e0','cpd27787', 'HMDB00973', 'HMDB01429', 'HMDB02105', 'HMDB02142', 'HMDB05947', 'NBIIXXVUZAFLBC-UHF']
    for each_pi in Pi_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_pi)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #print(len(metlink_model.reactions.get_by_id('DM_pyr_c').products))  
    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))
    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)
    print(substrate_id,product_id)
    if substrate_id.endswith('_e'):
        tmp_p="_e".join(substrate_id.split('_e')[:-1])+'_p'
        tmp_c="_e".join(substrate_id.split('_e')[:-1])+'_c'
        if substrate_id=='glc__D_e' or substrate_id=='glc_D_e':
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p   
        else:
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p 
 
    if product_id.endswith('_e'):
        tmp_p="_e".join(product_id.split('_e')[:-1])+'_p'
        tmp_c="_e".join(product_id.split('_e')[:-1])+'_c'
        if tmp_c in model_met_list:
            product_id = tmp_c    
        elif tmp_p in model_met_list:
            product_id = tmp_p 
  
    #print(substrate_id,product_id)
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
            
        try:
            inimodel.get_metabolite_by_id(eachmet.id)
        except:
            eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__D_e' 
            try:
                inimodel.get_metabolite_by_id(eachmet_id)
            except:
                eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__L_e' 
                try:
                    inimodel.get_metabolite_by_id(eachmet_id)
                except:
                    pass
                else:
                    metabolite_infor=inimodel.get_metabolite_by_id(eachmet_id)
            else:
                metabolite_infor=inimodel.get_metabolite_by_id(eachmet_id)
        else:
            metabolite_infor=inimodel.get_metabolite_by_id(eachmet.id)
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+str(metabolite_infor.name)+'<hr/>'+\
                        'Metabolite formula: <br/>'+str(metabolite_infor.formula)+'<hr/>'
        eachmet.notes['map_info'] = {}
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A" 
        #print(eachmet.id,substrate_id) 
#         if eachmet.id=='glc__D_e':
#              #print(substrate_id)
#             eachmet.notes['map_info']['substrate'] ='glc__D_e' 
        if eachmet.id==substrate_id:
            #print(substrate_id)
            eachmet.notes['map_info']['substrate'] =substrate_id          
            #print(eachmet.id,substrate_id)    
        if eachmet.id==product_id:
            eachmet.notes['map_info']['product'] =product_id
        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])
        if eachmet.id.endswith('_c') or eachmet.id.endswith('_C') or eachmet.id.endswith('_p') or eachmet.id.endswith('_P') or eachmet.id.endswith('_e') or eachmet.id.endswith('_E'):
            pass
        else:
            eachmet.id=eachmet.id+'_c'
    return metlink_model

def build_model_from_FBA_result_plus(inimodel,model_pfba_solution,rxn_switch,model_change_json,coordinates,substrate_id,product_id):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    #print(substrate_id,product_id)
    metlink_model = Model('metlink_model')
    #print(substrate_id,product_id)
    #Demand反应 
    if product_id.startswith('DM_'):
        product_id=product_id.replace('DM_','')
        product_add = inimodel.metabolites.get_by_id(product_id)
        demand =inimodel.add_boundary(product_add, type='demand')
    product_id =product_id.replace('__D_e','_e').replace('__L_e','_e')    
        
    #Sink反应
    if substrate_id.startswith('SK_'):
        substrate_id=substrate_id.replace('SK_','')
        substract_add=inimodel.metabolites.get_by_id(substrate_id)
        sink_reaction=inimodel.add_boundary(substract_add, type="sink")
        #substrate_id=sink_reaction.id
    elif re.search('_reverse',substrate_id):
        substrate_id =substrate_id.replace('_reverse','')
    substrate_id =substrate_id.replace('__D_e','_e').replace('__L_e','_e')    
    #    
    substrate_sub=substrate_id.replace('_e','').replace('_p','').replace('_c','')
    product_sub=product_id.replace('_e','').replace('_p','').replace('_c','')
    #print(substrate_id,product_id)
    for index, row in model_pfba_solution.iterrows():
        new_equ=row['equ'].replace('__D_e','_e').replace('__L_e','_e')
        uni_str_line=get_reaction_uni_met(new_equ)    
        if len(uni_str_line)==1 and substrate_sub== uni_str_line[0]:  
            pass
        elif len(uni_str_line)==1 and product_sub== uni_str_line[0]: 
            pass         
        else:
            metlink_model=get_model_met_from_reaction_string(metlink_model,new_equ,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            try:
                reaction_infor=inimodel.reactions.get_by_id(index)
            except:
                reaction.name=index
                reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(index)+'<hr/>'+\
                        'Reaction subsystem: <br/>'+'none'+'<hr/>'+'Reaction equation: <br/>'+\
                        str(new_equ)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+'none'+'<hr/>'
            else:
                reaction.name=reaction_infor.name
                reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(reaction_infor.name)+'<hr/>'+\
                            'Reaction subsystem: <br/>'+str(reaction_infor.subsystem)+'<hr/>'+'Reaction equation: <br/>'+\
                            str(reaction_infor.reaction)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+reaction_infor.gene_reaction_rule+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                #reaction.notes['map_info']['display_name'] = reaction.id + " "+ str("%.2f" %(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f" %(abs(row['fluxes'])))
            elif rxn_switch == 'F':
                #reaction.notes['map_info']['display_name'] = reaction.name + " "+str("%.2f"%(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f"%(abs(row['fluxes'])))
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])

            # if index in model_change_json['exchange_infor'].keys():
            #     reaction.notes['map_info']['group']=4
            if index in model_change_json['add_reaction_infor'].keys():
                reaction.notes['map_info']['group']=5       
            elif index in model_change_json['del_reaction_infor'].keys():
                reaction.notes['map_info']['group']=1    
            elif index in model_change_json['del_gene_infor'].keys():
                reaction.notes['map_info']['group']=1 
            else:
                reaction.notes['map_info']['group']=3
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(new_equ) 

    #去除H
    H_list=['h_c','h_p','h_e','C00080','C00080[c]', 'C00080[p]','C00080[e]','PROTON', 'MNXM1', 'cpd00067', 'cpd00067_c0','cpd00067_p0','cpd00067_e0', 'HMDB59597', 'GPRLSGONYQIRFK-UHFFFAOYSA-N']
    for each_h in H_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #去除H2O
    H2O_list=['h2o_c','h2o_p','h2o_e','C00001', 'C01328','C00001[c]','C00001[p]','C00001[e]', 'C01328[c]','C01328[p]','C01328[e]', 'CPD-15815', 'HYDROXYL-GROUP', 'WATER', 'MNXM2', 'cpd00001', 'cpd00001_c0', 'cpd00001_p0','cpd00001_e0','cpd15275', 'cpd27222', 'HMDB01039', 'HMDB02111', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    for each_h2o in H2O_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h2o)
        except:
            pass
        else:
            hc_id.remove_from_model()
        #去除Pi
    Pi_list=['pi_c','pi_p','pi_e','C00009', 'C00009[c]','C00009[p]','C00009[e]','CPD-16459', 'PHOSPHATE-GROUP', 'Pi', 'MNXM9', 'cpd00009', 'cpd00009_c0', 'cpd00009_p0','cpd00009_e0','cpd27787', 'HMDB00973', 'HMDB01429', 'HMDB02105', 'HMDB02142', 'HMDB05947', 'NBIIXXVUZAFLBC-UHF']
    for each_pi in Pi_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_pi)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #print(len(metlink_model.reactions.get_by_id('DM_pyr_c').products))  
    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))
    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)
    print(substrate_id,product_id)
    if substrate_id.endswith('_e'):
        tmp_p="_e".join(substrate_id.split('_e')[:-1])+'_p'
        tmp_c="_e".join(substrate_id.split('_e')[:-1])+'_c'
        if substrate_id=='glc__D_e' or substrate_id=='glc_D_e':
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p   
        else:
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p 
 
    if product_id.endswith('_e'):
        tmp_p="_e".join(product_id.split('_e')[:-1])+'_p'
        tmp_c="_e".join(product_id.split('_e')[:-1])+'_c'
        if tmp_c in model_met_list:
            product_id = tmp_c    
        elif tmp_p in model_met_list:
            product_id = tmp_p 
  
    #print(substrate_id,product_id)
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
            
        try:
            inimodel.metabolites.get_by_id(eachmet.id)
        except:
            eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__D_e' 
            try:
                inimodel.metabolites.get_by_id(eachmet_id)
            except:
                eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__L_e' 
                try:
                    inimodel.metabolites.get_by_id(eachmet_id)
                except:
                    pass
                else:
                    metabolite_infor=inimodel.metabolites.get_by_id(eachmet_id)
            else:
                metabolite_infor=inimodel.metabolites.get_by_id(eachmet_id)
        else:
            metabolite_infor=inimodel.metabolites.get_by_id(eachmet.id)
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+str(metabolite_infor.name)+'<hr/>'+\
                        'Metabolite formula: <br/>'+str(metabolite_infor.formula)+'<hr/>'
        eachmet.notes['map_info'] = {}
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A" 
        #print(eachmet.id,substrate_id) 
#         if eachmet.id=='glc__D_e':
#              #print(substrate_id)
#             eachmet.notes['map_info']['substrate'] ='glc__D_e' 
        if eachmet.id==substrate_id:
            #print(substrate_id)
            eachmet.notes['map_info']['substrate'] =substrate_id          
            #print(eachmet.id,substrate_id)    
        if eachmet.id==product_id:
            eachmet.notes['map_info']['product'] =product_id
        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])
        if eachmet.id.endswith('_c') or eachmet.id.endswith('_C') or eachmet.id.endswith('_p') or eachmet.id.endswith('_P') or eachmet.id.endswith('_e') or eachmet.id.endswith('_E'):
            pass
        else:
            eachmet.id=eachmet.id+'_c'
    return metlink_model

def build_model_from_FBA_result_plusold(inimodel,model_pfba_solution,rxn_switch,model_change_json,coordinates,substrate_id,product_id):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    #print(substrate_id,product_id)
    metlink_model = Model('metlink_model')
    
    #Demand反应
    if product_id.startswith('DM_'):
        product_id=product_id.replace('DM_','')
        product_add = inimodel.metabolites.get_by_id(product_id)
        demand =inimodel.add_boundary(product_add, type='demand')
    for index, row in model_pfba_solution.iterrows():
        if re.search('_e ',row['equ']):
            pass   
        elif re.search('_p ',row['equ']) and not re.search('glc__D_p',row['equ']):
            pass        
        else:
            #print(row['equ'])
            metlink_model=get_model_met_from_reaction_string(metlink_model,row['equ'],fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            reaction_infor=inimodel.reactions.get_by_id(index)
            reaction.name=reaction_infor.name
            reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(reaction_infor.name)+'<hr/>'+\
                        'Reaction subsystem: <br/>'+str(reaction_infor.subsystem)+'<hr/>'+'Reaction equation: <br/>'+\
                        str(reaction_infor.reaction)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+reaction_infor.gene_reaction_rule+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                #reaction.notes['map_info']['display_name'] = reaction.id + " "+ str("%.2f" %(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f" %(abs(row['fluxes'])))
            elif rxn_switch == 'F':
                #reaction.notes['map_info']['display_name'] = reaction.name + " "+str("%.2f"%(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f"%(abs(row['fluxes'])))
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])

            # if index in model_change_json['exchange_infor'].keys():
            #     reaction.notes['map_info']['group']=4
            if index in model_change_json['add_reaction_infor'].keys():
                reaction.notes['map_info']['group']=5       
            elif index in model_change_json['del_reaction_infor'].keys():
                reaction.notes['map_info']['group']=1    
            elif index in model_change_json['del_gene_infor'].keys():
                reaction.notes['map_info']['group']=1 
            else:
                reaction.notes['map_info']['group']=3
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(row['equ'])  
        
    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))
    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)
    #print(substrate_id,product_id)
    if substrate_id.endswith('_e'):
        if substrate_id=='glc__D_e' or substrate_id=='glc_D_e':
            tmp_p="_e".join(substrate_id.split('_e')[:-1])+'_p'
            tmp_c="_e".join(substrate_id.split('_e')[:-1])+'_c'
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p   
        else:
            substrate_id = "_e".join(substrate_id.split('_e')[:-1])+'_c'
    if product_id.endswith('_e'):
        product_id = "_e".join(product_id.split('_e')[:-1])+'_c'    
    #print(substrate_id,product_id)
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
        metabolite_infor=inimodel.metabolites.get_by_id(eachmet.id)
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+str(metabolite_infor.name)+'<hr/>'+\
                        'Metabolite formula: <br/>'+str(metabolite_infor.formula)+'<hr/>'
        eachmet.notes['map_info'] = {}
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A" 
        if eachmet.id==substrate_id:
            #print(substrate_id)
            eachmet.notes['map_info']['substrate'] =substrate_id
        if eachmet.id==product_id:
            eachmet.notes['map_info']['product'] =product_id
        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])

    return metlink_model

#模型的顺序要与模型的名字相同
##by yangchunhe
def get_fluxes_detail_in_model_comparison(model_pfba_solution,model_comparison_outfile='./'):
    modelCount=0
    
#     model_name=list(model_pfba_solution.keys())
    model_name= dict(zip(list(range(1,len(model_pfba_solution.keys())+1)),list(model_pfba_solution.keys())))
    model_pfba_solution = list(model_pfba_solution.values())
   
     #模型列表
    model_select_reaction_list=[]
    for index in range(len(model_pfba_solution)):
        modelCount = modelCount + 1
        model_select_reaction_list.append(model_pfba_solution[index][abs(model_pfba_solution[index]['fluxes'])>0])

    model_comparison_pfba_solution_detail = {}
    #先合并，再去重
    #总列表
    total_model_index_list = []       
    for i in range(modelCount):
        total_model_index_list.extend(list(model_select_reaction_list[i].index))
    #去重得到并集
    reaction_total = list(set(total_model_index_list))

    #给反应分类
    eachreaction_dict={}
    for eachreaction in reaction_total: 
        li = ''
        for i in range(len(model_select_reaction_list)):
            if eachreaction in (model_select_reaction_list[i].index):
                li = li + f'{i+1}'+' '
            else:
                continue
        eachreaction_dict[eachreaction]=li
    #统计分类
    flipped = {} 
    for key, value in eachreaction_dict.items(): 
        if value not in flipped: 
            flipped[value] = [key] 
        else: 
            flipped[value].append(key)      
    model_comparison_pfba_solution_detail = {}
   
    for key, value in flipped.items():
        models = key.strip().split()
        reactions = list(value)
        for reaction_item in reactions:
            model_comparison_pfba_solution_detail[reaction_item]={}
#             ecode ={}
            fluxes={}
            gene={}
            for item in models:
                model_comparison_pfba_solution_detail[reaction_item][f'flux_in_model{item}'] = model_select_reaction_list[int(item)-1].loc[reaction_item]['fluxes']
                gene[f'{model_name[int(item)]}']=model_select_reaction_list[int(item)-1].loc[reaction_item]['gene id']
    #           ecode[f'{model_name[i]}']=model_select_reaction_list[int(models[i])-1].loc[reaction_item]['ec-code']
                fluxes[f'{model_name[int(item)]}']=model_select_reaction_list[int(item)-1].loc[reaction_item]['fluxes']
  
#             if len(ecode)>1:
#                 st=''
#                 for k,v in ecode.items():
#                     st=st+k+":"+v+'<hr/>'
#                 model_comparison_pfba_solution_detail[reaction_item]['ec-code'] = st.replace("'",'')
            model_comparison_pfba_solution_detail[reaction_item]['equ'] = model_select_reaction_list[int(item)-1].loc[reaction_item]['equ']  
            model_comparison_pfba_solution_detail[reaction_item]['reaction_name'] = model_select_reaction_list[int(models[0])-1].loc[reaction_item]['reaction_name'] 
            model_comparison_pfba_solution_detail[reaction_item]['gene'] = str(gene).replace('{',"").replace('}','').replace("'",'')
            model_comparison_pfba_solution_detail[reaction_item]['fluxes'] = str(fluxes).replace('{',"").replace('}','').replace("'",'')
            if len(key)>2:
                kk = []
                for i in range(len(key)):
                    if key[i] != ' ':
#                         kk =kk + ' model' + key[i]
                        kk.append(model_name[int(key[i])])
                kk = (' and ').join(kk)   
                model_comparison_pfba_solution_detail[reaction_item]['reaction type']=f'exsit in {kk}'
            else:
                model_comparison_pfba_solution_detail[reaction_item]['reaction type']=f'exsit in {model_name[int(key)]}'
#                 model_comparison_pfba_solution_detail[reaction_item]['ec-code'] = str(ecode).replace('{',"").replace('}','').replace("'",'').replace('"','')
            
    model_comparison_pfba_solution_detail_df = pd.DataFrame(model_comparison_pfba_solution_detail).T
    model_comparison_pfba_solution_detail_df['gene']=model_comparison_pfba_solution_detail_df['gene'].str.replace(',',';')
    model_comparison_pfba_solution_detail_df['fluxes']=model_comparison_pfba_solution_detail_df['fluxes'].str.replace(',',';')
    model_comparison_pfba_solution_detail_df['reaction_id']=list(model_comparison_pfba_solution_detail_df.index)  
    model_comparison_pfba_solution_detail_df.reset_index(drop=True,inplace=True)
    column_order = ['reaction_id', 'reaction_name', 'equ','gene','fluxes','reaction type']
    model_comparison_pfba_solution_detail_df = model_comparison_pfba_solution_detail_df[column_order + [col for col in model_comparison_pfba_solution_detail_df.columns if col not in column_order]]
    model_comparison_pfba_solution_detail_df.to_csv(model_comparison_outfile,index=False,sep='\t')  
    return model_comparison_pfba_solution_detail_df   

def build_model_from_FBA_result_cb_plus(model_comparison_results,typelist,rxn_switch,bigg_models_metabolites,coordinates,substrate_id_list,product_id_list):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    #print(substrate_id_list,product_id_list)
    metlink_model = Model('metlink_model')
    new_substrate_id_list=[]
    new_product_id_list=[]
    substrate_sub_list=[]
    product_sub_list=[]
    for eachproduct in product_id_list:
        #Demand反应
        eachproduct =eachproduct.replace('__D_e','_e').replace('__L_e','_e')
        if eachproduct.startswith('DM_'):
            eachproduct=eachproduct.replace('DM_','')
        product_sub=eachproduct.replace('_e','').replace('_p','').replace('_c','')
        product_sub_list.append(product_sub)
        new_product_id_list.append(eachproduct)
        
    for eachsubstrate in substrate_id_list:   
        eachsubstrate =eachsubstrate.replace('__D_e','_e').replace('__L_e','_e')
        #Sink
        if eachsubstrate.startswith('SK_'):
            eachsubstrate=eachsubstrate.replace('SK_','')
        elif re.search('_reverse',eachsubstrate):
            eachsubstrate =eachsubstrate.replace('_reverse','')
        substrate_sub=eachsubstrate.replace('_e','').replace('_p','').replace('_c','')
        new_substrate_id_list.append(eachsubstrate)
        substrate_sub_list.append(substrate_sub)

    model_comparison_results=model_comparison_results.fillna('')
    #print(model_comparison_results)
    for index, row in model_comparison_results.iterrows():
        new_equ=row['equ'].replace('__D_e','_e').replace('__L_e','_e')
        #print(substrate_id,product_id)
        uni_str_line=get_reaction_uni_met(new_equ)    
        if len(uni_str_line)==1 and uni_str_line[0] in substrate_sub_list:  
            pass
        elif len(uni_str_line)==1 and uni_str_line[0] in product_sub_list:  
            pass        
        else:   
            metlink_model=get_model_met_from_reaction_string(metlink_model,new_equ,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            reaction.name=row['reaction_name']
            reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(row['reaction_name'])+'<hr/>'+\
                        'Reaction fluxes: <br/>'+str(row['fluxes'])+'<hr/>'+'Reaction equation: <br/>'+\
                        str(row['equ'])+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+str(row['gene'])+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                reaction.notes['map_info']['display_name'] = index
            elif rxn_switch == 'F':
                reaction.notes['map_info']['display_name'] = row['reaction_name']

            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            for i in range(len(typelist)):           
                if row['reaction type'] ==typelist[i]:
                    reaction.notes['map_info']['group']=i+1
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(new_equ)  
   #去除H
    H_list=['h_c','h_p','h_e','C00080','C00080[c]', 'C00080[p]','C00080[e]','PROTON', 'MNXM1', 'cpd00067', 'cpd00067_c0','cpd00067_p0','cpd00067_e0', 'HMDB59597', 'GPRLSGONYQIRFK-UHFFFAOYSA-N']
    for each_h in H_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #去除H2O
    H2O_list=['h2o_c','h2o_p','h2o_e','C00001', 'C01328','C00001[c]','C00001[p]','C00001[e]', 'C01328[c]','C01328[p]','C01328[e]', 'CPD-15815', 'HYDROXYL-GROUP', 'WATER', 'MNXM2', 'cpd00001', 'cpd00001_c0', 'cpd00001_p0','cpd00001_e0','cpd15275', 'cpd27222', 'HMDB01039', 'HMDB02111', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    for each_h2o in H2O_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h2o)
        except:
            pass
        else:
            hc_id.remove_from_model()
        #去除Pi
    Pi_list=['pi_c','pi_p','pi_e','C00009', 'C00009[c]','C00009[p]','C00009[e]','CPD-16459', 'PHOSPHATE-GROUP', 'Pi', 'MNXM9', 'cpd00009', 'cpd00009_c0', 'cpd00009_p0','cpd00009_e0','cpd27787', 'HMDB00973', 'HMDB01429', 'HMDB02105', 'HMDB02142', 'HMDB05947', 'NBIIXXVUZAFLBC-UHF']
    for each_pi in Pi_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_pi)
        except:
            pass
        else:
            hc_id.remove_from_model()

    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)           
    #print(substrate_id_list,product_id_list) 

    for i in range(len(substrate_id_list)):
        if substrate_id_list[i].endswith('_e'):
            tmp_p="_e".join(substrate_id_list[i].split('_e')[:-1])+'_p'
            tmp_c="_e".join(substrate_id_list[i].split('_e')[:-1])+'_c'
            if substrate_id_list[i]=='glc__D_e' or substrate_id_list[i]=='glc_D_e':
                if tmp_c in model_met_list:
                    substrate_id_list[i] = tmp_c    
                elif tmp_p in model_met_list:
                    substrate_id_list[i] = tmp_p   
            else:
                if tmp_c in model_met_list:
                    substrate_id_list[i] = tmp_c    
                elif tmp_p in model_met_list:
                    substrate_id_list[i] = tmp_p 

    for j in range(len(product_id_list)):                
        if product_id_list[j].endswith('_e'):
            tmp_p="_e".join(product_id_list[j].split('_e')[:-1])+'_p'
            tmp_c="_e".join(product_id_list[j].split('_e')[:-1])+'_c'
            if tmp_c in model_met_list:
                product_id_list[j] = tmp_c    
            elif tmp_p in model_met_list:
                product_id_list[j] = tmp_p 

    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))    
    #bigg_models_metabolites=bigg_models_metabolites.fillna('')
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
        eachmet.notes['map_info'] = {}
        show_met=''
        if eachmet.id in list(bigg_models_metabolites.index):
            show_met=bigg_models_metabolites.loc[eachmet.id,'name']
        else:
            new_met="_e".join(eachmet.id.split('_e')[:-1])+'__D_e'  
            if new_met in list(bigg_models_metabolites.index):
                show_met=bigg_models_metabolites.loc[new_met,'name']
            else:
                new_met="_e".join(eachmet.id.split('_e')[:-1])+'__L_e' 
                if new_met in list(bigg_models_metabolites.index):
                    show_met=bigg_models_metabolites.loc[new_met,'name']
                else:
                    print(eachmet.id)
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+str(show_met)+'<hr/>'
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A"  
        #print(eachmet.id) 
        if eachmet.id in new_substrate_id_list:
            eachmet.notes['map_info']['substrate'] =eachmet.id
        if eachmet.id in new_product_id_list:
            eachmet.notes['map_info']['product'] =eachmet.id
 
        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])
        if eachmet.id.endswith('_c') or eachmet.id.endswith('_C') or eachmet.id.endswith('_p') or eachmet.id.endswith('_P') or eachmet.id.endswith('_e') or eachmet.id.endswith('_E'):
            pass
        else:
            eachmet.id=eachmet.id+'_c'
    return metlink_model

def build_model_from_FBA_result_cb_plusold(model_comparison_results,typelist,rxn_switch,bigg_models_metabolites,coordinates,substrate_id_list,product_id_list):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    metlink_model = Model('metlink_model')
    for index, row in model_comparison_results.iterrows():
        if re.search('_e ',row['equ']):
            pass   
        elif re.search('_p ',row['equ']) and not re.search('glc__D_p',row['equ']):
            pass        
        else:
            metlink_model=get_model_met_from_reaction_string(metlink_model,row['equ'],fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            reaction.name=row['reaction_name']
            reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(row['reaction_name'])+'<hr/>'+\
                        'Reaction fluxes: <br/>'+str(row['fluxes'])+'<hr/>'+'Reaction equation: <br/>'+\
                        str(row['equ'])+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+str(row['gene'])+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                reaction.notes['map_info']['display_name'] = index
            elif rxn_switch == 'F':
                reaction.notes['map_info']['display_name'] = row['reaction_name']

            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            for i in range(len(typelist)):           
                if row['reaction type'] ==typelist[i]:
                    reaction.notes['map_info']['group']=i+1
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(row['equ'])  
    #print(substrate_id_list,product_id_list)
    for i in range(len(substrate_id_list)):
        if substrate_id_list[i].endswith('_e'):
            if substrate_id_list[i]=='glc__D_e':
                substrate_id_list[i] = "_e".join(substrate_id_list[i].split('_e')[:-1])+'_p'    
            else:
                substrate_id_list[i] = "_e".join(substrate_id_list[i].split('_e')[:-1])+'_c'
    for j in range(len(product_id_list)):                
        if product_id_list[j].endswith('_e'):
            product_id_list[j] = "_e".join(product_id_list[j].split('_e')[:-1])+'_c'
    #print(substrate_id_list,product_id_list)    
    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))    
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
        eachmet.notes['map_info'] = {}
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+\
                        'Metabolite name: <br/>'+str(bigg_models_metabolites.loc[eachmet.id,'name'])+'<hr/>'
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        if eachmet.id in ["atp_c", "nadph_c",'nadh_c', 'nadp_c', 'nad_c','adp_c']:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A"  
        if eachmet.id in substrate_id_list:
            eachmet.notes['map_info']['substrate'] =eachmet.id
        if eachmet.id in product_id_list:
            eachmet.notes['map_info']['product'] =eachmet.id
        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])
    return metlink_model

def reduction_degree(met):
    rd_element={'C':4,'H':1,'O':-2,'N':-3,'S':6,'P':5}
    #calculate the reduction degree for a metabolite object, special process for coa, ACP containing molecules
    me=met.elements
    metid=met.id
    rd=0
    for element, count in me.items():
        if element in rd_element:
            rd+=count *rd_element[element]
    rd-=met.charge #need to correct with charge
    return rd

def met_generation_reaction(model,reaction_id,met_list,met_name,flux_solution_df,result_static):
    stat_name='%s generation reaction'%met_name
    for eachm in model.get_reaction_by_id(reaction_id).metabolites.keys():
        if model.get_reaction_by_id(reaction_id).metabolites[eachm]>0:
            if eachm in met_list and flux_solution_df.loc[reaction_id,'fluxes']>0:
                result_static[stat_name][reaction_id]={}
                result_static[stat_name][reaction_id]['reaction_id']=reaction_id
                result_static[stat_name][reaction_id]['reaction_name']=model.get_reaction_by_id(reaction_id).name
                result_static[stat_name][reaction_id]['fluxes']=abs(round(flux_solution_df.loc[reaction_id,'fluxes'],7))
                result_static[stat_name][reaction_id]['gene id']=model.get_reaction_by_id(reaction_id).gene_reaction_rule
                result_static[stat_name][reaction_id]['equ']=model.get_reaction_by_id(reaction_id).build_reaction_string()
                result_static[stat_name][reaction_id]['equ_name']=model.get_reaction_by_id(reaction_id).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
                break
        elif model.get_reaction_by_id(reaction_id).metabolites[eachm]<0:
            if eachm in met_list and flux_solution_df.loc[reaction_id,'fluxes']<0:
                result_static[stat_name][reaction_id]={}
                result_static[stat_name][reaction_id]['reaction_id']=reaction_id
                result_static[stat_name][reaction_id]['reaction_name']=model.get_reaction_by_id(reaction_id).name
                result_static[stat_name][reaction_id]['fluxes']=abs(round(flux_solution_df.loc[reaction_id,'fluxes'],7))
                result_static[stat_name][reaction_id]['gene id']=model.get_reaction_by_id(reaction_id).gene_reaction_rule
                result_static[stat_name][reaction_id]['equ']=model.get_reaction_by_id(reaction_id).build_reaction_string()
                result_static[stat_name][reaction_id]['equ_name']=model.get_reaction_by_id(reaction_id).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
                break
         
    #return result_static

def met_consumption_reaction(model,reaction_id,met_list,met_name,flux_solution_df,result_static):
    stat_name='%s consumption reaction'%met_name
    for eachm in model.get_reaction_by_id(reaction_id).metabolites.keys():
        if model.get_reaction_by_id(reaction_id).metabolites[eachm]>0:
            if eachm in met_list and flux_solution_df.loc[reaction_id,'fluxes']>0:
                result_static[stat_name][reaction_id]={}
                result_static[stat_name][reaction_id]['reaction_id']=reaction_id
                result_static[stat_name][reaction_id]['reaction_name']=model.get_reaction_by_id(reaction_id).name
                result_static[stat_name][reaction_id]['fluxes']=abs(round(flux_solution_df.loc[reaction_id,'fluxes'],7))
                result_static[stat_name][reaction_id]['gene id']=model.get_reaction_by_id(reaction_id).gene_reaction_rule
                result_static[stat_name][reaction_id]['equ']=model.get_reaction_by_id(reaction_id).build_reaction_string()
                result_static[stat_name][reaction_id]['equ_name']=model.get_reaction_by_id(reaction_id).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
                break
        elif model.get_reaction_by_id(reaction_id).metabolites[eachm]<0:
            if eachm in met_list and flux_solution_df.loc[reaction_id,'fluxes']<0:
                result_static[stat_name][reaction_id]={}
                result_static[stat_name][reaction_id]['reaction_id']=reaction_id
                result_static[stat_name][reaction_id]['reaction_name']=model.get_reaction_by_id(reaction_id).name
                result_static[stat_name][reaction_id]['fluxes']=abs(round(flux_solution_df.loc[reaction_id,'fluxes'],7))
                result_static[stat_name][reaction_id]['gene id']=model.get_reaction_by_id(reaction_id).gene_reaction_rule
                result_static[stat_name][reaction_id]['equ']=model.get_reaction_by_id(reaction_id).build_reaction_string()
                result_static[stat_name][reaction_id]['equ_name']=model.get_reaction_by_id(reaction_id).build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
                break

    #return result_static
 
    
def get_composition_from_formula(met_formula):
    element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")
    composition = {}
    parsed = element_re.findall(met_formula)
    for (element, count) in parsed:
        if count == "":
            count = 1
        else:
            try:
                count = float(count)
                int_count = int(count)
                if count == int_count:
                    count = int_count
                else:
                    warn(
                            "%s is not an integer (in formula %s)"
                            % (count, self.formula)
                        )
            except ValueError:
                warn("failed to parse %s (in formula %s)" % (count, self.formula))
                self.elements = {}
        if element in composition:
            composition[element] += count
        else:
            composition[element] = count
    return (composition)



def get_preview_default_jsonold(norm_model,met_C):

    model_preview_json={}
    default_setting_json={}

    model_preview_json['Reactions']=len(norm_model.reactions)
    model_preview_json['Metabolites']=len(norm_model.metabolites)
    model_preview_json['Genes']=len(norm_model.genes)
    
    default_setting_json['Oxygen']={}
    default_setting_json['Oxygen']['model_id']='None'
    default_setting_json['Oxygen']['model_correspondind']='None'
    default_setting_json['Oxygen']['value']=0 
    default_setting_json['Carbon_source']={}
    default_setting_json['Carbon_source']['substrate_id']='None'
    default_setting_json['Carbon_source']['substrate_model_correspondind']='None'
    default_setting_json['Carbon_source']['substrate_name']='None'
    default_setting_json['Carbon_source']['uptake_rate']=0 
    default_setting_json['Nitrogen_source']={}
    default_setting_json['Nitrogen_source']['substrate_id']='None'
    default_setting_json['Nitrogen_source']['substrate_model_correspondind']='None'
    default_setting_json['Nitrogen_source']['substrate_name']='None' 
    default_setting_json['Nitrogen_source']['uptake_rate']=0 
    default_setting_json['ATPM']={}
    default_setting_json['ATPM']['model_id']='None'
    default_setting_json['ATPM']['model_correspondind']='None'
    default_setting_json['ATPM']['value']=0 

    O2_list=['EX_o2_e_reverse','EX_o2_e','EX_cpd00007_b','EX_cpd00007_e0','EX_o2(e)','EX_MNXM735438@BOUNDARY','mnxr102090c2b','Ex_o2','Ex_o2_reverse','R1287']
    #ATPM_list=['ATPM','rxn00062_c0','DM_atp_c_','mnxr153054c3']
    ATPM_list=['ATPM','rxn00062','rxn00062_c0','DM_atp_c_','mnxr153054c3','R00086','META:ATPASE-RXN','ATPASE-RXN','MNXR96131','rxn09694', 'rxn11300','rxn09694_c0', 'rxn11300_c0','ATPM_c', 'ATPM_LPAREN_NGAM_RPAREN_', 'DM_atp_c', 'DM_atp[c]']

    #输出培养基
    medium_name_list=[]
    carbon_source_list=[]
    nitrogen_source_list=[]
    try:
        find_boundary_types(norm_model,'exchange')
    except:
        print('No medium information!')
    else:       
        for eachm in medium(norm_model):
            is_C_met='No'
            is_N_met='No'
            medium_name=norm_model.reactions.get_by_id(eachm).name.replace(' exchange','').replace('exchange','').replace('EXCHANGE','')
            if len(medium_name)==0:
                medium_name=eachm
            medium_name_list.append(medium_name)
            try:
                reactant_id=norm_model.reactions.get_by_id(eachm).reactants[0].id
            except:
                reactant_id=norm_model.reactions.get_by_id(eachm).products[0].id   
            try:
                product_C_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['C']
            except:
                try:
                    product_C_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['C']
                except:
                    print('No formula for %s!'%reactant_id)
                else: 
                    is_C_met='Yes'
            else:
                is_C_met='Yes'

            if is_C_met=='Yes':
                if abs(norm_model.reactions.get_by_id(eachm).lower_bound) <1000 and abs(norm_model.reactions.get_by_id(eachm).lower_bound) >0:
                    if 'Carbon_source' in default_setting_json.keys():
                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).lower_bound):
                            default_setting_json['Carbon_source']['substrate_id']=reactant_id
                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                            default_setting_json['Carbon_source']['substrate_name']=medium_name
                            default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)             
                    else:
                        default_setting_json['Carbon_source']={}
                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                        default_setting_json['Carbon_source']['substrate_id']=reactant_id
                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                        default_setting_json['Carbon_source']['substrate_name']=medium_name
                        default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                elif norm_model.reactions.get_by_id(eachm).upper_bound>0 and norm_model.reactions.get_by_id(eachm).upper_bound<1000:
                    if 'Carbon_source' in default_setting_json.keys():
                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).upper_bound):
                            default_setting_json['Carbon_source']['substrate_id']=reactant_id
                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                            default_setting_json['Carbon_source']['substrate_name']=medium_name
                            default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).upper_bound)             
                    else:
                        default_setting_json['Carbon_source']={}
                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                        default_setting_json['Carbon_source']['substrate_id']=reactant_id
                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                        default_setting_json['Carbon_source']['substrate_name']=medium_name
                        default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).upper_bound)                
                carbon_source_list.append(medium_name)
            try:
                product_N_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['N']
            except:
                try:
                    product_N_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['N']
                except:
                    print('No formula for %s!'%reactant_id)
                    if reactant_id=='nh4_e':
                        default_setting_json['Nitrogen_source']={}
                        #default_setting_json['Nitrogen_source']['substrate_model_correspondind']='EX_'+reactant_id
                        default_setting_json['Nitrogen_source']['substrate_id']=reactant_id
                        default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachm
                        default_setting_json['Nitrogen_source']['substrate_name']=medium_name
                        default_setting_json['Nitrogen_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                        nitrogen_source_list.append(medium_name)  
                else: 
                    is_N_met='Yes'
            else:
                is_N_met='Yes'

            if is_N_met=='Yes':
                if reactant_id=='nh3_e':
                    default_setting_json['Nitrogen_source']={}
                    #default_setting_json['Nitrogen_source']['substrate_model_correspondind']='EX_'+reactant_id
                    default_setting_json['Nitrogen_source']['substrate_id']=reactant_id
                    default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachm
                    default_setting_json['Nitrogen_source']['substrate_name']=medium_name
                    default_setting_json['Nitrogen_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                    nitrogen_source_list.append(medium_name)  
                else:
                    nitrogen_source_list.append(medium_name) 
    #ATPM
    for eachatpm in ATPM_list:
        try:                
            ATPM=norm_model.reactions.get_by_id(eachatpm).lower_bound
            print(eachatpm,ATPM)
        except:
            print('Model does not contain ATP Maintenance Requirement')
            #default_setting_json['ATPM']='none'
        else:
            default_setting_json['ATPM']={}
            default_setting_json['ATPM']['model_id']=eachatpm
            default_setting_json['ATPM']['model_correspondind']=eachatpm
            default_setting_json['ATPM']['value']=ATPM 
            break

    #有氧厌氧
    default_setting_json['Oxygen']={}
    default_setting_json['Oxygen']['model_id']='None'
    default_setting_json['Oxygen']['model_correspondind']='None'
    default_setting_json['Oxygen']['value']=0 
    for eachatpm in O2_list:
        try:
            O2=norm_model.reactions.get_by_id(eachatpm).lower_bound
        except:
            print('No O2 reaction!')         
        else:
            default_setting_json['Oxygen']={}
            default_setting_json['Oxygen']['model_id']=eachatpm
            default_setting_json['Oxygen']['model_correspondind']=eachatpm
            default_setting_json['Oxygen']['value']=abs(O2)
            break

    model_preview_json['Medium']=medium_name_list
    model_preview_json['Carbon_source']=carbon_source_list
    model_preview_json['Nitrogen_source']=nitrogen_source_list    
    #输出目标函数
    if len(str(norm_model.objective.expression))>300:
        model_preview_json['Optimization Objectives']='Biomass' 
        regex = re.compile('BIOMASS', flags=re.IGNORECASE)
        if len(norm_model.reactions.query(regex))>1:
            for eglc in norm_model.reactions.query(regex):
                with norm_model as model:
                    model.objective=eglc.id
                    try:
                        opt_r=model.optimize()
                    except:
                        pass
                    else:
                        if opt_r[eglc.id]>0:
                            default_setting_json['Object']={}
                            default_setting_json['Object']['model_id']=eglc.id.replace('EX_','')
                            default_setting_json['Object']['model_correspondind']=eglc.id
                            default_setting_json['Object']['substrate_name']='Biomass_core'
                            #default_setting_json['Object']['model_reaction']=eglc.reaction
                            model_preview_json['Optimization Objectives']=eglc.name
        else: 
            for eglc in norm_model.reactions.query(regex):
                default_setting_json['Object']={}
                default_setting_json['Object']['model_id']=eglc.id.replace('EX_','')
                default_setting_json['Object']['model_correspondind']=eglc.id
                default_setting_json['Object']['substrate_name']='Biomass'
                #default_setting_json['Object']['model_reaction']=eglc.reaction
                model_preview_json['Optimization Objectives']=eglc.name
    elif str(norm_model.objective.expression)=='0':
        model_preview_json['Optimization Objectives']='Not set' 
        default_setting_json['Object']={}
        default_setting_json['Object']['model_id']='None'
        default_setting_json['Object']['model_correspondind']='None'
        #default_setting_json['Object']['model_reaction']='None'
        default_setting_json['Object']['substrate_name']='None'
    else:
        obj_id=str(norm_model.objective.expression).split(' -')[0].split('*')[1]
        obj_name=norm_model.reactions.get_by_id(obj_id).name.replace(' exchange','')
        model_preview_json['Optimization Objectives']=obj_name 
        if re.search('BIOMASS',obj_id, re.I):
            default_setting_json['Object']={}
            default_setting_json['Object']['model_id']=obj_id.replace('EX_','')
            default_setting_json['Object']['model_correspondind']=obj_id 
            default_setting_json['Object']['substrate_name']='Biomass'
            #default_setting_json['Object']['model_reaction']=norm_model.reactions.get_by_id(obj_id).reaction
        else:
            default_setting_json['Object']={}
            default_setting_json['Object']['model_id']=obj_id.replace('EX_','') 
            default_setting_json['Object']['model_correspondind']=obj_id
            default_setting_json['Object']['substrate_name']=obj_name
            #default_setting_json['Object']['model_reaction']=norm_model.reactions.get_by_id(obj_id).reaction

    return model_preview_json,default_setting_json

def get_preview_default_json(norm_model,met_C):
    model_preview_json={}
    default_setting_json={}

    default_setting_json['Oxygen']={}
    default_setting_json['Oxygen']['model_id']='None'
    default_setting_json['Oxygen']['model_correspondind']='None'
    default_setting_json['Oxygen']['value']=0 
    default_setting_json['Carbon_source']={}
    default_setting_json['Carbon_source']['substrate_id']='None'
    default_setting_json['Carbon_source']['substrate_model_correspondind']='None'
    default_setting_json['Carbon_source']['substrate_name']='None'
    default_setting_json['Carbon_source']['uptake_rate']=0 
    default_setting_json['Nitrogen_source']={}
    default_setting_json['Nitrogen_source']['substrate_id']='None'
    default_setting_json['Nitrogen_source']['substrate_model_correspondind']='None'
    default_setting_json['Nitrogen_source']['substrate_name']='None' 
    default_setting_json['Nitrogen_source']['uptake_rate']=0 
    default_setting_json['ATPM']={}
    default_setting_json['ATPM']['model_id']='None'
    default_setting_json['ATPM']['model_correspondind']='None'
    default_setting_json['ATPM']['value']=0 

    O2_list=['EX_o2_e_reverse','EX_o2_e','EX_cpd00007_b','EX_cpd00007_e0','EX_o2(e)','EX_MNXM735438@BOUNDARY','mnxr102090c2b','Ex_o2','Ex_o2_reverse','R1287']
    #ATPM_list=['ATPM','rxn00062_c0','DM_atp_c_','mnxr153054c3']
    ATPM_list=['ATPM','rxn00062','rxn00062_c0','DM_atp_c_','mnxr153054c3','R00086','META:ATPASE-RXN','ATPASE-RXN','MNXR96131','rxn09694', 'rxn11300','rxn09694_c0', 'rxn11300_c0','ATPM_c', 'ATPM_LPAREN_NGAM_RPAREN_', 'DM_atp_c', 'DM_atp[c]']

    #输出培养基
    medium_name_list=[]
    carbon_source_list=[]
    nitrogen_source_list=[]
    try:
        find_boundary_types(norm_model,'exchange')
    except:
        print('No medium information!')
    else:       
        for eachm in medium(norm_model):
            is_C_met='No'
            is_N_met='No'
            is_O_met='No'
            is_P_met='No'
            is_S_met='No'            
            medium_name=norm_model.reactions.get_by_id(eachm).name.replace(' exchange','').replace('exchange','').replace('EXCHANGE','')
            if eachm not in model_preview_json.keys():
                model_preview_json[eachm]={}            
            if len(medium_name)==0:
                medium_name=eachm
            model_preview_json[eachm]['Metabolite']=medium_name
            model_preview_json[eachm]['Lower Bound']=norm_model.reactions.get_by_id(eachm).lower_bound
            model_preview_json[eachm]['Upper Bound']=norm_model.reactions.get_by_id(eachm).upper_bound
            model_preview_json[eachm]['Medium']='Yes'
            #medium_name_list.append(medium_name)
            try:
                reactant_id=norm_model.reactions.get_by_id(eachm).reactants[0].id
            except:
                reactant_id=norm_model.reactions.get_by_id(eachm).products[0].id   
            #C
            try:
                product_C_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['C']
            except:
                try:
                    product_C_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['C']
                except:
                    print('No formula for %s!'%reactant_id)
                else: 
                    is_C_met='Yes'
            else:
                is_C_met='Yes'

            if is_C_met=='Yes':
                if abs(norm_model.reactions.get_by_id(eachm).lower_bound) <1000 and abs(norm_model.reactions.get_by_id(eachm).lower_bound) >0:
                    if 'Carbon_source' in default_setting_json.keys():
                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).lower_bound):
                            default_setting_json['Carbon_source']['substrate_id']=reactant_id
                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                            default_setting_json['Carbon_source']['substrate_name']=medium_name
                            default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)             
                    else:
                        default_setting_json['Carbon_source']={}
                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                        default_setting_json['Carbon_source']['substrate_id']=reactant_id
                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                        default_setting_json['Carbon_source']['substrate_name']=medium_name
                        default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                elif norm_model.reactions.get_by_id(eachm).upper_bound>0 and norm_model.reactions.get_by_id(eachm).upper_bound<1000:
                    if 'Carbon_source' in default_setting_json.keys():
                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).upper_bound):
                            default_setting_json['Carbon_source']['substrate_id']=reactant_id
                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                            default_setting_json['Carbon_source']['substrate_name']=medium_name
                            default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).upper_bound)             
                    else:
                        default_setting_json['Carbon_source']={}
                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                        default_setting_json['Carbon_source']['substrate_id']=reactant_id
                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachm                    
                        default_setting_json['Carbon_source']['substrate_name']=medium_name
                        default_setting_json['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).upper_bound)                
                #carbon_source_list.append(medium_name)
                if eachm not in model_preview_json.keys():
                    model_preview_json[eachm]={}
                if 'Lower Bound' not in model_preview_json[eachm].keys():
                    model_preview_json[eachm]['Lower Bound']=norm_model.reactions.get_by_id(eachm).lower_bound
                    model_preview_json[eachm]['Upper Bound']=norm_model.reactions.get_by_id(eachm).upper_bound
                model_preview_json[eachm]['Carbon']='Yes'
            #N
            try:
                product_N_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['N']
            except:
                try:
                    product_N_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['N']
                except:
                    print('No formula for %s!'%reactant_id)
                    if reactant_id=='nh4_e':
                        default_setting_json['Nitrogen_source']={}
                        #default_setting_json['Nitrogen_source']['substrate_model_correspondind']='EX_'+reactant_id
                        default_setting_json['Nitrogen_source']['substrate_id']=reactant_id
                        default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachm
                        default_setting_json['Nitrogen_source']['substrate_name']=medium_name
                        default_setting_json['Nitrogen_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                        nitrogen_source_list.append(medium_name)  
                else: 
                    is_N_met='Yes'
            else:
                is_N_met='Yes'

            if is_N_met=='Yes':
                if reactant_id=='nh3_e':
                    default_setting_json['Nitrogen_source']={}
                    #default_setting_json['Nitrogen_source']['substrate_model_correspondind']='EX_'+reactant_id
                    default_setting_json['Nitrogen_source']['substrate_id']=reactant_id
                    default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachm
                    default_setting_json['Nitrogen_source']['substrate_name']=medium_name
                    default_setting_json['Nitrogen_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                    nitrogen_source_list.append(medium_name)  
                else:
                    nitrogen_source_list.append(medium_name) 
                    
                if eachm not in model_preview_json.keys():
                    model_preview_json[eachm]={}
                if 'Lower Bound' not in model_preview_json[eachm].keys():
                    model_preview_json[eachm]['Lower Bound']=norm_model.reactions.get_by_id(eachm).lower_bound
                    model_preview_json[eachm]['Upper Bound']=norm_model.reactions.get_by_id(eachm).upper_bound
                model_preview_json[eachm]['Nitrogen']='Yes'
            #O
            try:
                product_O_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['O']
            except:
                try:
                    product_O_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['O']
                except:
                    print('No formula for %s!'%reactant_id)
                else: 
                    is_O_met='Yes'
            else:
                is_O_met='Yes'

            if is_O_met=='Yes':              
                #carbon_source_list.append(medium_name)
                if eachm not in model_preview_json.keys():
                    model_preview_json[eachm]={}
                if 'Lower Bound' not in model_preview_json[eachm].keys():
                    model_preview_json[eachm]['Lower Bound']=norm_model.reactions.get_by_id(eachm).lower_bound
                    model_preview_json[eachm]['Upper Bound']=norm_model.reactions.get_by_id(eachm).upper_bound
                model_preview_json[eachm]['Oxygen']='Yes'
                
            #P
            try:
                product_P_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['P']
            except:
                try:
                    product_P_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['P']
                except:
                    print('No formula for %s!'%reactant_id)
                else: 
                    is_P_met='Yes'
            else:
                is_P_met='Yes'

            if is_P_met=='Yes':              
                #carbon_source_list.append(medium_name)
                if eachm not in model_preview_json.keys():
                    model_preview_json[eachm]={}
                if 'Lower Bound' not in model_preview_json[eachm].keys():
                    model_preview_json[eachm]['Lower Bound']=norm_model.reactions.get_by_id(eachm).lower_bound
                    model_preview_json[eachm]['Upper Bound']=norm_model.reactions.get_by_id(eachm).upper_bound
                model_preview_json[eachm]['Phosphorus']='Yes'
                
            #S
            try:
                product_S_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['S']
            except:
                try:
                    product_S_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['S']
                except:
                    print('No formula for %s!'%reactant_id)
                else: 
                    is_S_met='Yes'
            else:
                is_S_met='Yes'

            if is_S_met=='Yes':              
                #carbon_source_list.append(medium_name)
                if eachm not in model_preview_json.keys():
                    model_preview_json[eachm]={}
                if 'Lower Bound' not in model_preview_json[eachm].keys():
                    model_preview_json[eachm]['Lower Bound']=norm_model.reactions.get_by_id(eachm).lower_bound
                    model_preview_json[eachm]['Upper Bound']=norm_model.reactions.get_by_id(eachm).upper_bound
                model_preview_json[eachm]['Sulfur']='Yes'
                
            if 'Medium' not in model_preview_json[eachm].keys():
                model_preview_json[eachm]['Medium']='No'
            if 'Carbon' not in model_preview_json[eachm].keys():
                model_preview_json[eachm]['Carbon']='No'
            if 'Nitrogen' not in model_preview_json[eachm].keys():
                model_preview_json[eachm]['Nitrogen']='No'
            if 'Oxygen' not in model_preview_json[eachm].keys():
                model_preview_json[eachm]['Oxygen']='No'
            if 'Phosphorus' not in model_preview_json[eachm].keys():
                model_preview_json[eachm]['Phosphorus']='No'
            if 'Sulfur' not in model_preview_json[eachm].keys():
                model_preview_json[eachm]['Sulfur']='No'   
                
    #ATPM
    for eachatpm in ATPM_list:
        try:                
            ATPM=norm_model.reactions.get_by_id(eachatpm).lower_bound
            print(eachatpm,ATPM)
        except:
            print('Model does not contain ATP Maintenance Requirement')
            #default_setting_json['ATPM']='none'
        else:
            default_setting_json['ATPM']={}
            default_setting_json['ATPM']['model_id']=eachatpm
            default_setting_json['ATPM']['model_correspondind']=eachatpm
            default_setting_json['ATPM']['value']=ATPM 
            break

    #有氧厌氧
    default_setting_json['Oxygen']={}
    default_setting_json['Oxygen']['model_id']='None'
    default_setting_json['Oxygen']['model_correspondind']='None'
    default_setting_json['Oxygen']['value']=0 
    for eachatpm in O2_list:
        try:
            O2=norm_model.reactions.get_by_id(eachatpm).lower_bound
        except:
            print('No O2 reaction!')         
        else:
            default_setting_json['Oxygen']={}
            default_setting_json['Oxygen']['model_id']=eachatpm
            default_setting_json['Oxygen']['model_correspondind']=eachatpm
            default_setting_json['Oxygen']['value']=abs(O2)
            break

#     model_preview_json['Medium']=medium_name_list
#     model_preview_json['Carbon_source']=carbon_source_list
#     model_preview_json['Nitrogen_source']=nitrogen_source_list    
    #输出目标函数
    if len(str(norm_model.objective.expression))>300:
        model_preview_json['Optimization Objectives']='Biomass' 
        regex = re.compile('BIOMASS', flags=re.IGNORECASE)
        if len(norm_model.reactions.query(regex))>1:
            for eglc in norm_model.reactions.query(regex):
                with norm_model as model:
                    model.objective=eglc.id
                    try:
                        opt_r=model.optimize()
                    except:
                        pass
                    else:
                        if opt_r[eglc.id]>0:
                            default_setting_json['Object']={}
                            default_setting_json['Object']['model_id']=eglc.id.replace('EX_','')
                            default_setting_json['Object']['model_correspondind']=eglc.id
                            default_setting_json['Object']['substrate_name']='Biomass_core'
                            #default_setting_json['Object']['model_reaction']=eglc.reaction
                            if eglc.name:
                                model_preview_json['Optimization Objectives']=eglc.name
                            else:
                                model_preview_json['Optimization Objectives']=eglc.id
        else: 
            for eglc in norm_model.reactions.query(regex):
                default_setting_json['Object']={}
                default_setting_json['Object']['model_id']=eglc.id.replace('EX_','')
                default_setting_json['Object']['model_correspondind']=eglc.id
                default_setting_json['Object']['substrate_name']='Biomass'
                #default_setting_json['Object']['model_reaction']=eglc.reaction
                if eglc.name:
                    model_preview_json['Optimization Objectives']=eglc.name
                else:
                    model_preview_json['Optimization Objectives']=eglc.id

    elif str(norm_model.objective.expression)=='0':
        model_preview_json['Optimization Objectives']='Not set' 
        default_setting_json['Object']={}
        default_setting_json['Object']['model_id']='None'
        default_setting_json['Object']['model_correspondind']='None'
        #default_setting_json['Object']['model_reaction']='None'
        default_setting_json['Object']['substrate_name']='None'
    else:
        obj_id=str(norm_model.objective.expression).split(' -')[0].split('*')[1]
        obj_name=norm_model.reactions.get_by_id(obj_id).name.replace(' exchange','')
        if len(obj_name)<1:
            obj_name=obj_name+'('+obj_id+')'
        model_preview_json['Optimization Objectives']=obj_name 
        # if obj_name:
        #     model_preview_json['Optimization Objectives']=obj_name
        # else:
        #     model_preview_json['Optimization Objectives']=obj_id
        if re.search('BIOMASS',obj_id, re.I):
            default_setting_json['Object']={}
            default_setting_json['Object']['model_id']=obj_id.replace('EX_','')
            default_setting_json['Object']['model_correspondind']=obj_id 
            default_setting_json['Object']['substrate_name']='Biomass'
            #default_setting_json['Object']['model_reaction']=norm_model.reactions.get_by_id(obj_id).reaction
        else:
            default_setting_json['Object']={}
            default_setting_json['Object']['model_id']=obj_id.replace('EX_','') 
            default_setting_json['Object']['model_correspondind']=obj_id
            default_setting_json['Object']['substrate_name']=obj_name
            #default_setting_json['Object']['model_reaction']=norm_model.reactions.get_by_id(obj_id).reaction
    if len(model_preview_json['Optimization Objectives'])==0:
        model_preview_json['Optimization Objectives']='not parsed successfully'
    model_preview_json['Preview_static']='This model has %s reactions, %s metabolites and %s genes. The optimization objective is %s. The input metabolites (including carbon sources, nitrogen sources, etc.) of the model are listed below and can be updated.'%(len(norm_model.reactions),len(norm_model.metabolites),len(norm_model.genes),model_preview_json['Optimization Objectives'])
        
    return model_preview_json,default_setting_json

def get_met_C(met_C_file):
    met_C = pd.read_csv(met_C_file,sep='\t', index_col=0)
    return met_C

def create_model_change_json(model):
    # 生成model_change
    model_change_json={}
    model_change_json['reaction_infor']={}
    model_change_json['gene_infor']={}
    model_change_json['exchange_infor']={}  #add_reaction后 model.exchanges没了
    model_change_json['bigg_metabolite_ids']={}  
    model_change_json['add_reaction_infor']={}    # 与modify生成的model_change_json_final.keys()保持一致
    model_change_json['del_reaction_infor']={}
    model_change_json['del_gene_infor']={}
    for eachr in model.reactions:
        model_change_json['reaction_infor'][eachr.id]={}
        model_change_json['reaction_infor'][eachr.id]['reaction_name']=eachr.name
        model_change_json['reaction_infor'][eachr.id]['equation']=eachr.reaction
        model_change_json['reaction_infor'][eachr.id]['equation by name']=eachr.build_reaction_string(True).replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
        model_change_json['reaction_infor'][eachr.id]['gene_reaction_rule']=eachr.gene_reaction_rule
        model_change_json['reaction_infor'][eachr.id]['lower_bound']=model.reactions.get_by_id(eachr.id).lower_bound
        model_change_json['reaction_infor'][eachr.id]['upper_bound']=model.reactions.get_by_id(eachr.id).upper_bound
    for eachg in model.genes:
        if eachg.name:
            model_change_json['gene_infor'][eachg.id]={}
            model_change_json['gene_infor'][eachg.id]['gene_name']=eachg.name
            reaction_id_list=[]
            reaction_name_list=[]
            reaction_equ_list=[]
            for eachr in model.genes.get_by_id(eachg.id).reactions:
                reaction_name_list.append(eachr.name)
                reaction_equ_list.append(eachr.reaction)
            model_change_json['gene_infor'][eachg.id]['reaction_name_list']='; '.join(reaction_name_list) 
            model_change_json['gene_infor'][eachg.id]['reaction_equ_list']='; '.join(reaction_equ_list)    
    if len(exchanges(model))>0:
        for eachr in exchanges(model):
            model_change_json['exchange_infor'][eachr]={}
            model_change_json['exchange_infor'][eachr]['reaction_name']=model.reactions.get_by_id(eachr).name
            model_change_json['exchange_infor'][eachr]['lower_bound']=model.reactions.get_by_id(eachr).lower_bound
            model_change_json['exchange_infor'][eachr]['upper_bound']=model.reactions.get_by_id(eachr).upper_bound
    else:
        model_change_json['exchange_infor']={}
    for eachm in model.metabolites:
        model_change_json['bigg_metabolite_ids'][eachm.id]={}
        model_change_json['bigg_metabolite_ids'][eachm.id]['bigg_id']=eachm.id
        model_change_json['bigg_metabolite_ids'][eachm.id]['name']=eachm.name
    return model_change_json



def exchanges(model):
    """Exchange reactions in model.
    Reactions that exchange mass with the exterior. Uses annotations
    and heuristics to exclude non-exchanges such as sink reactions.
    """
    exchange_list=[]
    try: 
        find_boundary_types(model,'exchange')
    except:
        pass
    else:
        for reachr in find_boundary_types(model,'exchange'):
            exchange_list.append(reachr.id)
    return exchange_list

def medium(model):
    def is_active(reaction):
        """Determine if a boundary reaction permits flux towards creating
        metabolites
        """

        return (bool(reaction.products) and (reaction.upper_bound > 0)) or (
            bool(reaction.reactants) and (reaction.lower_bound < 0)
        )

    def get_active_bound(reaction):
        """For an active boundary reaction, return the relevant bound"""
        if reaction.reactants:
            return -reaction.lower_bound
        elif reaction.products:
            return reaction.upper_bound

    return {
        rxn.id: get_active_bound(rxn) for rxn in find_boundary_types(model,'exchange') if is_active(rxn)
    }

def exchanges_substrateold(model,substrate_use_dict):
    """Exchange reactions in model.
    Reactions that exchange mass with the exterior. Uses annotations
    and heuristics to exclude non-exchanges such as sink reactions.
    """
    substrate_name_list=[]
    for reachr in model.reactions:
        r_id_upper=reachr.id.upper()
        r_name_upper=reachr.name.upper()
        if re.search('EXCHANGE',r_name_upper, re.I):
            try:
                r_id_upper.split('EX_')[1]
            except:
                pass
            else:
                substrate_use_dict[reachr.id]={}
                substrate_use_dict[reachr.id]['substrate_id']=r_id_upper.split('EX_')[1]
                substrate_use_dict[reachr.id]['substrate_model_correspondind']=reachr.id
                substrate_use_dict[reachr.id]['substrate_name']=r_id_upper.replace('EXCHANGE','')+'('+r_id_upper.split('EX_')[1]+')'
                substrate_use_dict[reachr.id]['type']="model exchange ori"     
                substrate_name_list.append(reachr.name)
        elif r_id_upper.startswith('EX_'):
            substrate_use_dict[reachr.id]={}
            substrate_use_dict[reachr.id]['substrate_id']=r_id_upper
            substrate_use_dict[reachr.id]['substrate_model_correspondind']=reachr.id
            substrate_use_dict[reachr.id]['substrate_name']=r_id_upper.replace('EXCHANGE','')+'('+r_id_upper.split('EX_')[1]+')'
            substrate_use_dict[reachr.id]['type']="model exchange ori"     
            substrate_name_list.append(reachr.name)
            
    return [substrate_use_dict,substrate_name_list]

def exchanges_substrate(norm_model,substrate_use_dict):
    """Exchange reactions in model.
    Reactions that exchange mass with the exterior. Uses annotations
    and heuristics to exclude non-exchanges such as sink reactions.
    """
    substrate_name_list=[]
    #exchange_r_list=exchanges(model)#反应ID
#    medium_r_list=medium(model)#反应ID
#     for ex_r in exchange_r_list:
#         exchange_r=norm_model.reactions.get_by_id(ex_r)
#         r_id_upper=exchange_r.id.upper()
#         r_name_upper=exchange_r.name.upper()
#         if re.search('EXCHANGE',r_name_upper, re.I):
#             try:
#                 r_id_upper.split('EX_')[1]
#             except:
#                 pass
#             else:
#                 substrate_use_dict[ex_r]={}
#                 substrate_use_dict[ex_r]['substrate_id']=r_id_upper.split('EX_')[1]
#                 substrate_use_dict[ex_r]['substrate_model_correspondind']=ex_r
#                 substrate_use_dict[ex_r]['substrate_name']=r_id_upper.replace('EXCHANGE','')+'('+r_id_upper.split('EX_')[1]+')'
#                 substrate_use_dict[ex_r]['type']="model exchange ori"     
#                 substrate_name_list.append(exchange_r.name)
#         elif r_id_upper.startswith('EX_'):
#             substrate_use_dict[ex_r]={}
#             substrate_use_dict[ex_r]['substrate_id']=r_id_upper
#             substrate_use_dict[ex_r]['substrate_model_correspondind']=ex_r
#             substrate_use_dict[ex_r]['substrate_name']=r_id_upper.replace('EXCHANGE','')+'('+r_id_upper.split('EX_')[1]+')'
#             substrate_use_dict[ex_r]['type']="model exchange ori"     
#             substrate_name_list.append(exchange_r.name)
#     for me_r in medium_r_list:
#         medium_r=norm_model.reactions.get_by_id(me_r)
#         r_id_upper=medium_r.id.upper()
#         r_name_upper=medium_r.name.upper()
#         if re.search('EXCHANGE',r_name_upper, re.I):
#             try:
#                 r_id_upper.split('EX_')[1]
#             except:
#                 pass
#             else:
#                 if me_r not in substrate_use_dict.keys():
#                     substrate_use_dict[me_r]={}
#                     substrate_use_dict[me_r]['substrate_id']=r_id_upper.split('EX_')[1]
#                     substrate_use_dict[me_r]['substrate_model_correspondind']=me_r
#                     substrate_use_dict[me_r]['substrate_name']=r_id_upper.replace('EXCHANGE','')+'('+r_id_upper.split('EX_')[1]+')'
#                     substrate_use_dict[me_r]['type']="model exchange ori"     
#                     substrate_name_list.append(medium_r.name)
#         elif r_id_upper.startswith('EX_'):
#             if me_r not in substrate_use_dict.keys():
#                 substrate_use_dict[me_r]={}
#                 substrate_use_dict[me_r]['substrate_id']=r_id_upper
#                 substrate_use_dict[me_r]['substrate_model_correspondind']=me_r
#                 substrate_use_dict[me_r]['substrate_name']=r_id_upper.replace('EXCHANGE','')+'('+r_id_upper.split('EX_')[1]+')'
#                 substrate_use_dict[me_r]['type']="model exchange ori"     
#                 substrate_name_list.append(medium_r.name)     
    for eachm in medium(norm_model):
        is_C_met='No'
        is_N_met='No'
        medium_name=norm_model.reactions.get_by_id(eachm).name.replace(' exchange','')
        try:
            reactant_id=norm_model.reactions.get_by_id(eachm).reactants[0].id
        except:
            reactant_id=norm_model.reactions.get_by_id(eachm).products[0].id   
        try:
            product_C_num=get_composition_from_formula(norm_model.metabolites.get_by_id(reactant_id).formula)['C']
        except:
            try:
                product_C_num=get_composition_from_formula(met_C.loc[reactant_id,'formula'])['C']
            except:
                print('No formula for %s!'%reactant_id)
            else: 
                is_C_met='Yes'
        else:
            is_C_met='Yes'

        if is_C_met=='Yes':
#             if abs(norm_model.reactions.get_by_id(eachm).lower_bound) <1000:
#                 if 'Carbon_source' in substrate_use_dict.keys():
#                     if substrate_use_dict['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).lower_bound):
#                         substrate_use_dict['Carbon_source']={}
#                         substrate_use_dict['Carbon_source']['substrate_id']=reactant_id
#                         substrate_use_dict['Carbon_source']['substrate_model_correspondind']=eachm
#                         substrate_use_dict['Carbon_source']['substrate_name']=medium_name
#                         substrate_use_dict['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
#                         substrate_use_dict['Carbon_source']['type']="default"     
#                         substrate_name_list.append(medium_name)            
#                 else:
#                     substrate_use_dict['Carbon_source']={}
#                     substrate_use_dict['Carbon_source']['substrate_id']=reactant_id
#                     substrate_use_dict['Carbon_source']['substrate_model_correspondind']=eachm
#                     substrate_use_dict['Carbon_source']['substrate_name']=medium_name
#                     substrate_use_dict['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
#                     substrate_use_dict['Carbon_source']['type']="default"     
#                     substrate_name_list.append(medium_name) 
            if abs(norm_model.reactions.get_by_id(eachm).lower_bound) <1000 and abs(norm_model.reactions.get_by_id(eachm).lower_bound) >0:
                if 'Carbon_source' in substrate_use_dict.keys():
                    if substrate_use_dict['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).lower_bound):
                        substrate_use_dict['Carbon_source']['substrate_id']=reactant_id
                        substrate_use_dict['Carbon_source']['substrate_model_correspondind']=eachm                    
                        substrate_use_dict['Carbon_source']['substrate_name']=medium_name
                        substrate_use_dict['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)    
                        substrate_use_dict['Carbon_source']['type']="default"     
                        substrate_name_list.append(medium_name) 
                else:
                    substrate_use_dict['Carbon_source']={}
                    #substrate_use_dict['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                    substrate_use_dict['Carbon_source']['substrate_id']=reactant_id
                    substrate_use_dict['Carbon_source']['substrate_model_correspondind']=eachm                    
                    substrate_use_dict['Carbon_source']['substrate_name']=medium_name
                    substrate_use_dict['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).lower_bound)
                    substrate_use_dict['Carbon_source']['type']="default"     
                    substrate_name_list.append(medium_name) 
            elif norm_model.reactions.get_by_id(eachm).upper_bound>0 and norm_model.reactions.get_by_id(eachm).upper_bound<1000:
                if 'Carbon_source' in substrate_use_dict.keys():
                    if substrate_use_dict['Carbon_source']['uptake_rate'] < abs(norm_model.reactions.get_by_id(eachm).upper_bound):
                        substrate_use_dict['Carbon_source']['substrate_id']=reactant_id
                        substrate_use_dict['Carbon_source']['substrate_model_correspondind']=eachm                    
                        substrate_use_dict['Carbon_source']['substrate_name']=medium_name
                        substrate_use_dict['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).upper_bound)     
                        substrate_use_dict['Carbon_source']['type']="default"     
                        substrate_name_list.append(medium_name) 
                else:
                    substrate_use_dict['Carbon_source']={}
                    #substrate_use_dict['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                    substrate_use_dict['Carbon_source']['substrate_id']=reactant_id
                    substrate_use_dict['Carbon_source']['substrate_model_correspondind']=eachm                    
                    substrate_use_dict['Carbon_source']['substrate_name']=medium_name
                    substrate_use_dict['Carbon_source']['uptake_rate']=abs(norm_model.reactions.get_by_id(eachm).upper_bound)
                    substrate_use_dict['Carbon_source']['type']="default"     
                    substrate_name_list.append(medium_name) 
    return [substrate_use_dict,substrate_name_list]

#新定义的函数
def close_model_default_carbon(model,met_C):
    #关闭碳源
    default_C_list=[]
    exclude_C_reaction_list=['EX_BIOTIN_reverse','EX_BIOTIN','EX_biotin_e', 'EX_btn', 'EX_btn_e', 'EX_btn_e_', 'EX_btn(e)', 'EX_btn[e]', 'EX_btn_LPAREN_e_RPAREN_','EX_biotin_e_reverse', 'EX_btn_reverse', 'EX_btn_e_reverse', 'EX_btn_e__reverse', 'EX_btn(e)_reverse', 'EX_btn[e]_reverse', 'EX_btn_LPAREN_e_RPAREN__reverse','META:TRANS-RXN0-240','MNXR96333']
    for eachr in model.reactions.keys():
        if "EX_" in eachr:
            #print(eachr)
            #r_bounds=norm_model.reactions.get_by_id(eachr).bounds
            product_equ=model.get_reaction_by_id(eachr)
            # 初始化计数器
            count_negative = 0
            # 迭代字典的值
            for value in product_equ.metabolites.values():
                if value < 0:
                    count_negative += 1
            if count_negative>1:
                pass
            else:
                for eachm in product_equ.metabolites.keys():
                    if product_equ.metabolites[eachm]<0:
                        #print(eachr,eachm.id,eachm.formula)
                        #met_coef=abs(product_equ.get_coefficient(eachm.id))
                        try:
                            met_formula=get_composition_from_formula(model.get_metabolite_by_id(eachm).formula)
                        except:
                            try:
                                met_formula=get_composition_from_formula(met_C.loc[eachm,'formula'])
                                #print(eachr,eachm.id,met_formula)s
                            except:
                                print('No formula!')
                            else:
                                if 'C' in met_formula.keys():#是碳源才关闭
        #                             if eachr not in exclude_C_reaction_list:
        #                                 model.reactions.get_by_id(eachr).bounds = (0, 1000) 
        #                             default_C_list.append(eachr)
                                    if model.get_reaction_by_id(eachr).lower_bound<0:
                                        if eachr not in exclude_C_reaction_list:
                                            model.get_reaction_by_id(eachr).lower_bound = 0
                                            model.get_reaction_by_id(eachr).upper_bound = 1000
                                        default_C_list.append(eachr)
                                    else:
                                        if eachr not in exclude_C_reaction_list:
                                            model.get_reaction_by_id(eachr).lower_bound = 0
                                            model.get_reaction_by_id(eachr).upper_bound = 0
                                        default_C_list.append(eachr)  
                        else:
                            if 'C' in met_formula.keys():#是碳源才关闭
        #                         if eachr not in exclude_C_reaction_list:
        #                             model.reactions.get_by_id(eachr).bounds = (0, 1000) 
        #                         default_C_list.append(eachr)
                                if model.get_reaction_by_id(eachr).lower_bound<0:
                                    if eachr not in exclude_C_reaction_list:
                                        model.get_reaction_by_id(eachr).lower_bound = 0
                                        model.get_reaction_by_id(eachr).upper_bound = 1000 
                                    default_C_list.append(eachr)
                                else:
                                    if eachr not in exclude_C_reaction_list:
                                        model.get_reaction_by_id(eachr).lower_bound = 0
                                        model.get_reaction_by_id(eachr).upper_bound = 0
                                    default_C_list.append(eachr)
            count_positive = 0
            # 迭代字典的值
            for value in product_equ.metabolites.values():
                if value > 0:
                    count_positive += 1
            if count_positive>1:
                pass
            else:
                for eachm in product_equ.metabolites.keys():
                    if product_equ.metabolites[eachm]>0:
                    #print(eachr,eachm.id,eachm.formula)
                    #met_coef=abs(product_equ.get_coefficient(eachm.id))
                        try:
                            met_formula=get_composition_from_formula(model.get_metabolite_by_id(eachm).formula)
                        except:
                            try:
                                met_formula=get_composition_from_formula(met_C.loc[eachm,'formula'])
                                #print(eachr,eachm.id,met_formula)
                            except:
                                print('No formula!')
                            else:
                                if 'C' in met_formula.keys():#是碳源才关闭
                                    if model.get_reaction_by_id(eachr).lower_bound<0:
                                        if eachr not in exclude_C_reaction_list:
                                            model.get_reaction_by_id(eachr).lower_bound = 0
                                            model.get_reaction_by_id(eachr).upper_bound = 1000
                                        default_C_list.append(eachr)
                                    else:
                                        if eachr not in exclude_C_reaction_list:
                                            model.get_reaction_by_id(eachr).lower_bound = 0
                                            model.get_reaction_by_id(eachr).upper_bound = 0 
                                        default_C_list.append(eachr)                                
                        else:
                            if 'C' in met_formula.keys():#是碳源才关闭
                                if model.get_reaction_by_id(eachr).lower_bound<0:
                                    if eachr not in exclude_C_reaction_list:
                                        model.get_reaction_by_id(eachr).lower_bound = 0
                                        model.get_reaction_by_id(eachr).upper_bound = 1000
                                    default_C_list.append(eachr)
                                else:
                                    if eachr not in exclude_C_reaction_list:
                                        model.get_reaction_by_id(eachr).lower_bound = 0
                                        model.get_reaction_by_id(eachr).upper_bound = 0 
                                    default_C_list.append(eachr)                                                            
                            
    return [model,default_C_list]

#新定义的函数
def get_QC_detail_df(model):
    pfba_solution = cobra.flux_analysis.pfba(model)
    model_pfba_solution = pfba_solution.to_frame()
    model_pfba_solution_detail = {}
    for index, row in model_pfba_solution.iterrows():
        reaction_detail = model.reactions.get_by_id(index)
        model_pfba_solution_detail[index]={}
        if abs(row['fluxes'])>1e-9:
            useflux=round(row['fluxes'],6)
        else:
            useflux=0
        model_pfba_solution_detail[index]['reaction_id'] = index
        model_pfba_solution_detail[index]['reaction_name'] = reaction_detail.name
        model_pfba_solution_detail[index]['fluxes'] = useflux    
        if reaction_detail.gene_reaction_rule:
            model_pfba_solution_detail[index]['gene id'] =reaction_detail.gene_reaction_rule
        else:
            model_pfba_solution_detail[index]['gene id'] = 'none'
        model_pfba_solution_detail[index]['equ'] = reaction_detail.reaction
        model_pfba_solution_detail[index]['equ_name'] = reaction_detail.build_reaction_string(True).replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
    model_pfba_solution_detail_df = pd.DataFrame(model_pfba_solution_detail).T
    model_pfba_solution_detail_df=model_pfba_solution_detail_df[abs(model_pfba_solution_detail_df['fluxes'])>0]
    return model_pfba_solution_detail_df

def model_glucose_information(model,model_check_json,default_C_list):
    #设置葡萄糖浓度为10时APT通量和途径结果
    GLU_R_list=['EX_glc__D_e','META:TRANS-RXN0-574','TRANS-RXN0-574','rxn08617','rxn09679','rxn09875','MNXR100188','EX_glc', 'EX_glc_D_e', 'EX_glc_D[e]', 'EX_glc_e', 'EX_glc_e_', 'EX_glc(e)', 'EX_glc_LPAREN_e_RPAREN_','R1070']
    GLU_R_use_list=[]
    for each_GLU in GLU_R_list:
        try:
            model.reactions.get_by_id(each_GLU)
        except:
            pass
        else:
            GLU_R_use_list.append(each_GLU)    
    if len(GLU_R_use_list)==0:
        GLU_R_use=default_C_list[0]
        try:
            C_id=model.reactions.get_by_id(GLU_R_use).reactants[0].name
            C_name=model.reactions.get_by_id(GLU_R_use).reactants[0].name
        except:
            C_id=model.reactions.get_by_id(GLU_R_use).products[0].id
            C_name=model.reactions.get_by_id(GLU_R_use).products[0].name
        if C_name:
            C_name=C_id
        model_check_json['Carbon Source Information']='There is no glucose. Use model default setting (%s)!'%C_name
    elif len(GLU_R_use_list) >1:
        GLU_R_use=GLU_R_use_list[0]
        try:
            C_id=model.reactions.get_by_id(GLU_R_use).reactants[0].name
            C_name=model.reactions.get_by_id(GLU_R_use).reactants[0].name
        except:
            C_id=model.reactions.get_by_id(GLU_R_use).products[0].id
            C_name=model.reactions.get_by_id(GLU_R_use).products[0].name
        if C_name:
            C_name=C_id
        model_check_json['Carbon Source Information']='There are multi glucose reaction('+(' and ').join(GLU_R_use_list)+'). Select %s(%s) as carbon source!'%(GLU_R_use,C_name)
    else:
        GLU_R_use=GLU_R_use_list[0]
        try:
            C_id=model.reactions.get_by_id(GLU_R_use).reactants[0].name
            C_name=model.reactions.get_by_id(GLU_R_use).reactants[0].name
        except:
            C_id=model.reactions.get_by_id(GLU_R_use).products[0].id
            C_name=model.reactions.get_by_id(GLU_R_use).products[0].name
        if C_name:
            C_name=C_id
        model_check_json['Carbon Source Information']='Use %s(%s) as carbon source!'%(GLU_R_use,C_name)
    return [GLU_R_use,model_check_json]

#新函数
def change_reaction_bound(model):
    for eachr in model.reactions:
        if eachr.lower_bound < 0 and eachr.upper_bound > 0:
                eachr.bounds=(-1000,1000)
        if eachr.lower_bound > 0 and eachr.upper_bound > 0:
                eachr.bounds=(0,1000)
        if eachr.lower_bound < 0 and eachr.upper_bound < 0:
                eachr.bounds=(-1000,0)    
    return model

def Draw_MDF_By_Growth_rate(MDF_list,xmin,xlab,xmax,step,ytick,png_file):
    plt.figure(figsize=(15, 10), dpi=300)
    MDF_list=MDF_list.sort_index(ascending=True) 
    plt.plot(MDF_list.index, MDF_list[MDF_list.columns[0]], color="cyan", linewidth=3.0, linestyle="-", label=MDF_list.columns[0])
    plt.plot(MDF_list.index, MDF_list[MDF_list.columns[1]], color="darkorange", linewidth=3.0, linestyle="-", label=MDF_list.columns[1])
    font1 = {#'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 23,
    }

    font2 = {#'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 30,
    }

    plt.ylabel("MDF of pathways (kJ/mol)",font2)
    #xlab = "Growth rate ($\mathregular{h^-1}$)"
    plt.xlabel(xlab,font2)

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_position(('data', xmin))
    ax.spines['bottom'].set_position(('data', 0))
    plt.legend(loc="lower left",prop=font1)
    plt.xlim(xmin, xmax)
    plt.ylim(np.min(ytick), np.max(ytick))

    plt.tick_params(labelsize=23)

    plt.xticks(generate_list(xmin, xmax, step))
    plt.yticks(ytick)

    plt.savefig(png_file)
    plt.show()

def generate_list(min_value, max_value, step):
    result = []
    current_value = min_value
    while current_value <= max_value:
        result.append(current_value)
        current_value += step
    return result

def find_gene_knockout_reactions(norm_model,gene_id):
    reaction_list = []
    for reaction in norm_model.reactions.values():
        gene_list = reaction.gene_reaction_rule.split(" ")
        if gene_id in gene_list:
            reaction_list.append(reaction.id)
    return reaction_list

def build_model_from_local_file(model_pfba_solution,rxn_switch,coordinates,col_equ,col_flux):#新加函数
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    metlink_model = Model('metlink_model')

    for index, row in model_pfba_solution.iterrows():
        new_equ=row[col_equ].replace('__D_e','_e').replace('__L_e','_e')
        uni_str_line=get_reaction_uni_met(new_equ)    
        if len(uni_str_line)==1:  
            pass       
        else:
            metlink_model=get_model_met_from_reaction_string(metlink_model,new_equ,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            reaction.name=index
            reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(index)+'<hr/>'+\
                    'Reaction subsystem: <br/>'+'none'+'<hr/>'+'Reaction equation: <br/>'+\
                    str(new_equ)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+'none'+'<hr/>'

            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                #reaction.notes['map_info']['display_name'] = reaction.id + " "+ str("%.2f" %(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f" %(abs(row[col_flux])))
            elif rxn_switch == 'F':
                #reaction.notes['map_info']['display_name'] = reaction.name + " "+str("%.2f"%(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f"%(abs(row[col_flux])))
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])


            reaction.notes['map_info']['group']=3
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(new_equ) 

    #去除H
    H_list=['h_c','h_p','h_e','C00080','C00080[c]', 'C00080[p]','C00080[e]','PROTON', 'MNXM1', 'cpd00067', 'cpd00067_c0','cpd00067_p0','cpd00067_e0', 'HMDB59597', 'GPRLSGONYQIRFK-UHFFFAOYSA-N']
    for each_h in H_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #去除H2O
    H2O_list=['h2o_c','h2o_p','h2o_e','C00001', 'C01328','C00001[c]','C00001[p]','C00001[e]', 'C01328[c]','C01328[p]','C01328[e]', 'CPD-15815', 'HYDROXYL-GROUP', 'WATER', 'MNXM2', 'cpd00001', 'cpd00001_c0', 'cpd00001_p0','cpd00001_e0','cpd15275', 'cpd27222', 'HMDB01039', 'HMDB02111', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    for each_h2o in H2O_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h2o)
        except:
            pass
        else:
            hc_id.remove_from_model()
        #去除Pi
    Pi_list=['pi_c','pi_p','pi_e','C00009', 'C00009[c]','C00009[p]','C00009[e]','CPD-16459', 'PHOSPHATE-GROUP', 'Pi', 'MNXM9', 'cpd00009', 'cpd00009_c0', 'cpd00009_p0','cpd00009_e0','cpd27787', 'HMDB00973', 'HMDB01429', 'HMDB02105', 'HMDB02142', 'HMDB05947', 'NBIIXXVUZAFLBC-UHF']
    for each_pi in Pi_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_pi)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #print(len(metlink_model.reactions.get_by_id('DM_pyr_c').products))  
    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))
    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)
  
    #print(substrate_id,product_id)
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
            
        eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__D_e'     
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+''+'<hr/>'+\
                        'Metabolite formula: <br/>'+''+'<hr/>'
        eachmet.notes['map_info'] = {}
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A" 

        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])
        if eachmet.id.endswith('_c') or eachmet.id.endswith('_C') or eachmet.id.endswith('_p') or eachmet.id.endswith('_P') or eachmet.id.endswith('_e') or eachmet.id.endswith('_E'):
            pass
        else:
            eachmet.id=eachmet.id+'_c'
    return metlink_model

def flux2map(flux_solution_df,rxn_switch,coordinates,cofactor_switch,id_or_name,flow,col_equ,col_flux,html_file):
    exclude_rxn = ['EX_co2_e','H2Otex','O2tpp','CO2tpp','O2tex','NH4tex','CO2tex','NH4tpp','H2Otpp','EX_o2_e','EX_h2o_e','EX_nh4_e','EX_o2_e','ADD_h_c','ADD_h_p','ADD_respiratory1']
    total_common_factor_list =  ['hco3_c','ni2_c','cdp_c', 'ag_c', 'dctp_c', 'dutp_c', 'ctp_c', 'gdp_c', 'gtp_c', 'ump_c', 'ca2_c', \
                                'h2o_c', 'datp_c', 'co2_c', 'no2_c', 'no_c', 'k_c', 'zn2_c', 'no3_c', 'o2_c', 'cl_c', 'udp_c', 'damp_c',\
                                'ditp_c', 'dump_c', 'q8h2_c', 'pppi_c', 'idp_c', 'dimp_c', 'pi_c', 'dttp_c', 'so4_c', 'adp_c', 'xtp_c',\
                                'dgtp_c', 'dadp_c', 'coa_c', 'ppi_c', 'h2_c', 'cmp_c', 'fe2_c', 'o2s_c', 'h_c', 'gmp_c', 'itp_c', 'q8_c', \
                                'cobalt2_c', 'n2o_c', 'xmp_c', 'xdp_c', 'nadph_c', 'cu_c', 'cu2_c', 'atp_c', 'dgmp_c', 'imp_c', 'h2s_c', 'utp_c',\
                                'dtmp_c', 'fadh2_c', 'so3_c', 'fad_c', 'cd2_c', 'dgdp_c', 'nad_c', 'nadh_c', 'hg2_c', 'dcmp_c', 'dudp_c', 'dtdp_c',\
                                'didp_c', 'mn2_c', 'dcdp_c', 'nh4_c', 'amp_c', 'fe3_c', 'nadp_c', 'so2_c', 'h2o2_c', 'mg2_c',\
                                    'hco3_p','ni2_p','cdp_p', 'ag_p', 'dctp_p', 'dutp_p', 'ctp_p', 'gdp_p', 'gtp_p', 'ump_p', 'ca2_p', \
                                'h2o_p', 'datp_p', 'co2_p', 'no2_p', 'no_p', 'k_p', 'zn2_p', 'no3_p', 'o2_p', 'cl_p', 'udp_p', 'damp_p',\
                                'ditp_p', 'dump_p', 'q8h2_p', 'pppi_p', 'idp_p', 'dimp_p', 'pi_p', 'dttp_p', 'so4_p', 'adp_p', 'xtp_p',\
                                'dgtp_p', 'dadp_p', 'coa_p', 'ppi_p', 'h2_p', 'cmp_p', 'fe2_p', 'o2s_p', 'h_p', 'gmp_p', 'itp_p', 'q8_p', \
                                'cobalt2_p', 'n2o_p', 'xmp_p', 'xdp_p', 'nadph_p', 'cu_p', 'cu2_p', 'atp_p', 'dgmp_p', 'imp_p', 'h2s_p', 'utp_p',\
                                'dtmp_p', 'fadh2_p', 'so3_p', 'fad_p', 'cd2_p', 'dgdp_p', 'nad_p', 'nadh_p', 'hg2_p', 'dcmp_p', 'dudp_p', 'dtdp_p',\
                                'didp_p', 'mn2_p', 'dcdp_p', 'nh4_p', 'amp_p', 'fe3_p', 'nadp_p', 'so2_p', 'h2o2_p', 'mg2_p',\
                                    'hco3_e','ni2_e','cdp_e', 'ag_e', 'dctp_e', 'dutp_e', 'ctp_e', 'gdp_e', 'gtp_e', 'ump_e', 'ca2_e', \
                                'h2o_e', 'datp_e', 'co2_e', 'no2_e', 'no_e', 'k_e', 'zn2_e', 'no3_e', 'o2_e', 'cl_e', 'udp_e', 'damp_e',\
                                'ditp_e', 'dump_e', 'q8h2_e', 'pppi_e', 'idp_e', 'dimp_e', 'pi_e', 'dttp_e', 'so4_e', 'adp_e', 'xtp_e',\
                                'dgtp_e', 'dadp_e', 'coa_e', 'ppi_e', 'h2_e', 'cmp_e', 'fe2_e', 'o2s_e', 'h_e', 'gmp_e', 'itp_e', 'q8_e', \
                                'cobalt2_e', 'n2o_e', 'xmp_e', 'xdp_e', 'nadph_e', 'cu_e', 'cu2_e', 'atp_e', 'dgmp_e', 'imp_e', 'h2s_e', 'utp_e',\
                                'dtmp_e', 'fadh2_e', 'so3_e', 'fad_e', 'cd2_e', 'dgdp_e', 'nad_e', 'nadh_e', 'hg2_e', 'dcmp_e', 'dudp_e', 'dtdp_e',\
                                'didp_e', 'mn2_e', 'dcdp_e', 'nh4_e', 'amp_e', 'fe3_e', 'nadp_e', 'so2_e', 'h2o2_e', 'mg2_e', 'flxso_c','flxr_c','WATER_c',\
                                'PROTON_c','NADH_c','NAD_c','NADPH_c','NADP_c','ATP_c','ADP_c','Pi_c','pqq_c','pqqh2_c']  
    
    metlink_model= build_model_from_local_file(flux_solution_df,rxn_switch,coordinates,col_equ,col_flux)

    type_list=[]
    for eachr in metlink_model.reactions:
        type_list.append(eachr.notes['map_info']['group'])
    uni_type_list=list(set(type_list))

    # 可视化
    All_met, All_rxn= get_model_information(metlink_model)
    #print(All_rxn)
    common_factor_list=list(set(total_common_factor_list).intersection(set(All_met)))          

    if cofactor_switch == "T":
        d3flux.update_cofactors(metlink_model, common_factor_list)
        
    html = d3flux.flux_map(metlink_model, overwrite_reversibility=True,
                            excluded_metabolites = common_factor_list,excluded_reactions =exclude_rxn,
                            figsize=(1500, 1500), display_name_format=id_or_name,
                            hide_unused="true", hide_unused_cofactors="true",
                            flux_dict=flux_solution_df['fluxes'],
                            flowLayout = flow)  # 是否要隐藏flux为0的

    with open(html_file, "w")as wf:
        wf.write(html.data) 

    return html

def get_substrate_product_list_from_model(norm_model,substrate_use_file,product_use_file):
    product_use_dict={}
    substrate_use_dict={}
    #substrate_name_list=[]
    EX_exclude_reaction_list=['EX_pi_e','EX_h_e','EX_fe3_e','EX_mn2_e','EX_co2_e','EX_fe2_e','EX_h2_e','EX_zn2_e',\
                                    'EX_mg2_e','EX_ca2_e','EX_so3_e','EX_ni2_e','EX_no_e','EX_cu2_e','EX_hg2_e','EX_cd2_e',\
                                    'EX_h2o2_e','EX_h2o_e','EX_no2_e','EX_nh4_e','EX_so4_e','EX_k_e','EX_na1_e','EX_o2_e',\
                                    'EX_o2s_e','EX_ag_e','EX_cu_e','EX_so2_e','EX_cl_e','EX_n2o_e','EX_cs1_e','EX_cobalt2_e']

    EX_exclude_reaction_list=EX_exclude_reaction_list+[i+'_reverse' for i in EX_exclude_reaction_list]
    [substrate_use_dict,substrate_name_list]=exchanges_substrate(norm_model,substrate_use_dict)
    for eachm in norm_model.metabolites: 
        if eachm.id not in substrate_use_dict.keys():
            try:
                new_reacion=norm_model.add_boundary(norm_model.metabolites.get_by_id(eachm.id), type="sink")
            except:
                r_id='SK_'+eachm.id
                try:
                    r_name=norm_model.reactions.get_by_id(r_id).name.replace('sink','').replace('exchange','').replace("\n","")
                except:
                    pass
                else:
                    if r_name not in substrate_name_list:
                        substrate_use_dict[r_id]={}
                        substrate_use_dict[r_id]['substrate_id']=eachm.id
                        substrate_use_dict[r_id]['substrate_model_correspondind']=r_id
                        substrate_use_dict[r_id]['substrate_name']=r_name+'('+eachm.id+')'
                        substrate_use_dict[r_id]['type']="sink"
            else:
                r_name=new_reacion.name.replace('sink','').replace("\n","")
                #print(r_name)
                if r_name not in substrate_name_list:
                    substrate_use_dict[new_reacion.id]={}
                    substrate_use_dict[new_reacion.id]['substrate_id']=eachm.id
                    substrate_use_dict[new_reacion.id]['substrate_model_correspondind']=new_reacion.id
                    substrate_use_dict[new_reacion.id]['substrate_name']=r_name+'('+eachm.id+')'
                    substrate_use_dict[new_reacion.id]['type']="sink"
        if re.search('_c$',eachm.id):
            #print(eachm.id)
            EX_id = 'EX_'+eachm.id         
            if EX_id not in product_use_dict.keys():
                DE_id = 'user_add_DM_'+eachm.id  
                product_use_dict[DE_id]={}
                product_use_dict[DE_id]['obj_id']=eachm.id 
                product_use_dict[DE_id]['obj_model_correspondind']=DE_id
                product_use_dict[DE_id]['obj_name']=eachm.name.replace('exchange','').replace("\n","")+'('+eachm.id+')'
                product_use_dict[DE_id]['obj_reaction']='DM_'+eachm.id
                product_use_dict[DE_id]['type']="demand"
        elif re.search('_p$',eachm.id):
            pass
        elif re.search('_e$',eachm.id):
            pass   
        else:
            EX_id = 'EX_'+eachm.id         
            if EX_id not in product_use_dict.keys():
                DE_id = 'user_add_DM_'+eachm.id  
                product_use_dict[DE_id]={}
                product_use_dict[DE_id]['obj_id']=eachm.id 
                product_use_dict[DE_id]['obj_model_correspondind']=DE_id
                product_use_dict[DE_id]['obj_name']=eachm.name.replace('exchange','').replace("\n","")+'('+eachm.id+')'
                product_use_dict[DE_id]['obj_reaction']='DM_'+eachm.id
                product_use_dict[DE_id]['type']="demand"
 
    product_use_dict['Object']={}
    obj_id = ''
    regex = re.compile('BIOMASS', flags=re.IGNORECASE)

    try:
        objective = norm_model.objective.expression
        reaction_id = list(norm_model.objective.variables)[0].name
    except:
        if len(norm_model.reactions.query(regex))>1:
            for eglc in norm_model.reactions.query(regex):
                with norm_model as model:
                    model.objective=eglc.id
                    try:
                        opt_r=model.optimize()
                    except:
                        pass
                    else:
                        if opt_r[eglc.id]>0:
                            if len(obj_id)<len(eglc.id):
                                obj_id = eglc.id

            product_use_dict['Object']['obj_id']=obj_id
            product_use_dict['Object']['obj_model_correspondind']=obj_id
            product_use_dict['Object']['obj_name']='Biomass'
            product_use_dict['Object']['obj_reaction']=obj_id
            product_use_dict['Object']['type']="default"
        else:
            for eglc in norm_model.reactions.query(regex):
                product_use_dict['Object']['obj_id']=eglc.id
                product_use_dict['Object']['obj_model_correspondind']=eglc.id
                product_use_dict['Object']['obj_name']='Biomass'
                product_use_dict['Object']['obj_reaction']=eglc.id
                product_use_dict['Object']['type']="default"
    else:
        try:
            reaction = norm_model.reactions.get_by_id(reaction_id)
        except:
            reaction_id = list(norm_model.objective.variables)[1].name
            reaction = norm_model.reactions.get_by_id(reaction_id)
            
        product_use_dict['Object']['obj_id']=reaction.id.replace('EX_','')
        product_use_dict['Object']['obj_model_correspondind']=reaction.id
        product_use_dict['Object']['obj_name']=reaction.name
        product_use_dict['Object']['obj_reaction']=reaction.id
        product_use_dict['Object']['type']="default"   

    substrate_use_list=pd.DataFrame(substrate_use_dict).T
    product_use_list=pd.DataFrame(product_use_dict).T
    # 保存substrate、product
    #substrate_use_file = './%s_substrate_use_list.tsv'%model_name
    #product_use_file  = './%s_product_use_list.tsv'%model_name
    substrate_use_list.to_csv(substrate_use_file,sep='\t',index=False)   
    product_use_list.to_csv(product_use_file,sep='\t',index=False)
    
def get_substrate_product_list_from_mcmodel(model,substrate_use_file,product_use_file):
    product_use_dict={}
    substrate_use_dict={}

    substrate_id = list(model.medium.keys())[0]
    substrate_use_dict[substrate_id]={}
    r_infor = model.get_reaction_by_id(substrate_id)
    substrate_use_dict[substrate_id]['substrate_id']=list(r_infor.metabolites.keys())[0]
    substrate_use_dict[substrate_id]['substrate_model_correspondind']=substrate_id                  
    substrate_use_dict[substrate_id]['substrate_name']=model.get_reaction_by_id(substrate_id).name
    substrate_use_dict[substrate_id]['type']="default"  
    for metabolite in model.metabolites.keys():
        eachm = model.get_metabolite_by_id(metabolite)
        if eachm.id not in substrate_use_dict.keys():
            try:
                r_id='SK_'+eachm.id
                model.add_sink_reaction(eachm.id)
                new_reacion=model.get_reaction_by_id(r_id)
            except:
                try:
                    r_name=model.get_reaction_by_id(r_id).name.replace('sink','').replace('exchange','').replace("\n","")
                except:
                    pass
                else:
                    substrate_use_dict[r_id]={}
                    substrate_use_dict[r_id]['substrate_id']=eachm.id
                    substrate_use_dict[r_id]['substrate_model_correspondind']=r_id
                    substrate_use_dict[r_id]['substrate_name']=r_name+'('+eachm.id+')'
                    substrate_use_dict[r_id]['type']="sink"
            else:
                r_name=eachm.name
                substrate_use_dict[new_reacion.id]={}
                substrate_use_dict[new_reacion.id]['substrate_id']=eachm.id
                substrate_use_dict[new_reacion.id]['substrate_model_correspondind']=new_reacion.id
                substrate_use_dict[new_reacion.id]['substrate_name']=r_name+'('+eachm.id+')'
                substrate_use_dict[new_reacion.id]['type']="sink"
        if re.search('_c$',eachm.id):
            #print(eachm.id)
            EX_id = 'EX_'+eachm.id         
            if EX_id not in product_use_dict.keys():
                DE_id = 'user_add_DM_'+eachm.id  
                product_use_dict[DE_id]={}
                product_use_dict[DE_id]['obj_id']=eachm.id 
                product_use_dict[DE_id]['obj_model_correspondind']=DE_id
                product_use_dict[DE_id]['obj_name']=eachm.name.replace('exchange','').replace("\n","")+'('+eachm.id+')'
                product_use_dict[DE_id]['obj_reaction']='DM_'+eachm.id
                product_use_dict[DE_id]['type']="demand"
        elif re.search('_p$',eachm.id):
            pass
        elif re.search('_e$',eachm.id):
            pass   
        else:
            EX_id = 'EX_'+eachm.id         
            if EX_id not in product_use_dict.keys():
                DE_id = 'user_add_DM_'+eachm.id  
                product_use_dict[DE_id]={}
                product_use_dict[DE_id]['obj_id']=eachm.id 
                product_use_dict[DE_id]['obj_model_correspondind']=DE_id
                product_use_dict[DE_id]['obj_name']=eachm.name.replace('exchange','').replace("\n","")+'('+eachm.id+')'
                product_use_dict[DE_id]['obj_reaction']='DM_'+eachm.id
                product_use_dict[DE_id]['type']="demand"
 
    product_use_dict['Object']={}     
    product_use_dict['Object']['obj_id']=model.obj_name
    product_use_dict['Object']['obj_model_correspondind']=model.obj_name
    product_use_dict['Object']['obj_name']=model.get_reaction_by_id(model.obj_name).name
    product_use_dict['Object']['obj_reaction']=model.obj_name
    product_use_dict['Object']['type']="default"   

    substrate_use_list=pd.DataFrame(substrate_use_dict).T
    product_use_list=pd.DataFrame(product_use_dict).T
    # 保存substrate、product
    #substrate_use_file = './%s_substrate_use_list.tsv'%model_name
    #product_use_file  = './%s_product_use_list.tsv'%model_name
    substrate_use_list.to_csv(substrate_use_file,sep='\t',index=False)   
    product_use_list.to_csv(product_use_file,sep='\t',index=False)

def create_mcmodel_change_json(model):
    # 生成model_change
    model_change_json={}
    model_change_json['reaction_infor']={}
    model_change_json['gene_infor']={}
    model_change_json['exchange_infor']={}  #add_reaction后 model.exchanges没了
    model_change_json['bigg_metabolite_ids']={}  
    model_change_json['add_reaction_infor']={}    # 与modify生成的model_change_json_final.keys()保持一致
    model_change_json['del_reaction_infor']={}
    model_change_json['del_gene_infor']={}
    for reaction in model.reactions.keys():
        eachr = model.get_reaction_by_id(reaction)
        model_change_json['reaction_infor'][eachr.id]={}
        model_change_json['reaction_infor'][eachr.id]['reaction_name']=eachr.name
        model_change_json['reaction_infor'][eachr.id]['equation']=eachr.build_reaction_string()
        model_change_json['reaction_infor'][eachr.id]['equation by name']=eachr.build_reaction_string().replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
        model_change_json['reaction_infor'][eachr.id]['gene_reaction_rule']=eachr.gene_reaction_rule
        model_change_json['reaction_infor'][eachr.id]['lower_bound']=eachr.lower_bound
        model_change_json['reaction_infor'][eachr.id]['upper_bound']=eachr.upper_bound
        model_change_json['reaction_infor'][eachr.id]['g0']=eachr.g0
        model_change_json['reaction_infor'][eachr.id]['mw']=eachr.mw
        model_change_json['reaction_infor'][eachr.id]['kcat']=eachr.kcat
    for gene in model.genes.keys():
        eachg = model.get_gene_by_id(gene)
        if eachg.name:
            model_change_json['gene_infor'][eachg.id]={}
            model_change_json['gene_infor'][eachg.id]['gene_name']=eachg.name
            reaction_id_list=[]
            reaction_name_list=[]
            reaction_equ_list=[]
            for reaction in model.reactions.keys():
                eachr = model.get_reaction_by_id(reaction)
                if eachg.id in eachr.gene_reaction_rule:
                    reaction_name_list.append(eachr.name)
                    reaction_equ_list.append(eachr.build_reaction_string())
            model_change_json['gene_infor'][eachg.id]['reaction_name_list']='; '.join(reaction_name_list) 
            model_change_json['gene_infor'][eachg.id]['reaction_equ_list']='; '.join(reaction_equ_list)    
    for metabolite in model.metabolites.keys():
        eachm = model.get_metabolite_by_id(metabolite)
        model_change_json['bigg_metabolite_ids'][eachm.id]={}
        model_change_json['bigg_metabolite_ids'][eachm.id]['bigg_id']=eachm.id
        model_change_json['bigg_metabolite_ids'][eachm.id]['name']=eachm.name
    return model_change_json

def get_preview_default_json_mcmodel(model,met_C):
    default_setting_json={}

    default_setting_json['Oxygen']={}
    default_setting_json['Oxygen']['model_id']='None'
    default_setting_json['Oxygen']['model_correspondind']='None'
    default_setting_json['Oxygen']['value']=0 
    default_setting_json['Carbon_source']={}
    default_setting_json['Carbon_source']['substrate_id']='None'
    default_setting_json['Carbon_source']['substrate_model_correspondind']='None'
    default_setting_json['Carbon_source']['substrate_name']='None'
    default_setting_json['Carbon_source']['uptake_rate']=0 
    default_setting_json['Nitrogen_source']={}
    default_setting_json['Nitrogen_source']['substrate_id']='None'
    default_setting_json['Nitrogen_source']['substrate_model_correspondind']='None'
    default_setting_json['Nitrogen_source']['substrate_name']='None' 
    default_setting_json['Nitrogen_source']['uptake_rate']=0 
    default_setting_json['ATPM']={}
    default_setting_json['ATPM']['model_id']='None'
    default_setting_json['ATPM']['model_correspondind']='None'
    default_setting_json['ATPM']['value']=0 
    default_setting_json['Object']={}
    default_setting_json['Object']['model_id']=model.obj_name
    default_setting_json['Object']['model_correspondind']=model.obj_name
    default_setting_json['Object']['substrate_name']= model.get_reaction_by_id(model.obj_name).name

    O2_list=['EX_o2_e_reverse','EX_o2_e','EX_cpd00007_b','EX_cpd00007_e0','EX_o2(e)','EX_MNXM735438@BOUNDARY','mnxr102090c2b','Ex_o2','Ex_o2_reverse','R1287']
    #ATPM_list=['ATPM','rxn00062_c0','DM_atp_c_','mnxr153054c3']
    ATPM_list=['ATPM','rxn00062','rxn00062_c0','DM_atp_c_','mnxr153054c3','R00086','META:ATPASE-RXN','ATPASE-RXN','MNXR96131','rxn09694', 'rxn11300','rxn09694_c0', 'rxn11300_c0','ATPM_c', 'ATPM_LPAREN_NGAM_RPAREN_', 'DM_atp_c', 'DM_atp[c]']

    #关闭碳源
    default_C_list=[]
    exclude_C_reaction_list=['EX_BIOTIN_reverse','EX_BIOTIN','EX_biotin_e', 'EX_btn', 'EX_btn_e', 'EX_btn_e_', 'EX_btn(e)', 'EX_btn[e]', 'EX_btn_LPAREN_e_RPAREN_','EX_biotin_e_reverse', 'EX_btn_reverse', 'EX_btn_e_reverse', 'EX_btn_e__reverse', 'EX_btn(e)_reverse', 'EX_btn[e]_reverse', 'EX_btn_LPAREN_e_RPAREN__reverse','META:TRANS-RXN0-240','MNXR96333']
    for eachr in model.reactions.keys():
        if "EX_" in eachr:
            #print(eachr)
            #r_bounds=norm_model.reactions.get_by_id(eachr).bounds
            product_equ=model.get_reaction_by_id(eachr)
            # 初始化计数器
            count_negative = 0
            # 迭代字典的值
            for value in product_equ.metabolites.values():
                if value < 0:
                    count_negative += 1
            if count_negative>1:
                pass
            else:
                for eachm in product_equ.metabolites.keys():
                    if product_equ.metabolites[eachm]<0:
                        #print(eachr,eachm.id,eachm.formula)
                        #met_coef=abs(product_equ.get_coefficient(eachm.id))
                        try:
                            met_formula=get_composition_from_formula(model.get_metabolite_by_id(eachm).formula)
                        except:
                            try:
                                met_formula=get_composition_from_formula(met_C.loc[eachm,'formula'])
                                #print(eachr,eachm.id,met_formula)s
                            except:
                                print('No formula!')
                            else:
                                if 'C' in met_formula.keys():#是碳源才关闭
                                    if abs(model.get_reaction_by_id(eachr).lower_bound) <1000 and abs(model.get_reaction_by_id(eachr).lower_bound) >0:
                                        if 'Carbon_source' in default_setting_json.keys():
                                            if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).lower_bound):
                                                default_setting_json['Carbon_source']['substrate_id']=eachm
                                                default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                                default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                                default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)             
                                        else:
                                            default_setting_json['Carbon_source']={}
                                            #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                   
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)
                                    elif model.get_reaction_by_id(eachr).upper_bound>0 and model.get_reaction_by_id(eachr).upper_bound<1000:
                                        if 'Carbon_source' in default_setting_json.keys():
                                            if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).upper_bound):
                                                default_setting_json['Carbon_source']['substrate_id']=eachm
                                                default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                    
                                                default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                                default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)             
                                        else:
                                            default_setting_json['Carbon_source']={}
                                            #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                    
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)  
                                    
                                if 'N' in met_formula.keys():#是
                                    default_setting_json['Nitrogen_source']['substrate_id']=eachm
                                    default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachr                   
                                    default_setting_json['Nitrogen_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                    if model.get_reaction_by_id(eachr).lower_bound<0:
                                        default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound) 
                                    else:
                                        default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)      
                                                                           
                        else:
                            if 'C' in met_formula.keys():#是碳源才关闭
                                if abs(model.get_reaction_by_id(eachr).lower_bound) <1000 and abs(model.get_reaction_by_id(eachr).lower_bound) >0:
                                    if 'Carbon_source' in default_setting_json.keys():
                                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).lower_bound):
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                    
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)             
                                    else:
                                        default_setting_json['Carbon_source']={}
                                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                        default_setting_json['Carbon_source']['substrate_id']=eachm
                                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                    
                                        default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                        default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)
                                elif model.get_reaction_by_id(eachr).upper_bound>0 and model.get_reaction_by_id(eachr).upper_bound<1000:
                                    if 'Carbon_source' in default_setting_json.keys():
                                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).upper_bound):
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                   
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)             
                                    else:
                                        default_setting_json['Carbon_source']={}
                                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                        default_setting_json['Carbon_source']['substrate_id']=eachm
                                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                    
                                        default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                        default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound) 
                            if 'N' in met_formula.keys():#是
                                default_setting_json['Nitrogen_source']['substrate_id']=eachm
                                default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachr                   
                                default_setting_json['Nitrogen_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                if model.get_reaction_by_id(eachr).lower_bound<0:
                                    default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound) 
                                else:
                                    default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound) 

            count_positive = 0
            # 迭代字典的值
            for value in product_equ.metabolites.values():
                if value > 0:
                    count_positive += 1
            if count_positive>1:
                pass
            else:
                for eachm in product_equ.metabolites.keys():
                    if product_equ.metabolites[eachm]>0:
                    #print(eachr,eachm.id,eachm.formula)
                    #met_coef=abs(product_equ.get_coefficient(eachm.id))
                        try:
                            met_formula=get_composition_from_formula(model.get_metabolite_by_id(eachm).formula)
                        except:
                            try:
                                met_formula=get_composition_from_formula(met_C.loc[eachm,'formula'])
                                #print(eachr,eachm.id,met_formula)
                            except:
                                print('No formula!')
                            else:
                                if 'C' in met_formula.keys():#是碳源才关闭
                                    if abs(model.get_reaction_by_id(eachr).lower_bound) <1000 and abs(model.get_reaction_by_id(eachr).lower_bound) >0:
                                        if 'Carbon_source' in default_setting_json.keys():
                                            if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).lower_bound):
                                                default_setting_json['Carbon_source']['substrate_id']=eachm
                                                default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                                default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                                default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)             
                                        else:
                                            default_setting_json['Carbon_source']={}
                                            #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)
                                    elif model.get_reaction_by_id(eachr).upper_bound>0 and model.get_reaction_by_id(eachr).upper_bound<1000:
                                        if 'Carbon_source' in default_setting_json.keys():
                                            if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).upper_bound):
                                                default_setting_json['Carbon_source']['substrate_id']=eachm
                                                default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                    
                                                default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                                default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)             
                                        else:
                                            default_setting_json['Carbon_source']={}
                                            #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)     
                                if 'N' in met_formula.keys():#是
                                    default_setting_json['Nitrogen_source']['substrate_id']=eachm
                                    default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachr                   
                                    default_setting_json['Nitrogen_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                    if model.get_reaction_by_id(eachr).lower_bound<0:
                                        default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound) 
                                    else:
                                        default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)                           
                        else:
                            if 'C' in met_formula.keys():#是碳源才关闭
                                if abs(model.get_reaction_by_id(eachr).lower_bound) <1000 and abs(model.get_reaction_by_id(eachr).lower_bound) >0:
                                    if 'Carbon_source' in default_setting_json.keys():
                                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).lower_bound):
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                   
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)             
                                    else:
                                        default_setting_json['Carbon_source']={}
                                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                        default_setting_json['Carbon_source']['substrate_id']=eachm
                                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                        default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                        default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound)
                                elif model.get_reaction_by_id(eachr).upper_bound>0 and model.get_reaction_by_id(eachr).upper_bound<1000:
                                    if 'Carbon_source' in default_setting_json.keys():
                                        if default_setting_json['Carbon_source']['uptake_rate'] < abs(model.get_reaction_by_id(eachr).upper_bound):
                                            default_setting_json['Carbon_source']['substrate_id']=eachm
                                            default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                            default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                            default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)             
                                    else:
                                        default_setting_json['Carbon_source']={}
                                        #default_setting_json['Carbon_source']['substrate_model_correspondind']='EX_'+reactant_id
                                        default_setting_json['Carbon_source']['substrate_id']=eachm
                                        default_setting_json['Carbon_source']['substrate_model_correspondind']=eachr                     
                                        default_setting_json['Carbon_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                        default_setting_json['Carbon_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)    
                            if 'N' in met_formula.keys():#是
                                default_setting_json['Nitrogen_source']['substrate_id']=eachm
                                default_setting_json['Nitrogen_source']['substrate_model_correspondind']=eachr                   
                                default_setting_json['Nitrogen_source']['substrate_name']=model.get_metabolite_by_id(eachm).name
                                if model.get_reaction_by_id(eachr).lower_bound<0:
                                    default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).lower_bound) 
                                else:
                                    default_setting_json['Nitrogen_source']['uptake_rate']=abs(model.get_reaction_by_id(eachr).upper_bound)                                                      
    #ATPM
    for eachatpm in ATPM_list:
        try:                
            ATPM=model.get_reaction_by_id(eachatpm).lower_bound
            print(eachatpm,ATPM)
        except:
            print('Model does not contain ATP Maintenance Requirement')
            #default_setting_json['ATPM']='none'
        else:
            default_setting_json['ATPM']={}
            default_setting_json['ATPM']['model_id']=eachatpm
            default_setting_json['ATPM']['model_correspondind']=eachatpm
            default_setting_json['ATPM']['substrate_name']='ATPM'
            default_setting_json['ATPM']['value']=ATPM 
            break

    #有氧厌氧
    for eachatpm in O2_list:
        try:
            O2=model.get_reaction_by_id(eachatpm).upper_bound
        except:
            print('No O2 reaction!')         
        else:
            default_setting_json['Oxygen']={}
            default_setting_json['Oxygen']['model_id']=eachatpm
            default_setting_json['Oxygen']['model_correspondind']=eachatpm
            default_setting_json['Oxygen']['substrate_name']='Oxygen'
            default_setting_json['Oxygen']['value']=abs(O2)
            break
               
    return default_setting_json