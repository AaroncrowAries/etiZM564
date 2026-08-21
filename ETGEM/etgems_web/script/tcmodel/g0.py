from ..core.model import Model 
from ..core.environment import Environment 
from ..core.reaction import Reaction
from ..core.metabolite import Metabolite
from ..core.gene import Gene
from equilibrator_api import ComponentContribution, Q_

def get_reaction_g0(model):
    """get ΔG'° use equilibrator_api.
    Arguments
    ----------
    * model: the dictionary of model
    * p_h: the pH of the environment
    * p_mg: the pMg of the environment
    * ionic_strength: the ionic strength of the environment
    * temperature: the temperature of the environment
    """
    p_h = model.environment.p_h
    p_mg = model.environment.p_mg
    ionic_strength = model.environment.ionic_strength
    temperature = model.environment.temperature
    cc = None
    while cc is None:
        try:
            cc = ComponentContribution()
        except JSONDecodeError:
            logger.warning('Waiting for zenodo.org... Retrying in 5s')
            sleep(5)

    #cc = ComponentContribution()
    cc.p_h = Q_(p_h)
    cc.p_mg = Q_(p_mg)
    cc.ionic_strength = Q_(ionic_strength)
    cc.temperature = Q_(temperature)
    # reaction_g0={}
    for eachr in model.reactions.values():
        if not eachr.id.startswith("EX_") and not eachr.id.startswith("DM_"):
            reaction_left=[]
            reaction_right=[]
            for k, v in eachr.metabolites.items():
                if str(k).endswith('_c'):
                    k_new = "_c".join(str(k).split('_c')[:-1])
                elif str(k).endswith('_p'):
                    k_new = "_p".join(str(k).split('_p')[:-1])
                elif str(k).endswith('_e'):
                    k_new = "_e".join(str(k).split('_e')[:-1]) 
                else:
                    k_new = str(k)[:-2]
                #kegg,chebi,metanetx;kegg:C00002 + CHEBI:15377 = metanetx.chemical:MNXM7 + bigg.metabolite:pi
                if v<0:                    
                    reaction_left.append(str(-v)+' bigg.metabolite:'+k_new)
                else:
                    reaction_right.append(str(v)+' bigg.metabolite:'+k_new)
            reaction_equ=(' + ').join(reaction_left)+' -> '+(' + ').join(reaction_right)
            #get ΔG'° use equilibrator_api
            try:
                equilibrator_api_reaction = cc.parse_reaction_formula(reaction_equ)
                #print("The reaction is " + ("" if equilibrator_api_reaction.is_balanced() else "not ") + "balanced")
                dG0_prime = cc.standard_dg_prime(equilibrator_api_reaction)
            except:
                pass
                # print('error reaction: '+eachr.id)
                # print(reaction_equ)
            else:  
                # reaction_g0[eachr.id]={}
                # reaction_g0[eachr.id]['reaction']=eachr.id
                # reaction_g0[eachr.id]['equ']=eachr.reaction
                # #dG0 (kilojoule/mole)
                # reaction_g0[eachr.id]['g0']=str(dG0_prime).split(') ')[0].split('(')[1].split(' +/- ')[0]
                eachr.g0 = float(str(dG0_prime).split(') ')[0].split('(')[1].split(' +/- ')[0])
    return model