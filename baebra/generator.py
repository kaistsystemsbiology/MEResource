import re
import os
from copy import deepcopy
import logging

import pandas as pd
from cobra import Reaction, Metabolite
from cobra.io import read_sbml_model, write_sbml_model

from tqdm import tqdm


'''
Module 0. Extract information from the summarized data
'''

def extract_bigg_info(bigg_met_data):
    p = re.compile('(CHEBI:\d+);')

    chebi2info = {}
    with open(bigg_met_data, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if 'CHEBI:' in line:
                data = line.split('\t')
                met_id = data[1]
                name = data[2]
                chebi_ids = p.findall(line)
                for each_chebi in chebi_ids:
                    if each_chebi not in chebi2info:
                        chebi2info[each_chebi] = {'name':[], 'met_id':[]}
                    if name not in chebi2info[each_chebi]['name']:
                        chebi2info[each_chebi]['name'].append(name)
                    if met_id not in chebi2info[each_chebi]['met_id']:
                        chebi2info[each_chebi]['met_id'].append(met_id)
    print(f'Number of CHEBI IDs in BiGG database: {len(chebi2info)}')
    return chebi2info


def extract_stoich(stoich_data):
    with open(stoich_data, 'r') as fp:
        n = 0
        lines = fp.readlines()
        stoich_info = {}
        for i, line in enumerate(lines):
            ind = i%3
            if ind == 0:
                rhea_id = line.split('\t')[1]
            elif ind == 1:
                metabolites = line.strip().split('\t')[1:]
            elif ind == 2:
                stoichs = line.strip().split('\t')[1:]
                tmp = {}
                for met_id, coeff in zip(metabolites, stoichs):
                    try:
                        tmp[met_id] = float(coeff)
                    except:
                        print(met_id, coeff, line)
                stoich_info[rhea_id] = tmp
    return stoich_info


def extract_metabolite_info(met_data):
    table = pd.read_table(met_data, index_col='CHEBI_ID')
    return table 


def extract_reaction_info(rxn_data):
    table = pd.read_table(rxn_data, index_col='Rhea')
    return table



'''
Module 1. Find all metabolites which participate in reactions for target chemical production path
'''

def _get_path_reactions(target_chemical_line): # target_chemical_line: CHEBI:XXXX butan-1-ol acetyl-CoA R1 R2 R3
    data = target_chemical_line.strip().split('\t')
    target_chem = data[0]
    path_reactions = data[3:]
    return target_chem, path_reactions


def _extract_metabolites(reaction_id, stoich_info):
    met_stoichs = stoich_info[reaction_id]
    return met_stoichs.keys()


def _extract_metabolite_ID(met_id, bigg_met_info):
    if met_id in bigg_met_info:
        met_bigg_id = bigg_met_info[met_id]['met_id']
        if len(met_bigg_id) > 1:
            print(f'{met_id} has multiple BiGG metabolite ids')
        bigg_id = met_bigg_id[0]
    else:
        p = re.compile('CHEBI:(\w+)')

        try:
            bigg_id = 'met_%s'%(p.match(met_id).group(1))
        except:
            print(met_id)
            bgg

    return bigg_id
   

def _get_metabolite_object(met_id, met_bigg_id, cobra_model, met_info_table, extracellular=False):
    if extracellular:
        met_compart = 'e'
    else:
        met_compart = 'c'
    metabolite_id = f'{met_bigg_id}_{met_compart}'
    model_metabolites = [met.id for met in cobra_model.metabolites]
    if metabolite_id in model_metabolites:
        metabolite = cobra_model.metabolites.get_by_id(metabolite_id)
    else:
        metabolite_info = met_info_table.loc[met_id]
        formula = metabolite_info['Formula']
        charge = metabolite_info['Charge']
        mass = metabolite_info['Mass']
        
        if not formula:
            formula ='None'
        if not charge:
            charge = 0
        else:
            charge = int(charge)

        metabolite = Metabolite(metabolite_id,
                                name=metabolite_info['Name'],
                                formula=formula,
                                charge=charge,
                                compartment=met_compart
                               )
    return metabolite


def getMetabolites(target_data, stoich_info, bigg_met_info, cobra_model, met_info_table):
    metabolite_dict = {}
    with open(target_data, 'r') as fp:
        target_path_lines = fp.readlines()
    for target_chemical_line in target_path_lines[1:]:
        _, path_reactions = _get_path_reactions(target_chemical_line)
        for each_path_reaction in path_reactions:
            path_metabolites = _extract_metabolites(each_path_reaction, stoich_info)
            for each_met in path_metabolites:
                if each_met in metabolite_dict:
                    continue
                each_met_bigg_id = _extract_metabolite_ID(each_met, bigg_met_info)
                metabolite_dict[each_met] = _get_metabolite_object(each_met,
                                                                   each_met_bigg_id,
                                                                   cobra_model,
                                                                   met_info_table)
    return metabolite_dict


'''
Module 2. Create (or get) reactions in the target chemical production pathway
'''

def getReactions(target_data, stoich_info, cobra_model, metabolite_dict, rxn_info_table):
    model_reactions = [rxn.id for rxn in cobra_model.reactions]
    with open(target_data, 'r') as fp:
        target_path_lines = fp.readlines()
    path2targets = {}
    wrong_dir = {}
    for target_chemical_line in target_path_lines[1:]:
        target_chem, path_reactions = _get_path_reactions(target_chemical_line)
        if target_chem not in path2targets:
            path2targets[target_chem] = []
        path_reaction_objects = []
        for each_path_reaction in path_reactions:
            try:
                reaction_bigg_id = rxn_info_table.loc[each_path_reaction]['BiGG_ID']
            except:
                print(target_chem, path_reactions)
                
            reaction_info = rxn_info_table.loc[each_path_reaction]
            direction = reaction_info['Direction']
            reaction_name = reaction_info['Enzyme']
            if direction == '<=>':
                rhea_dir = (-1000, 1000)
            elif direction == '>':
                rhea_dir = (0, 1000)
            else:
                print('Reaction %s: Wrong direction (%s)'%(each_path_reaction, direction))
            if reaction_bigg_id in model_reactions:
                tmp_reaction = cobra_model.reactions.get_by_id(reaction_bigg_id)
                if rhea_dir != tmp_reaction.bounds:
                    if reaction_bigg_id not in wrong_dir:
                        wrong_dir[reaction_bigg_id] = []
                    if each_path_reaction not in wrong_dir[reaction_bigg_id]:
                        wrong_dir[reaction_bigg_id].append(each_path_reaction)      
                    tmp_reaction.bounds = tmp_reaction.bounds

                continue
                
            reaction = Reaction(reaction_bigg_id)
            reaction.name = reaction_name
            reaction.bounds = rhea_dir
            coeff_dict = {}
            tmp_coeff_dict = stoich_info[each_path_reaction]
            coeff_dict = {metabolite_dict[item]:tmp_coeff_dict[item] for item in tmp_coeff_dict}
            reaction.add_metabolites(coeff_dict)
            
            path_reaction_objects.append(reaction)
        path2targets[target_chem].append(path_reaction_objects)
        
    return path2targets, wrong_dir


'''
Module 3. Create exchange reaction for target chemicals
'''
def _make_external_rxn(met_bigg_id, target_chem_cytosol):
    extractRxn = Reaction(f'EXT_{met_bigg_id}_c')
    extractRxn.name = f'{target_chem_cytosol.name} to extracellular space'
    extractRxn.bounds = (0, 1000)
    extractRxn.add_metabolites({
        target_chem_cytosol:-1
    })
    return [extractRxn]


def _curate_model(met_bigg_id, model):
    curate_dict = {
        'ala_B': {'BALAt2pp': (-1000, 1000)},
        'but': {'BUTCT': (-1000, 1000)},
        'gam6p': {'GAM6Pt6_2pp': (-1000, 1000)},
        'gcald': {'GCALDD': (-1000, 1000)},
        'gln__L': {'GLNabcpp': (-1000, 1000)},
        'pser__L': {'PSP_Lpp': (-1000, 1000)},
        'octa': {'FACOAL80t2pp': (-1000, 1000)},
        'skm': {'SKMt2pp': (-1000, 1000)},
        'acetone': {'ACACCT': (-1000, 1000)},
        'met_30742': {'GCALDD': (-1000, 1000)},
        'itacon': {'CAD': (-1000, 1000)},
        '2ppoh': {'ACACCT': (-1000, 1000)},
        'ppa': {'PPAt4pp': (-1000, 1000)},
        'for': {'FORtex_a': (0, 1000)},
        'hxa': {'FACOAE60': (0, 1000)},
        'met_37830': {'FACOAE60': (0, 1000)},
        'asp__L': {'ASPtex_a': (0, 1000)},
        'glyc__R': {'GLYCK': (-1000, 1000)},
        'actn__R': {'ACTNtpp': (-1000, 1000)},
        '4hoxoh': {'HOPNTAL2': (-1000, 1000)},
    }
    
    if met_bigg_id not in curate_dict:
        return 
    
    model_reaction_list = [rxn.id for rxn in model.reactions]
    for item in curate_dict[met_bigg_id]:
        if item in model_reaction_list:
            model.reactions.get_by_id(item).bounds = curate_dict[met_bigg_id][item]
    return 


def _curate_reaction(add_reaction_list):
    curate_dict = {
        'PPPGP': (-1000, 1000),
        'G3PD5': (-1000, 1000),
        'GLYCK': (-1000, 1000),
        'r0393': (-1000, 1000),
        'PROD2': (-1000, 1000),
        'P5CR': (-1000, 1000),
        'ACGAMK': (-1000, 1000),
        'GLYCTO1': (-1000, 1000),
        'BUTCT': (-1000, 1000),
        'LCADi': (-1000, 1000),
    }
    
    for rxn in add_reaction_list:
        if rxn.id in curate_dict:
            rxn.bounds = curate_dict[rxn.id]
    return



def save_chemical_models(path_reaction_dict, bigg_met_info, input_model,
                         met_dict, met_info_table, output_dir, model_name):
    cobra_model = deepcopy(input_model)
    added_reactions = {}
    added_reactions_name = {}
    metabolite_list = [met.id for met in cobra_model.metabolites]
    reaction_list = [rxn.id for rxn in cobra_model.reactions]
    # pbar = tqdm(path_reaction_dict, position=0, leave=True)
    target_count = {}
    with tqdm(len(path_reaction_dict), position=0, leave=True) as pbar:
        for i, target_chem in enumerate(tqdm(path_reaction_dict, position=0, leave=True)):
            candidate_reactions_set = path_reaction_dict[target_chem]
            target_count[target_chem] = 0
            for j, candidate_reactions in enumerate(candidate_reactions_set):
                with cobra_model as model:
                    target_chem_object = met_dict[target_chem]
                    met_bigg_id = _extract_metabolite_ID(target_chem, bigg_met_info)
                    
                    # instead of using chem_c --> chem_e system, make a reaction which maeks chem --> None
                    candidate_reactions += _make_external_rxn(met_bigg_id, target_chem_object)
                
                    _curate_reaction(candidate_reactions)
                    
                    model.add_reactions(candidate_reactions)
                    target_chem_name = model.metabolites.get_by_id(met_bigg_id+'_c').name
                    
                    _curate_model(met_bigg_id, model)
                    
                    if met_bigg_id in added_reactions:
                        met_bigg_id += '_path_%d'%target_count[target_chem]
                            
                    added_reactions[met_bigg_id] = [added_rxn.id for added_rxn in candidate_reactions]
                    # added_reactions_name[target_chem_name] = [added_rxn.name for added_rxn in candidate_reactions]
                    added_reactions_name[model_name+'_'+met_bigg_id] = [target_chem_name]+[added_rxn.name for added_rxn in candidate_reactions]

                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    if output_dir[-1] == '/':
                        output_directory = f'{output_dir}/{model_name}_{met_bigg_id}.xml'
                    else:
                        output_directory = f'{output_dir}/{model_name}_{met_bigg_id}.xml'
                        
                    write_sbml_model(model, output_directory)
                    target_count[target_chem] += 1
                    
            pbar.set_description(f'Model number\tTarget chemical: {met_bigg_id}')
    added_table = pd.DataFrame.from_dict(added_reactions, orient='index')
    added_table_name = pd.DataFrame.from_dict(added_reactions_name, orient='index')
    return added_table, added_table_name




def save_chemical_model(target_chem, chem_paths, chebi2info, input_model, met_dict):

    added_reactions = {}
    added_reactions_name = {}
    # pbar = tqdm(path_reaction_dict, position=0, leave=True)
    target_count = {}
    target_models = {}

    for i, candidate_reactions in enumerate(chem_paths):
        model = deepcopy(input_model)
        target_count[target_chem] = i
        target_chem_object = met_dict[target_chem]
        met_bigg_id = _extract_metabolite_ID(target_chem, chebi2info)
        
        # instead of using chem_c --> chem_e system, make a reaction which maeks chem --> None
        candidate_reactions += _make_external_rxn(met_bigg_id, target_chem_object)

        _curate_reaction(candidate_reactions)
        
        model.add_reactions(candidate_reactions)
        target_chem_name = model.metabolites.get_by_id(met_bigg_id+'_c').name
        
        _curate_model(met_bigg_id, model)
        
        model_name = f'{model.id}_{met_bigg_id}'
        if model_name in added_reactions:
            model_name += f'_path_{target_count[target_chem]}'
                
        # added_reactions[met_bigg_id] = [added_rxn.id for added_rxn in candidate_reactions]
        # added_reactions_name[target_chem_name] = [added_rxn.name for added_rxn in candidate_reactions]
        # added_reactions_name[model.id+'_'+met_bigg_id] = [target_chem_name]+[added_rxn.name for added_rxn in candidate_reactions]

        # model_name = f'{model.id}_{met_bigg_id}'
        model.id = model_name

        target_models[model_name] = model
        added_reactions[model_name] = [added_rxn.id for added_rxn in candidate_reactions]
        added_reactions_name[model_name] = [added_rxn.name for added_rxn in candidate_reactions]
                    
    return target_models, added_reactions, added_reactions_name


def make_addition_table(add_dict):
    path_max_length = max([len(val) for val in add_dict.values()])
    path_max_length
    for _, val in add_dict.items():
        val += ['']*(path_max_length-len(val))
    df = pd.DataFrame.from_dict(add_dict, ).T
    df.index.name = 'Model_name'
    return df


def chebi_manual_curation(chebi2info):
    chebi2info['CHEBI:17968'] = {'name':['Butyrate (n-C4:0)'], 'met_id':['but']} # butyrate
    chebi2info['CHEBI:17272'] = {'name':['Propionate (n-C3:0)'], 'met_id':['ppa']} # propionate
    chebi2info['CHEBI:57776'] = {'name':['N-Acetyl-D-glucosamine 1-phosphate'], 'met_id':['acgam1p']} # acgam1p
    chebi2info['CHEBI:78811'] = {'name':['Heptanoyl-COA'], 'met_id':['hptcoa']} # heptanoyl-CoA
    chebi2info['CHEBI:43074'] = {'name':['Hydroxymethylglutaryl CoA C27H39N7O20P3S'], 'met_id':['hmgcoa']} # heptanoyl-CoA
    chebi2info['CHEBI:32362'] = {'name':['Heptanoate'], 'met_id':['hpta']} # heptanoate
    chebi2info['CHEBI:15377'] = {'name':['H2O H2O'], 'met_id':['h2o']}
    chebi2info['CHEBI:28885'] = {'name':['Butanol'], 'met_id':['btoh']}
    chebi2info['CHEBI:28938'] = {'name':['Ammonium'], 'met_id':['nh4']}
    chebi2info['CHEBI:57367'] = {'name':['Propenoyl-CoA'], 'met_id':['prpncoa']}
    chebi2info['CHEBI:16982'] = {'name':['R R 2 3 Butanediol C4H10O2'], 'met_id':['btd_RR']}
    chebi2info['CHEBI:73142'] = {'name':['4-Hydroxy-2-oxohexanoic acid'], 'met_id':['4hoxoh']}

    chebi2info['CHEBI:57316'] = {'name':['(S)-3-Hydroxybutanoyl-CoA'], 'met_id':['3hbcoa']}
    chebi2info['CHEBI:46911'] = {'name':['Ornithine'], 'met_id':['orn']}
    return chebi2info


def model_manual_curation(model):
    model_reaction_list = [rxn.id for rxn in model.reactions]
    curation_reverse = [
        'HSDy', 'G3PD2', 'SULR_1', 'ALCD2y', 'ACOAD3',
        'ACOAD4', 'ACOAD5', 'ACOAD6', 'ACOAD7', 'PYDXO'
    ]
    curation_forward = [
        'SHK3Dr', 'GTHOr', 'DHFR', 'TRSARr', '13PPDH2',
        'DURADx', 'ALDD19xr', '13PPDH2_1', 'SULR', 'PUNP1',
        'PUNP2', 'PUNP3', 'PUNP4', 'PUNP5', 'PUNP6', 
        'PUNP7', 'PUNP8', 'ALCD19_L', 'ALCD19y', 'SHK3Dr',
        'ACOAD1', 'ACOAD2', 'DGK1', 'MAN1PG', 'PDX5POi',
        'FACOAL140', 'DHORDi', 'AHSERL4'
    ]
    for curate_rxn in curation_reverse:
        if curate_rxn not in model_reaction_list:
            continue
        tmp_rxn = model.reactions.get_by_id(curate_rxn)
        if tmp_rxn.bounds != (-1000, 0):
            tmp_rxn.bounds = (-1000, 0)
            print('%s (%s) direction curated'%(curate_rxn, tmp_rxn.name))
    for curate_rxn in curation_forward:
        if curate_rxn not in model_reaction_list:
            continue
        tmp_rxn = model.reactions.get_by_id(curate_rxn)
        if tmp_rxn.bounds != (0, 1000):
            tmp_rxn.bounds = (0, 1000)
            print('%s (%s) direction curated'%(curate_rxn, tmp_rxn.name))

    return model



def model_path_construction(model, target_chem, BIGG_MET_DATA, STOICH_DATA, MET_DATA, RXN_DATA, TARGET_DATA):
    chebi2info = extract_bigg_info(BIGG_MET_DATA) # {CHEBI_ID: {'name': -, 'met_id': -}}
    # universal metabolite curation
    chebi2info = chebi_manual_curation(chebi2info)

    multiple_met_ids = []
    for item in chebi2info:
        met_list = chebi2info[item]['met_id']
        if len(met_list) > 1:
            multiple_met_ids.append(item)
    logging.info(f'Number of CHEBI IDs with multiple BiGG metabolite IDs: {len(multiple_met_ids)}')


    stoich_info = extract_stoich(STOICH_DATA) # {RHEA_ID: {met_id: coeff}}
    met_info_table = extract_metabolite_info(MET_DATA)
    met_dict = getMetabolites(TARGET_DATA, stoich_info, chebi2info, model, met_info_table)
        # {met_id: cobra.Metabolite}
    logging.info('Number of stoich info: %d\tmetabolite info: %d\t metabolites object: %d'\
                %(len(stoich_info), len(met_info_table), len(met_dict)))


    rxn_info_table = extract_reaction_info(RXN_DATA)
    path2targets, _ = getReactions(TARGET_DATA, stoich_info, model, met_dict, rxn_info_table)
        # {CHEBI_ID: [[cobra.Reaction, cobra.Reaction], [cobra.Reaction]]}

    chem_paths = path2targets[target_chem] # [[cobra.Reaction, cobra.Reaction], [cobra.Reaction]]
    target_models, added_rxns, added_rxn_names = save_chemical_model(
        target_chem, chem_paths, chebi2info, model, met_dict
    )

    return target_models, added_rxns, added_rxn_names