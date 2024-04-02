import re
import os
from copy import deepcopy
import logging

import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction, Metabolite

from baebra.generator import *
from baebra.utils import argument_parser


bigg_met_data = './inputs/bigg_models_metabolites.txt'
stoich_data = './SummarizedData/ChemMap_Table - Stoichiometry.tsv'
met_data= './SummarizedData/ChemMap_Table - Metabolite_Info.tsv'
rxn_data = './SummarizedData/ChemMap_Table - Reaction_Info.tsv'
target_data = './SummarizedData/ChemMap_Table - TargetPath.tsv'

logger = logging.getLogger('HeteroFinder')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s-%(name)s-%(levelname)s-%(message)s')


if __name__ == '__main__':
    parser = argument_parser()
    options = parser.parse_args()
    model_name = options.input_model

    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)

    cobra_model = './inputs/%s.xml'%model_name
    model = read_sbml_model(cobra_model)
    # manual curation of target model
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
            logger.info('%s (%s) direction curated'%(curate_rxn, tmp_rxn.name))
    for curate_rxn in curation_forward:
        if curate_rxn not in model_reaction_list:
            continue
        tmp_rxn = model.reactions.get_by_id(curate_rxn)
        if tmp_rxn.bounds != (0, 1000):
            tmp_rxn.bounds = (0, 1000)
            logger.info('%s (%s) direction curated'%(curate_rxn, tmp_rxn.name))

    chebi2info = extract_bigg_info(bigg_met_data)

    multiple_met_ids = []
    for item in chebi2info:
        met_list = chebi2info[item]['met_id']
        if len(met_list) > 1:
            multiple_met_ids.append(item)
    logger.info(f'Number of CHEBI IDs with multiple BiGG metabolite IDs: {len(multiple_met_ids)}')


    # universal metabolite curation
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

    stoich_info = extract_stoich(stoich_data)
    met_info_table = extract_metabolite_info(met_data)
    met_dict = getMetabolites(target_data, stoich_info, chebi2info, model, met_info_table)
    logger.info('Number of stoich info: %d\tmetabolite info: %d\t metabolites object: %d'\
                %(len(stoich_info), len(met_info_table), len(met_dict)))


    rxn_info_table = extract_reaction_info(rxn_data)
    path2targets, wrong_dir = getReactions(target_data, stoich_info, model, met_dict, rxn_info_table)
    # manually_deleted_targets = ['CHEBI:350546', 'CHEBI:16796', 'CHEBI:5MTOTPAM']
    # for delete_target in manually_deleted_targets:
    #     if delete_target in path2targets:
    #         del path2targets[delete_target]
    logger.info(f'Number of target chemicals: {len(path2targets)}')

    if not os.path.exists('./TargetChemicalModels/%s'%model_name):
        os.makedirs('./TargetChemicalModels/%s'%model_name)
    added_table, added_table_name = save_chemical_models(path2targets, chebi2info, model, met_dict,
                        met_info_table, './TargetChemicalModels/%s/'%model_name, model_name=model_name)


    added_table = added_table.sort_index()
    added_table_name = added_table_name.sort_index()
    if not os.path.exists('./outputs/%s'%model_name):
        os.makedirs('./outputs/%s'%model_name)
    added_table.to_csv('./outputs/%s/%s_added_reaction_list.csv'%(model_name, model_name))
    added_table_name.to_csv('./outputs/%s/%s_added_reaction_name_list.csv'%(model_name, model_name))


    # if not os.path.exists('./outputs/models/%s'%model_name):
    #     os.makedirs('./outputs/models/%s'%model_name)
    # added_table, added_table_name = save_chemical_models(path2targets, chebi2info, model, met_dict,
    #                     met_info_table, './outputs/models/%s/'%model_name, model_name=model_name)

    # added_table = added_table.sort_index()
    # added_table_name = added_table_name.sort_index()
    # if not os.path.exists('./outputs/additional/%s'%model_name):
    #     os.makedirs('./outputs/additional/%s'%model_name)
    # added_table.to_csv('./outputs/additional/%s/%s_added_reaction_list.csv'%(model_name, model_name))
    # added_table_name.to_csv('./outputs/additional/%s/%s_added_reaction_name_list.csv'%(model_name, model_name))