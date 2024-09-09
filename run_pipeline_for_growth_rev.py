import re
import os
from glob import glob
from math import isnan
from copy import deepcopy

import logging

import yaml
import pandas as pd
from tqdm import tqdm
from cobra.io import read_sbml_model, write_sbml_model, load_json_model

from baebra.utils import argument_parser_pipeline
from baebra.generator import make_addition_table
from baebra.generator import model_manual_curation
from baebra.generator import model_path_construction
from baebra.calculator import calculate_yield_growth
from baebra.yieldfinder import set_model_condition
from baebra.yieldfinder import validate_hetero_target_growth
from baebra.yieldfinder import predict_heterologous_reactions_growth



logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    handlers=[logging.StreamHandler()]
                    )

logger = logging.getLogger('YieldPipeline')


if __name__ == '__main__':
    parser = argument_parser_pipeline()
    options = parser.parse_args()

    with open(options.base_config, 'r') as f:
        config = yaml.full_load(f)
    
    if options.config:
        with open(f'./config/{options.config}.yaml', 'r') as f:
            config_tmp = yaml.full_load(f)
        config.update(config_tmp)

    if options.input_model:
        config['model_path'] = options.input_model
    if options.target_chem:
        config['target_id'] = options.target_chem
    if options.output_dir:
        config['output_path'] = options.output_dir


    target_chem = config['target_id']
    output_dir = config['output_path']

    c_sources = config['carbon_sources']
    air_conditions = config['air_conditions']
    simulations_types = config['targeting_simulation']
    maintenance = config['maintenance']

    universal_model_dir = config['univ_path']
    num_cpu = config['cpu_number']
    limit_reaction_num = config['num_hetero']

    limit_reaction_num = 1


    if output_dir[-1] == '/':
        output_dir = output_dir[:-1]

    output_dir = output_dir + '/' + target_chem
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    target_models = {}
    for model_dir in glob(f'{output_dir}/*.xml'):
        model_name = model_dir.split('/')[-1].split('.xml')[0]
        each_model = read_sbml_model(model_dir)
        target_models[model_name] = each_model



    logger.info('Calculating maximum yields ... ')
    result_yield = []
    for model_name, each_model in tqdm(target_models.items()):
        for c_source in c_sources:
            for air in air_conditions:
                molar_yield, gram_yield, molar_C_yield = calculate_yield_growth(
                    model=each_model, c_source=c_source, air=air,
                )
                tmp = [
                    model_name, c_source, air, 
                    'YT', molar_yield, gram_yield, molar_C_yield
                ]
                result_yield.append(tmp)


    df_yield = pd.DataFrame(result_yield,)
    df_yield.columns = [
        'Model', 'Carbon_source', 'Air_condition', 'Yield_type', 
        'Molar_yield', 'Gram_yield', 'Molar_C_yield'
    ]
    df_yield.to_csv(output_dir + '/result_yield.txt', sep='\t')
    logger.info('Calculated maximum yields are saved')

    # heterologous reaction finder
    logger.info('Finding Heterologous reaction targets ...')
    manually_curated_rxns = [
        'MDHy', 'HADPCOADH', 'NAD_H2',
        'HYDFDN', 'HYDFDN2r', 'HYDFDi',
        'FAO1', 'FAO2', 'FAO3', 'FAO4', 'FAO10', 'FAO11', 
    ]

    universal_model = load_json_model(universal_model_dir)
    rm_reactions = []
    for each_reaction in universal_model.reactions:
        if each_reaction.check_mass_balance() != {}:
            rm_reactions.append(each_reaction)
            continue
        compartments = []
        for met in each_reaction.reactants + each_reaction.products:
            compartments.append(met.id[-1])
        compartments = list(set(compartments))
        if len(compartments) == 1 and compartments[0] == 'c':
            pass
        else:
            rm_reactions.append(each_reaction)

        if each_reaction.id in manually_curated_rxns:
            rm_reactions.append(each_reaction)
    

    logger.info('No. of removed reactions in universal model: %d'%(len(rm_reactions)))
    universal_model.remove_reactions(rm_reactions)


    ###

    for rxn in universal_model.reactions:
        rxn.bounds = (-1000, 1000)

    ###

    result_hetero = []
    for model_name, each_model in tqdm(target_models.items()):
        for c_source in c_sources:
            for air in air_conditions:

                df_tmp = df_yield[df_yield['Model']==model_name]
                df_tmp = df_tmp[df_tmp['Carbon_source']==c_source]
                df_tmp = df_tmp[df_tmp['Air_condition']==air]
                df_tmp = df_tmp[df_tmp['Yield_type']=='YT']
                yield_original = df_tmp['Molar_yield'].iloc[0]

                target_model = deepcopy(each_model)
                target_model = set_model_condition(
                    target_model, c_source, 
                    prev_c_source=c_source, air=air,
                    maintenance=None
                )
                target_identified = predict_heterologous_reactions_growth(
                    target_model, universal_model, f'EX_{c_source}_e',
                    limit_reaction_num, num_cpu, 
                )

                if target_identified == False:
                    print('NOOO')
                    continue
                
                valid_targets = validate_hetero_target_growth(
                    target_model, universal_model, 
                    target_identified, f'EX_{c_source}_e'
                )

                for key, val in valid_targets.items():
                    if val - yield_original > 1e-4: # 1e-3 --> 1e-4
                        tmp = [
                            model_name, c_source, air,
                            key, val, yield_original, (val/yield_original-1)*100
                        ]
                        result_hetero.append(tmp)


    if len(result_hetero) > 0:
        df_hetero = pd.DataFrame(result_hetero,)
        df_hetero.columns = [
            'Model', 'Carbon_source', 'Air_condition', 'Target_reactions',
            'Improved_yield', 'Previous_yield', 'Increased_percentage'
        ]
        df_hetero.to_csv(output_dir + '/result_hetero.txt', sep='\t')

        logger.info('Heterologous targets are saved')

    else:
        logger.info('No heterologous reactions were predicted')