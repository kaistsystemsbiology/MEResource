import re
import os
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
from baebra.calculator import calculate_yield
from baebra.yieldfinder import validate_target
from baebra.yieldfinder import set_model_condition
from baebra.yieldfinder import predict_heterologous_reactions


logger = logging.getLogger('YieldPipeline')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s-%(name)s-%(levelname)s-%(message)s')


if __name__ == '__main__':
    parser = argument_parser_pipeline()
    options = parser.parse_args()

    with open(options.base_config, 'r') as f:
        config = yaml.full_load(f)
    
    if options.config:
        with open(options.config, 'r') as f:
            config_tmp = yaml.full_load(f)
        config.update(config_tmp)

    if options.input_model:
        config['model_path'] = options.input_model
    if options.target_chem:
        config['target_id'] = options.target_chem
    if options.output_dir:
        config['output_path'] = options.output_dir

    BIGG_MET_DATA = config['BIGG_MET_DATA']
    STOICH_DATA = config['STOICH_DATA']
    MET_DATA = config['MET_DATA']
    RXN_DATA = config['RXN_DATA']
    TARGET_DATA = config['TARGET_DATA']

    model_dir = config['model_path']
    target_chem = config['target_id']
    output_dir = config['output_path']

    c_sources = config['carbon_sources']
    air_conditions = config['air_conditions']
    simulations_types = config['targeting_simulation']
    maintenance = config['maintenance']

    universal_model_dir = config['univ_path']
    num_cpu = config['cpu_number']
    limit_reaction_num = config['num_hetero']


    if output_dir[-1] == '/':
        output_dir = output_dir[:-1]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)

    model = read_sbml_model(model_dir)
    # manual curation of target model
    model = model_manual_curation(model)

    logger.info('Constructing models for paths ... ')
    target_models, added_rxns, added_rxn_names = model_path_construction(
        model, target_chem, BIGG_MET_DATA, STOICH_DATA, MET_DATA, RXN_DATA, TARGET_DATA
    )

    df_added_rxn = make_addition_table(added_rxns)
    df_added_rxn_names = make_addition_table(added_rxn_names)
    
    df_added_rxn.to_csv(output_dir + '/added_reaction.txt', sep='\t')
    df_added_rxn_names.to_csv(output_dir + '/added_reaction_name.txt', sep='\t')
    for model_name, each_model in target_models.items():
        write_sbml_model(each_model, f'{output_dir}/{model_name}.xml')
    logger.info('Constructed models are saved')

    logger.info('Calculating maximum yields ... ')
    result_yield = []
    for model_name, each_model in tqdm(target_models.items()):
        for c_source in c_sources:
            for air in air_conditions:
                molar_yield, gram_yield, molar_C_yield = calculate_yield(
                    model=each_model, c_source=c_source, air=air, achievable=False
                )
                tmp = [
                    model_name, c_source, air, 
                    'YT', molar_yield, gram_yield, molar_C_yield
                ]
                result_yield.append(tmp)

                molar_yield, gram_yield, molar_C_yield = calculate_yield(
                    model=each_model, c_source=c_source, air=air, achievable=True
                )
                tmp = [
                    model_name, c_source, air, 
                    'YA', molar_yield, gram_yield, molar_C_yield
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
    universal_model = load_json_model(universal_model_dir)

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
                    prev_c_source='glc__D', air=air,
                    maintenance=maintenance
                )
                target_identified = predict_heterologous_reactions(
                    target_model, universal_model, f'EX_{c_source}_e',
                    limit_reaction_num, num_cpu, 
                )

                if target_identified == False:
                    continue
                
                valid_targets = validate_target(
                    target_model, universal_model, 
                    target_identified, f'EX_{c_source}_e'
                )

                for key, val in valid_targets.items():
                    if val - yield_original > 1e-3:
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




            