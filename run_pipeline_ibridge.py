import os
from glob import glob
from copy import deepcopy
import logging

import yaml
import pandas as pd
from tqdm import tqdm
from cobra.io import read_sbml_model, load_json_model

from baebra.utils import argument_parser_pipeline
from baebra.yieldfinder import set_model_condition
from baebra.ibridge import predict_ibridge_targets


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
    maintenance = config['maintenance']

    num_cpu = config['cpu_number']

    exclude_met_dir = './inputs/excluded_metabolites.txt'


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
    
    manual_curation = ['EX_hdcea_e']
    for model_name, each_model in target_models.items():
        rxn_ids = [rxn.id for rxn in each_model.reactions]
        for rxn_id in manual_curation:
            if rxn_id in rxn_ids:
                each_model.reactions.get_by_id(rxn_id).bounds = (0, 0)
    
    logger.info('Constructed models are loaded')


    result_ibridge = []
    for model_name, each_model in tqdm(target_models.items()):
        for c_source in c_sources:
            for air in air_conditions:

                prev_c_source = 'glc__D'

                target_model = deepcopy(each_model)
                target_model = set_model_condition(
                    target_model, c_source, 
                    prev_c_source=prev_c_source, air=air,
                    maintenance=None
                )

                target_identified = predict_ibridge_targets(
                    target_model, exclude_met_dir, num_step=10,
                )

                if len(target_identified) == 0:
                    continue
                
                target_identified['Model'] = [model_name] * len(target_identified)
                target_identified['Carbon_source'] = [c_source] * len(target_identified)
                target_identified['Air_condition'] = [air] * len(target_identified)

                result_ibridge.append(target_identified)


    if len(result_ibridge) > 0:
        df_ibridge = pd.concat(result_ibridge, sort=False).reset_index(drop=True)
        df_ibridge.to_csv(output_dir + '/result_ibridge.txt', sep='\t')
        df_ibridge = df_ibridge[[
            'Model', 'Carbon_source', 'Air_condition', 'Target_reactions',
            'Regulation_type', 'Positive_metabolite', 'Negative_metabolite', 
            'Normalized_SoC', 'SoC',
        ]]
        df_ibridge.to_csv(output_dir + '/result_ibridge.txt', sep='\t')

        logger.info('iBridge targets are saved')

    else:
        logger.info('No iBridge target reactions were predicted')