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
from baebra.yieldfinder import validate_hetero_target_fva
from baebra.yieldfinder import predict_heterologous_reactions


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    handlers=[logging.StreamHandler()]
                    )

logger = logging.getLogger('YieldPipeline')
# logger.setLevel(logging.INFO)


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

    universal_model_dir = config['univ_path']
    num_cpu = config['cpu_number']
    limit_reaction_num = config['num_hetero']


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


    df_yield = pd.read_csv(output_dir + '/result_yield.txt', sep='\t', index_col='Unnamed: 0')
    logger.info('Calculated yields are loaded')

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

    result_hetero = []
    for model_name, each_model in tqdm(target_models.items()):
        for c_source in c_sources:
            for air in air_conditions:

                df_tmp = df_yield[df_yield['Model']==model_name]
                df_tmp = df_tmp[df_tmp['Carbon_source']==c_source]
                df_tmp = df_tmp[df_tmp['Air_condition']==air]
                df_tmp = df_tmp[df_tmp['Yield_type']=='YT']
                yield_original = df_tmp['Molar_yield'].iloc[0]

                if yield_original == 0:
                    continue

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
                
                valid_targets = validate_hetero_target_fva(
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