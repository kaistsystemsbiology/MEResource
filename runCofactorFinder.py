import os
import re
import copy
import logging
import warnings

import numpy as np
import pandas as pd

from cobra.io import read_sbml_model, load_json_model
from cobra.io import write_sbml_model
from cobra import Reaction

from yieldFinder.utils import argument_parser
from yieldFinder.Simulator import Simulator
from yieldFinder.Cofactor import SwapCofactor



def _recordRxnStoich(reaction):
    met_ids = [met.id for met in reaction.metabolites]
    met_ids.sort()
    met_coeffs = [reaction.get_coefficient(met_id) for met_id in met_ids]
    inverse_coeffs = [-reaction.get_coefficient(met_id) for met_id in met_ids]
    return met_ids, met_coeffs, inverse_coeffs


def _swap_cofactor(model, rxn, prev_cofactor, new_cofactor):
    cfc_pair = {
        'nadh_c': 'nad_c',
        'nadph_c': 'nadp_c',
        # 'mql8_c': 'mqn8_c',
        # 'q8h2_c': 'q8_c'
    }
    rxn_stoich = {met: rxn.get_coefficient(met.id) for met in rxn.metabolites}
    new_rxn = Reaction('%s_swap_%s'%(rxn.id, new_cofactor))
    new_rxn.name = '%s_%s'%(rxn.name, new_cofactor)
    new_rxn.bounds = rxn.bounds
    new_rxn.add_metabolites(rxn_stoich)
    prev_coeff = new_rxn.get_coefficient(prev_cofactor)
    new_rxn.add_metabolites({
        model.metabolites.get_by_id(prev_cofactor): -prev_coeff,
        model.metabolites.get_by_id(new_cofactor): prev_coeff,
        model.metabolites.get_by_id(cfc_pair[prev_cofactor]): prev_coeff,
        model.metabolites.get_by_id(cfc_pair[new_cofactor]): -prev_coeff
    })
    return new_rxn


def _create_swap_reactions(model, rxn, cofactor_list):
    for met_id in cofactor_list:
        if met_id in [met.id for met in rxn.metabolites]:
            current_cfc = met_id
    new_reactions = []
    for met_id in cofactor_list:
        if met_id == current_cfc:
            continue
        new_rxn = _swap_cofactor(model, rxn, current_cfc, met_id)
        new_reactions.append(new_rxn)
    model.add_reactions(new_reactions)
    return new_reactions
    


def predict_cofactor_reactions(target_model, 
                                target_reaction, limit_reaction_num, num_cpu,
                                output_file, c_source='EX_glc__D_e'):
    num_step = 51
    original_cobra_model = copy.deepcopy(target_model)
    target_cobra_model_reactions = [reaction.id for reaction in target_model.reactions]

    obj = Simulator()
    obj.load_cobra_model(original_cobra_model)
    a,b,c = obj.run_FBA(new_objective=target_reaction)
    original_max_yield = c[target_reaction]
    
    if abs(c[c_source]) < 1e-3:
        logger.info('Chemical not produced')
        return False, None
    prev_theoretical_yield = -c[target_reaction]/c[c_source]
    logger.info('Original max yield: %0.6f'%prev_theoretical_yield)
    if original_max_yield < 1e-3:
        logger.info('Chemical not produced')
        return False, None

    # cofactor_list = ['nadh_c','nadph_c','mql8_c','q8h2_c']
    cofactor_list = ['nadh_c','nadph_c']
    met_lists = [original_cobra_model.metabolites.get_by_id(met) for met in cofactor_list]
    cofactor_use_rxns = []
    for rxn in original_cobra_model.reactions:
        if len(set(met_lists) & set(rxn.metabolites)) == 1:
            if rxn.gene_reaction_rule:
                cofactor_use_rxns.append(rxn)
    logger.info('Number of cofactor using reactions: %d'%(len(cofactor_use_rxns)))

    candidate_reactions = {}
    added_reactions = []
    for rxn in cofactor_use_rxns:
        gen_reactions = _create_swap_reactions(original_cobra_model, rxn, cofactor_list)
        candidate_reactions[rxn.id] = [each_rxn.id for each_rxn in gen_reactions]
        added_reactions += gen_reactions

    obj = Simulator()
    obj.load_cobra_model(original_cobra_model)
    a,b,c = obj.run_FBA(new_objective=target_reaction)
    max_yield = c[target_reaction]
    if abs(c[c_source]) < 1e-3:
        logger.info('Unable to uptake carbon source')
        return False, None
    
    theoretical_yield = -c[target_reaction]/c[c_source]
    if theoretical_yield - prev_theoretical_yield < 1e-3:
        logger.info('no yield improvement')
        return False, None
    elif theoretical_yield < 1e-3:
        logger.info('zero yield')
        return False, None
    
    logger.info('Improved max yield: %0.6f'%theoretical_yield)

    start_yield = original_max_yield
    limit_yield = original_max_yield * 1.5
    logger.info('From %s to %s'%(start_yield, limit_yield))
    yield_list = np.linspace(start_yield, limit_yield, num_step)
    target_identified = False
    identified_targets = []
    
    obj = SwapCofactor()
    obj.load_cobra_model(original_cobra_model)
    obj.set_limit_number(limit_reaction_num)
    obj.set_cpu_number(num_cpu)

    fp = open(output_file, 'w')
    fp.write('%s\t%s\t%s\t%s\n'%('Target reaction', 'Solution', 'Yield constraint', 'Candidate reactions'))
    for each_constraint in yield_list[1:]:
        logger.info('Yield constraint : %0.3fX\t%0.6f'%(each_constraint/original_max_yield, each_constraint))

        obj.set_threshold(each_constraint)
        obj.set_additional_constr(identified_targets)
        result_info = obj.run_swap_cofactor(original_cobra_model, target_reaction, candidate_reactions)
        if len(result_info) == 0:
            break
        for each_solution in result_info:
            target_reactions = result_info[each_solution]
            if len(target_reactions) > 0:
                target_identified = True
            target_reactions.sort()
            identified_targets.append(target_reactions)
            fp.write('%s\t%s\t%s\t%s\n'%(target_reaction, each_solution, each_constraint, ';'.join(target_reactions)))
    fp.close()
    return target_identified, original_cobra_model


def validate_results(target_model, target_reaction, result_file, output_file, c_source='EX_glc__D_e'):
    df = pd.read_table(result_file)

    fp = open(output_file, 'w')
    fp.write('%s\t%s\t%s\n'%('Target', 'Candidate reactions', 'Maximum yield'))
    for each_group, each_df in df.groupby('Candidate reactions'):
        reactions = each_group.split(';')

        added_reactions = []
        for each_reaction in reactions:
            obj = target_model.reactions.get_by_id(each_reaction)
            added_reactions.append(obj.id)

        added_rxn_bounds = {}
        for rxn_id in added_reactions:
            added_rxn_bounds[rxn_id] = target_model.reactions.get_by_id(rxn_id).bounds

        with target_model as temp_model:
            for rxn in temp_model.reactions:
                if '_swap' in rxn.id:
                    rxn.bounds = (0,0)
            
            obj = Simulator()
            obj.load_cobra_model(temp_model)
            a,b,c = obj.run_FBA(new_objective=target_reaction)
            prev_yield = abs(c[target_reaction]/c[c_source])
            for rxn_id in added_reactions:
                temp_model.reactions.get_by_id(rxn_id).bounds = added_rxn_bounds[rxn_id]

            obj = Simulator()
            obj.load_cobra_model(temp_model)
            a,b,c = obj.run_FBA(new_objective=target_reaction)
            improved_yield = abs(c[target_reaction]/c[c_source])
            if improved_yield < prev_yield:
                continue
        fp.write('%s\t%s\t%s\n'%(target_reaction, each_group, abs(c[target_reaction]/c[c_source])))
    fp.close()
    return

logger = logging.getLogger('CofactorFinder')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s-%(name)s-%(levelname)s-%(message)s')


if __name__ == '__main__':
    parser = argument_parser()
    options = parser.parse_args()

    output_dir = options.output_dir
    model_dir = options.input_model
    model_name = options.model_name
    log_dir = options.log_dir
    limit_reaction_num = options.num_reaction
    atpm = options.maintenance
    num_cpu = options.cpu_number
    aerobic = options.oxygen_availability
    carbon_source_type = options.carbon_source

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    p_target = re.compile('%s_(.*?)(_path_\d+)?[.]xml'%model_name)
    target_chem = p_target.search(model_dir).group(1)
    multi_path = p_target.search(model_dir).group(2)
    # target_reaction = 'EX_%s_e'%target_chem
    target_reaction = 'EXT_%s_c'%target_chem

    stream_handler = logging.StreamHandler()

    if multi_path:
        file_handler = logging.FileHandler('%s/%s%s_%s'%(output_dir, target_chem, multi_path, log_dir))
    else:
        file_handler = logging.FileHandler('%s/%s_%s'%(output_dir, target_chem, log_dir))
    
    file_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info('Use %d cpu'%num_cpu)
    warnings.filterwarnings(action='ignore')

    model = read_sbml_model(model_dir)
    logger.info('Read target model: %s'%model_dir)
    logger.info('Target reaction: %s'%target_reaction)
    logger.info('Target chemical: %s (%s)'%(model.metabolites.get_by_id(target_chem+'_c').name, target_chem))

    model.reactions.get_by_id(atpm).bounds = (0, 0)
    logger.info('Non growth associated ATP requirement knocked out: %s'%atpm)

    if aerobic:
        logger.info('Aerobic condition')
    else:
        model.reactions.EX_o2_e.bounds = (0, 0)
        logger.info('Anaerobic condition')
    # model.reactions.EX_o2_e.bounds = (-0.5, 0)
    # logger.info('Microaerobic condition')
    
    prev_c_lb = model.reactions.EX_glc__D_e.lower_bound
    model.reactions.EX_glc__D_e.lower_bound = 0
    
    model_reactions = [rxn.id for rxn in model.reactions]
    if f'EX_{carbon_source_type}_e' not in model_reactions:
        logger.info(f'No carbon source {carbon_source_type}')
        raise

    if carbon_source_type == 'glc__D':
        logger.info('Glucose as the carbon source')
        model.reactions.EX_glc__D_e.lower_bound = prev_c_lb
    elif carbon_source_type == 'fru':
        logger.info('Fructose as the carbon source')
        model.reactions.EX_fru_e.lower_bound = prev_c_lb
    elif carbon_source_type == 'gal':
        logger.info('Galactose as the carbon source')
        model.reactions.EX_gal_e.lower_bound = prev_c_lb
    elif carbon_source_type == 'glyc':
        logger.info('Glycerol as the carbon source')
        model.reactions.EX_glyc_e.lower_bound = prev_c_lb * 2
    elif carbon_source_type == 'sucr':
        logger.info('Sucrose as the carbon source')
        model.reactions.EX_sucr_e.lower_bound = prev_c_lb / 2
    elif carbon_source_type == 'xyl__D':
        logger.info('Xylose as the carbon source')
        model.reactions.EX_xyl__D_e.lower_bound = prev_c_lb*6/5
    elif carbon_source_type == 'arab__L':
        logger.info('Arabinose as the carbon source')
        model.reactions.EX_arab__L_e.lower_bound = prev_c_lb*6/5
    else:
        logger.warnings('Wrong carbon source')
        raise
    

    warnings.filterwarnings(action='default')

    if multi_path:
        output_file = '%s/%s%s.txt'%(output_dir, target_chem, multi_path) 
        validation_output_file = '%s/%s%s_validated.txt'%(output_dir, target_chem, multi_path)
    else:
        output_file = '%s/%s.txt'%(output_dir, target_chem) 
        validation_output_file = '%s/%s_validated.txt'%(output_dir, target_chem)
    
    prediction_success, analysis_model = predict_cofactor_reactions(model, target_reaction, limit_reaction_num, num_cpu, output_file, 'EX_%s_e'%carbon_source_type)
    if prediction_success:
        logger.info('prediction success')
        validate_results(analysis_model, target_reaction, output_file, validation_output_file, 'EX_%s_e'%carbon_source_type)
