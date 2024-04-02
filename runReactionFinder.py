import os
import re
import copy
import logging
import warnings

import numpy as np
import pandas as pd

from cobra.io import read_sbml_model, load_json_model
from cobra.io import write_sbml_model

from yieldFinder.utils import argument_parser
from yieldFinder.Simulator import Simulator
from yieldFinder.Gap import GapFilling
from baebra.calculator import _make_anaerobic_condition



def _recordRxnStoich(reaction):
    met_ids = [met.id for met in reaction.metabolites]
    met_ids.sort()
    met_coeffs = [reaction.get_coefficient(met_id) for met_id in met_ids]
    inverse_coeffs = [-reaction.get_coefficient(met_id) for met_id in met_ids]
    return met_ids, met_coeffs, inverse_coeffs


def _filter_same_reactions(model, universal_model):
    model_met_assembly = []
    model_stoich_assembly = []
    model_inv_stoich_assembly = []
    for rxn in model.reactions:
        met_id_info, met_coeff_info, met_inv_coeff = _recordRxnStoich(rxn)
        model_met_assembly.append(met_id_info)
        model_stoich_assembly.append(met_coeff_info)
        model_inv_stoich_assembly.append(met_inv_coeff)

    add_reactions = []
    target_model_rxns = [rxn.id for rxn in model.reactions]

    cfc_pair = {'nadh_c': 'nad_c', 'nadph_c': 'nadp_c'}
    new_cfc_pair = {'nadh_c':'nadph_c', 'nadph_c':'nadh_c'}

    for univ_rxn in universal_model.reactions:
        if univ_rxn.id in target_model_rxns:
            continue

        # delete insufficiently annotated data
        if univ_rxn.name == '':
            continue

        cnt = 0
        univ_met_id_info, univ_met_stoich_info, _ = _recordRxnStoich(univ_rxn)

        # check whether cofactor exchanged reaction exists in the prev model
        cfc_mets = set(['nadh_c', 'nadph_c']) & set(univ_met_id_info)
        if len(cfc_mets) == 1:
            current_cfc = list(cfc_mets)[0]
            current_pair = cfc_pair[current_cfc]
            cfc_exchanged_met_id_info = univ_met_id_info.copy()
            for i, met in enumerate(cfc_exchanged_met_id_info):
                if met == current_cfc:
                    cfc_exchanged_met_id_info[i] = new_cfc_pair[current_cfc]
                elif met == cfc_pair[current_cfc]:
                    cfc_exchanged_met_id_info[i] = cfc_pair[new_cfc_pair[current_cfc]]
            cfc_participate = True

        else:
            cfc_participate = False
        
   
        for i, model_met_id_info in enumerate(model_met_assembly):
            if univ_met_id_info == model_met_id_info:
                if univ_met_stoich_info == model_stoich_assembly[i]:
                    cnt += 1
                elif univ_met_stoich_info == model_inv_stoich_assembly[i]:
                    cnt += 1
            elif cfc_participate:
                if cfc_exchanged_met_id_info == model_met_id_info:
                    if cfc_exchanged_met_id_info == model_stoich_assembly[i]:
                        cnt += 1
                    elif cfc_exchanged_met_id_info == model_inv_stoich_assembly[i]:
                        cnt += 1

        if not cnt:
            add_reactions.append(univ_rxn)

    remove_reactions = [rxn for rxn in universal_model.reactions if rxn not in add_reactions]

    return add_reactions, remove_reactions


def predict_heterologous_reactions(target_model, universal_model, 
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
        return False
    prev_theoretical_yield = -c[target_reaction]/c[c_source]
    logger.info('Original max yield: %0.6f'%prev_theoretical_yield)
    if prev_theoretical_yield < 1e-3:
        logger.info('Chemical not produced')
        return False

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
            if each_reaction.id in target_cobra_model_reactions:
                rm_reactions.append(each_reaction)
            pass
        else:
            rm_reactions.append(each_reaction)

    logger.info('No. of removed reactions in universal model: %d'%(len(rm_reactions)))
    universal_model.remove_reactions(rm_reactions)

    added_reactions, remove_reactions = _filter_same_reactions(original_cobra_model, universal_model)
    logger.info('No. of added reactions: %d'%(len(added_reactions)))
    universal_model.remove_reactions(remove_reactions)

    original_cobra_model.add_reactions(added_reactions)


    obj = Simulator()
    obj.load_cobra_model(original_cobra_model)
    a,b,c = obj.run_FBA(new_objective=target_reaction)
    max_yield = c[target_reaction]
    if abs(c[c_source]) < 1e-3:
        logger.info('Unable to uptake carbon source')
        return False
    theoretical_yield = -c[target_reaction]/c[c_source]
    if theoretical_yield - prev_theoretical_yield < 1e-3:
        logger.info('no yield improvement')
        return False
    elif theoretical_yield < 1e-3:
        logger.info('zero yield')
        return False
    
    logger.info('Improved max yield: %0.6f'%theoretical_yield)

    start_yield = original_max_yield
    limit_yield = original_max_yield * 1.5
    logger.info('From %s to %s'%(start_yield, limit_yield))
    yield_list = np.linspace(start_yield, limit_yield, num_step)
    target_identified = False
    identified_targets = []

    tested_model = copy.deepcopy(target_model)
    obj = GapFilling()
    obj.set_limit_number(limit_reaction_num)
    obj.set_cpu_number(num_cpu)

    fp = open(output_file, 'w')
    fp.write('%s\t%s\t%s\t%s\n'%('Target reaction', 'Solution', 'Yield constraint', 'Candidate reactions'))
    for each_constraint in yield_list[1:]:
        logger.info('Yield constraint : %0.3fX\t%0.6f'%(each_constraint/original_max_yield, each_constraint))
        obj.set_threshold(each_constraint)
        obj.set_additional_constr(identified_targets)

        result_info = obj.run_gap_filling(universal_model, tested_model, target_reaction)
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
    return target_identified


def validate_results(target_model, universal_model, target_reaction, result_file, output_file, c_source='EX_glc__D_e'):
    df = pd.read_table(result_file)

    fp = open(output_file, 'w')
    fp.write('%s\t%s\t%s\n'%('Target', 'Candidate reactions', 'Maximum yield'))
    for each_group, each_df in df.groupby('Candidate reactions'):
        reactions = each_group.split(';')

        added_reactions = []
        for each_reaction in reactions:
            obj = universal_model.reactions.get_by_id(each_reaction)
            added_reactions.append(obj)

        temp_model = copy.deepcopy(target_model)
        temp_model.add_reactions(added_reactions)

        obj = Simulator()
        obj.load_cobra_model(temp_model)
        a,b,c = obj.run_FBA(new_objective=target_reaction)
        fp.write('%s\t%s\t%s\n'%(target_reaction, each_group, abs(c[target_reaction]/c[c_source])))
    fp.close()
    return


logger = logging.getLogger('HeteroFinder')
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
    universal_model_dir = options.universal_model
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

    # if aerobic:
    #     logger.info('Aerobic condition')
    # else:
    #     model.reactions.EX_o2_e.bounds = (0, 0)
    #     if 'yeast850' in model_name:
    #         model = _make_anaerobic_condition(model)
    #         logger.info('Yeast anaerobic configuration')
    #     logger.info('Anaerobic condition')

    # model.reactions.EX_o2_e.bounds = (-0.5, 0)
    if 'yeast850' in model_name:
        model = _make_anaerobic_condition(model)
        model.reactions.EX_o2_e.bounds = (-0.5, 0)
        model.reactions.GLYCDy.upper_bound = 1000
        logger.info('Yeast microaerobic configuration')
    # logger.info('Microaerobic condition')
    
    prev_c_lb = model.reactions.EX_glc__D_e.lower_bound
    model.reactions.EX_glc__D_e.lower_bound = 0
    
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
    

    if 'xml' in universal_model_dir:
        universal_model = read_sbml_model(universal_model_dir)
    elif 'json' in universal_model_dir:
        universal_model = load_json_model(universal_model_dir)
    logger.info('Read universal model: %s'%universal_model_dir)
    logger.info('Number of reactions in universal model: %d'%(len(universal_model.reactions)))

    warnings.filterwarnings(action='default')

    if multi_path:
        output_file = '%s/%s%s.txt'%(output_dir, target_chem, multi_path) 
        validation_output_file = '%s/%s%s_validated.txt'%(output_dir, target_chem, multi_path)
    else:
        output_file = '%s/%s.txt'%(output_dir, target_chem) 
        validation_output_file = '%s/%s_validated.txt'%(output_dir, target_chem)
    
    prediction_success = predict_heterologous_reactions(model, universal_model, target_reaction, limit_reaction_num, num_cpu, output_file, 'EX_%s_e'%carbon_source_type)
    if prediction_success:
        validate_results(model, universal_model, target_reaction, output_file, validation_output_file, 'EX_%s_e'%carbon_source_type)
