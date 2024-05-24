import os
import re
from math import isnan
import pandas as pd
import numpy as np
from cobra.io import read_sbml_model
from cobra.flux_analysis.loopless import add_loopless
from cobra.flux_analysis import flux_variability_analysis

from tqdm import tqdm


def read_target_models(model_dir, model_name_list):
    target_models = {}
    with tqdm(len(model_name_list), position=0, leave=True) as pbar:
        for i, each_model in enumerate(tqdm(model_name_list, position=0, leave=True)):
            try:
                model = read_sbml_model(model_dir+each_model)
                target_models[each_model] = model
            except:
                print('reading model failed')
                raise NotImplementedError
                
            pbar.set_description('Model number %d\tTarget model: %s'%(i, each_model))
    return target_models


def initial_setting(target_models, constraints):
    for each_model in target_models:
        model = target_models[each_model]
        for each_rxn in constraints:
            model.reactions.get_by_id(each_rxn).bounds = constraints[each_rxn]
    return


def _make_anaerobic_condition(model):
    model.reactions.EX_o2_e.lower_bound = 0
    
    model.reactions.ATPM.lower_bound = 0 # NGAM = 0 Jouthen et al. 2012
    
    new_GAM = 30.49 # GAM = 30.49 Nissen et al. 1997
    for met_id in ['atp_c', 'adp_c', 'h_c', 'pi_c', 'h2o_c',]:
        prev_coeff = model.reactions.r_4041.get_coefficient(met_id)
        model.reactions.r_4041.subtract_metabolites({model.metabolites.get_by_id(met_id): prev_coeff})
        if prev_coeff < 0:
            model.reactions.r_4041.add_metabolites({model.metabolites.get_by_id(met_id): -new_GAM})
        else:
            model.reactions.r_4041.add_metabolites({model.metabolites.get_by_id(met_id): new_GAM})

    model.reactions.r_4041.subtract_metabolites({
            model.metabolites.protein_c: 1
        })
    model.reactions.r_4041.add_metabolites({
        model.metabolites.protein_c: 0.461 # P = 0.461 Nissen et al. 1997
    })
    
    for met_id in ['hemeA_c', 'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'coa_c']:
        prev_coeff = model.reactions.r_4598.get_coefficient(met_id)
        model.reactions.r_4598.subtract_metabolites({model.metabolites.get_by_id(met_id): prev_coeff})
        
    for rxn_id in ['EX_ergst_e', 'EX_lanost_e', 'EX_hdcea_e', 'EX_zymst_e', 'r_2134', 'r_2137', 'EX_ocdcea_e']:
        model.reactions.get_by_id(rxn_id).lower_bound = -1000
        
    for rxn_id in ['MDHm', 'MDH']:
        model.reactions.get_by_id(rxn_id).lower_bound = 0
    for rxn_id in ['GLUSx']:
        model.reactions.get_by_id(rxn_id).upper_bound = 0
        
    return model



def _optimize_model(model, carbon_source, prev_lb, carbon_num, target_chem_id, achievable):
    with model as m:
        m.reactions.get_by_id(carbon_source).lower_bound = prev_lb * carbon_num / 6
        if achievable:
            max_bio = m.slim_optimize()
            if isnan(max_bio):
                return 0, 0
            elif abs(max_bio) < 1e-3:
                return 0, 0
            else:
                for rxn in m.reactions:
                    if rxn.objective_coefficient == 1:
                        biomass_rxn = rxn.id 
                m.reactions.get_by_id(biomass_rxn).lower_bound = max_bio*0.1
        else:
            m.reactions.ATPM.lower_bound = 0
        obj_reaction = 'EXT_%s_c'%target_chem_id
        m.objective = obj_reaction
        m.objective_direction = 'max'
        sol = m.optimize()

    if sol.status == 'infeasible':
        target_flux = 0
        carbon_flux = 0
    elif abs(sol.fluxes[carbon_source]) < 1e-3:
        target_flux = 0
        carbon_flux = 0
    else:
        target_flux = max(0, sol.fluxes[obj_reaction])
        carbon_flux = -sol.fluxes[carbon_source]
        
    return target_flux, carbon_flux



def _yield_calculate(target_models, model_name, carbon_source, carbon_num, prev_lb, achievable=False):
    results = {}
    multi_paths = []
    p_multi = re.compile('%s_(\w+)_path_(\d+).xml'%model_name)
    p = re.compile('%s_(\w+).xml'%model_name)
    with tqdm(len(target_models), position=0, leave=True) as pbar:
        for i, each_model in enumerate(tqdm(target_models, position=0, leave=True) ):
            multipath_check = p_multi.match(each_model)
            if multipath_check:
                multi_paths.append(multipath_check.group(1))
                continue
            model = target_models[each_model]
            target_chem_id = p.match(each_model).group(1)
            target_chem = model.metabolites.get_by_id(target_chem_id+'_c')
            target_chem_name = target_chem.name

            target_chem_MW = target_chem.formula_weight
            carbon_source_MW = model.metabolites.get_by_id(carbon_source[3:]).formula_weight

            target_flux, carbon_flux = _optimize_model(model, carbon_source, prev_lb, carbon_num, target_chem_id, achievable)
            
            if carbon_flux == 0:
                molar_yield = 0
            else:
                molar_yield = np.divide(target_flux, carbon_flux)
            molar_c_yield = molar_yield / carbon_num
            gram_yield = molar_yield * target_chem_MW / carbon_source_MW
            chem_yield = (target_chem_id, molar_yield, molar_c_yield, gram_yield)
            results[target_chem_name] = chem_yield

            pbar.set_description('Model number %s\tTarget chemical: %s'%(i, target_chem_id))
        
    yield_table = pd.DataFrame.from_dict(results).T
    yield_table.columns = ['target ID', 'molar yield (mol/mol)', 'molar yield (mol/C mol)', 'gram yield (g/g)']
    yield_table = yield_table.round(6)
    yield_table = yield_table.sort_index()
    
    multi_paths = list(set(multi_paths))
    multi_paths.sort()
    return yield_table, multi_paths


def _multipath_yield_calculate(multi_paths, target_models, model_name, carbon_source, carbon_num, prev_lb, achievable=False):
    results = {}
    p_multi = re.compile('%s_(\w+)_path_(\d+).xml'%model_name)
    p = re.compile('%s_(\w+).xml'%model_name)
    with tqdm(len(multi_paths), position=0, leave=True) as pbar:
        for i, each_target in enumerate(tqdm(multi_paths, position=0, leave=True)):
            results_target = {}
            for each_model in target_models:
                
                p_compile = p_multi.match(each_model)
                if not p_compile:
                    target_chem_id = p.match(each_model).group(1)
                    target_num = 0
                else:
                    target_chem_id = p_compile.group(1)
                    target_num = p_compile.group(2)
                
                if each_target != target_chem_id:
                    continue

                model = target_models[each_model]
                    
                target_chem = model.metabolites.get_by_id(target_chem_id+'_c')
                target_chem_name = target_chem.name
                target_chem_MW = target_chem.formula_weight
                carbon_source_MW = model.metabolites.get_by_id(carbon_source[3:]).formula_weight
                
                target_flux, carbon_flux = _optimize_model(model, carbon_source, prev_lb, carbon_num, target_chem_id, achievable)

                if carbon_flux == 0:
                    molar_yield = 0
                else:
                    molar_yield = np.divide(target_flux, carbon_flux)
                molar_c_yield = molar_yield / carbon_num
                gram_yield = molar_yield * target_chem_MW / carbon_source_MW
                
                chem_yield = (target_chem_id, molar_yield, molar_c_yield, gram_yield)
                results_target[f'{target_chem_name}_path_{target_num}'] = chem_yield
                results_target_table = pd.DataFrame.from_dict(results_target).T
                results_target_table.columns = ['target ID', 'molar yield (mol/mol)', 'molar yield (mol/C mol)', 'gram yield (g/g)']
                results_target_table = results_target_table.round(6)
                results_target_table.sort_index()
            
            results[each_target] = results_target_table
            pbar.set_description('Model number %s\tTarget chemical: %s'%(i, target_chem_id))
        
    return results


def saveYieldCalculationResult(target_models, model_name, carbon_source, carbon_num, prev_lb, ngam, achievable=False):
    yield_table, multi_paths = _yield_calculate(target_models, model_name, carbon_source, carbon_num, prev_lb, achievable)
    multipath_yield_tables = _multipath_yield_calculate(multi_paths, target_models, model_name, carbon_source, carbon_num, prev_lb, achievable)
    multi_max_yield = {}
    for item in multipath_yield_tables:
        tmp_dir = './outputs/%s/multi_paths/%s'%(model_name, item)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
            
        multipath_yield_tables[item].to_csv('%s/%s_%s.csv'%(tmp_dir, carbon_source[3:-2], ngam))
        df = multipath_yield_tables[item]
        mask = df['molar yield (mol/mol)'] != 0
        df = df[mask].sort_values('molar yield (mol/mol)', ascending=False)
        if len(df) == 0:
            continue
        multi_max_yield[item] = df.iloc[0]
    
    return yield_table, multi_max_yield



def calculate_yield(model, c_source='glc__D', prev_c_source='glc__D', air='aerobic', achievable=False, loopless=True):
    target_chem_id = '_'.join(model.id.split('_')[1:])
    target_chem_id = target_chem_id.split('_path_')[0]
    target_chem = model.metabolites.get_by_id(target_chem_id + '_c')
    target_chem_name = target_chem.name
    target_chem_MW = target_chem.formula_weight
    obj_reaction = f'EXT_{target_chem_id}_c'


    p_C_num = re.compile('C(\d*)')

    c_source_chem = model.metabolites.get_by_id(c_source+'_e')
    c_source_MW = c_source_chem.formula_weight
    c_source_C_num = float(p_C_num.match(c_source_chem.formula).group(1))
    c_source_rxn_id = f'EX_{c_source}_e'


    prev_c_source_chem = model.metabolites.get_by_id(prev_c_source+'_e')

    prev_c_source_C_num = float(p_C_num.match(prev_c_source_chem.formula).group(1))
    prev_c_source_rxn_id = f'EX_{prev_c_source}_e'

    carbon_ratio = c_source_C_num / prev_c_source_C_num

    for rxn in model.reactions:
        if rxn.objective_coefficient == 1:
            biomass_rxn = rxn.id 

    with model as m:
        if air == 'aerobic':
            m.reactions.EX_o2_e.lower_bound = -1000
        elif air == 'anaerobic':
            if 'yeast' in m.id:
                m = _make_anaerobic_condition(m)
            m.reactions.EX_o2_e.lower_bound = 0
        elif air == 'microaerobic':
            if 'yeast' in m.id:
                m = _make_anaerobic_condition(m)
            m.reactions.EX_o2_e.lower_bound = -0.5
        else:
            raise Exception("Wrong air condition")
        

        if loopless:
            add_loopless(m)

        
        prev_c_source_rxn = m.reactions.get_by_id(prev_c_source_rxn_id)
        prev_lb = prev_c_source_rxn.lower_bound
        prev_c_source_rxn.lower_bound = 0
        m.reactions.get_by_id(c_source_rxn_id).lower_bound = prev_lb / carbon_ratio
        if achievable:
            max_bio = m.slim_optimize()
            if isnan(max_bio):
                return 0, 0, 0
            elif abs(max_bio) < 1e-3:
                return 0, 0, 0
            else:
                m.reactions.get_by_id(biomass_rxn).lower_bound = max_bio * 0.1 
        else:
            m.reactions.ATPM.bounds = (0, 0)

        m.objective = obj_reaction
        m.objective_direction = 'max'
        sol = m.optimize()
        max_f = sol.fluxes[obj_reaction]
        min_c = sol.fluxes[c_source_rxn_id]

        # max_f = m.slim_optimize()
        if isnan(max_f):
            return 0, 0, 0
        elif abs(max_f) < 1e-3:
            return 0, 0, 0
        else:
            m.reactions.get_by_id(obj_reaction).bounds = (max_f, max_f)
            m.objective = c_source_rxn_id
            m.objective_direction = 'max'
            opt_c = m.slim_optimize()

            if not isnan(opt_c):
                min_c = opt_c
        
        if isnan(min_c):
            return 0, 0, 0
        elif abs(min_c) < 1e-3:
            return 0, 0, 0
        
        molar_yield = max_f / abs(min_c)
        gram_yield = (max_f * target_chem_MW) / (abs(min_c) * c_source_MW)
        molar_C_yield = molar_yield / c_source_C_num
        
        return molar_yield, gram_yield, molar_C_yield
    


def calculate_yield_fva(model, c_source='glc__D', prev_c_source='glc__D', air='aerobic', achievable=False, loopless=True):
    target_chem_id = '_'.join(model.id.split('_')[1:])
    target_chem_id = target_chem_id.split('_path_')[0]
    target_chem = model.metabolites.get_by_id(target_chem_id + '_c')
    target_chem_name = target_chem.name
    target_chem_MW = target_chem.formula_weight
    obj_reaction = f'EXT_{target_chem_id}_c'


    p_C_num = re.compile('C(\d*)')

    c_source_chem = model.metabolites.get_by_id(c_source+'_e')
    c_source_MW = c_source_chem.formula_weight
    c_source_C_num = float(p_C_num.match(c_source_chem.formula).group(1))
    c_source_rxn_id = f'EX_{c_source}_e'


    prev_c_source_chem = model.metabolites.get_by_id(prev_c_source+'_e')

    prev_c_source_C_num = float(p_C_num.match(prev_c_source_chem.formula).group(1))
    prev_c_source_rxn_id = f'EX_{prev_c_source}_e'

    carbon_ratio = c_source_C_num / prev_c_source_C_num

    for rxn in model.reactions:
        if rxn.objective_coefficient == 1:
            biomass_rxn = rxn.id 

    with model as m:
        if air == 'aerobic':
            m.reactions.EX_o2_e.lower_bound = -1000
        elif air == 'anaerobic':
            if 'yeast' in m.id:
                m = _make_anaerobic_condition(m)
            m.reactions.EX_o2_e.lower_bound = 0
        elif air == 'microaerobic':
            if 'yeast' in m.id:
                m = _make_anaerobic_condition(m)
            m.reactions.EX_o2_e.lower_bound = -0.5
        else:
            raise Exception("Wrong air condition")
        
        
        prev_c_source_rxn = m.reactions.get_by_id(prev_c_source_rxn_id)
        prev_lb = prev_c_source_rxn.lower_bound
        prev_c_source_rxn.lower_bound = 0
        m.reactions.get_by_id(c_source_rxn_id).lower_bound = prev_lb / carbon_ratio

        max_bio = m.slim_optimize()
        if achievable:
            m.reactions.get_by_id(biomass_rxn).lower_bound = max_bio * 0.1
        else:
            m.reactions.ATPM.knock_out()
        
        fva = flux_variability_analysis(
            m, reaction_list=[obj_reaction], 
            loopless=loopless, fraction_of_optimum=0.0,
            processes=1
        )

        max_f = fva['maximum'][0]

        if isnan(max_f):
            return 0, 0, 0
        elif abs(max_f) < 1e-3:
            return 0, 0, 0
        
        m.reactions.get_by_id(obj_reaction).bounds = (max_f, max_f)
        m.objective = c_source_rxn_id
        min_c = m.slim_optimize()
                
        if isnan(min_c):
            return 0, 0, 0
        elif abs(min_c) < 1e-3:
            return 0, 0, 0
        
        molar_yield = max_f / abs(min_c)
        gram_yield = (max_f * target_chem_MW) / (abs(min_c) * c_source_MW)
        molar_C_yield = molar_yield / c_source_C_num
        
        return molar_yield, gram_yield, molar_C_yield