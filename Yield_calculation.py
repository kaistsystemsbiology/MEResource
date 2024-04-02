import os
import re

import pandas as pd
import numpy as np
from cobra.io import read_sbml_model

from baebra.calculator import *
from baebra.utils import argument_parser
from baebra.calculator import _make_anaerobic_condition

if __name__ == '__main__':
    parser = argument_parser()
    options = parser.parse_args()
    model_name = options.input_model


    cobra_model = './inputs/%s.xml'%model_name
    model_dir = './TargetChemicalModels/%s/'%model_name

    ATPM_reaction = 'ATPM'
    prev_carbon_source = 'EX_glc__D_e'
    oxygen_uptake_reaction = 'EX_o2_e'


    base_model = read_sbml_model('./inputs/%s.xml'%model_name)
    base_reactions = [rxn.id for rxn in base_model.reactions]
    prev_lb = base_model.reactions.get_by_id(prev_carbon_source).lower_bound
    prev_ngam = base_model.reactions.get_by_id(ATPM_reaction).lower_bound


    model_name_list = [item for item in os.listdir(model_dir) if item[-3:]=='xml']
    model_name_list = [item for item in model_name_list if model_name in item]
    print('Number of model for the analysis: %s'%(len(model_name_list)))

    target_models = read_target_models(model_dir, model_name_list)


    carbon_source_list = ['EX_glc__D_e', 'EX_fru_e', 'EX_gal_e', 'EX_glyc_e', 'EX_sucr_e', 'EX_xyl__D_e', 'EX_arab__L_e']
    carbon_numbers = [6, 6, 6, 3, 12, 5, 5]


    # Maximum achievable yield, aerobic
    initial_setting(target_models, {prev_carbon_source:(0,0)})
    ngam = 'achievable_aerobic'

    c_yield_table_ngam = {}
    for c_source, c_num in zip(carbon_source_list, carbon_numbers):
        if c_source not in base_reactions:
            print('%s not in base model'%c_source)
            continue
        tmp_table, multi_max_yield = saveYieldCalculationResult(target_models, model_name, c_source, c_num, prev_lb, ngam, achievable=True)
        for multi_chem in multi_max_yield:
            for i, data in enumerate(tmp_table['target ID']):
                if data == multi_chem:
                    tmp_table.iloc[i] = multi_max_yield[multi_chem]
        
        tmp_table.to_csv('./outputs/%s/%s_%s_yield_resource_%s.csv'%(model_name, model_name, c_source[3:-2], ngam))
        tmp_table.sort_values('molar yield (mol/mol)')
        c_yield_table_ngam[c_source] = tmp_table

        

    # Maximum theoretical yield, aerobic
    initial_setting(target_models, {prev_carbon_source:(0,0)})
    ngam = 'theoretical_aerobic'

    c_yield_table = {}
    for c_source, c_num in zip(carbon_source_list, carbon_numbers):
        if c_source not in base_reactions:
            print('%s not in base model'%c_source)
            continue
        tmp_table, multi_max_yield = saveYieldCalculationResult(target_models, model_name, c_source, c_num, prev_lb, ngam)
        for multi_chem in multi_max_yield:
            for i, data in enumerate(tmp_table['target ID']):
                if data == multi_chem:
                    tmp_table.iloc[i] = multi_max_yield[multi_chem]
                    
        tmp_table.to_csv('./outputs/%s/%s_%s_yield_resource_%s.csv'%(model_name, model_name, c_source[3:-2], ngam))
        tmp_table.sort_values('molar yield (mol/mol)')
        c_yield_table[c_source] = tmp_table


    if model_name == 'yeast850':
        for model_id, model in target_models.items():
            model = _make_anaerobic_condition(model)

    # Maximum achievable yield, anaerobic
    initial_setting(target_models, {prev_carbon_source:(0,0), oxygen_uptake_reaction:(0,0)})
    ngam = 'achievable_anaerobic'

    c_yield_table_ngam = {}
    for c_source, c_num in zip(carbon_source_list, carbon_numbers):
        if c_source not in base_reactions:
            print('%s not in base model'%c_source)
            continue
        tmp_table, multi_max_yield = saveYieldCalculationResult(target_models, model_name, c_source, c_num, prev_lb, ngam, achievable=True)
        for multi_chem in multi_max_yield:
            for i, data in enumerate(tmp_table['target ID']):
                if data == multi_chem:
                    tmp_table.iloc[i] = multi_max_yield[multi_chem]
        
        tmp_table.to_csv('./outputs/%s/%s_%s_yield_resource_%s.csv'%(model_name, model_name, c_source[3:-2], ngam))
        tmp_table.sort_values('molar yield (mol/mol)')
        c_yield_table_ngam[c_source] = tmp_table



    # Maximum theoretical yield, anaerobic
    initial_setting(target_models, {prev_carbon_source:(0,0), oxygen_uptake_reaction:(0,0)})
    ngam = 'theoretical_anaerobic'

    c_yield_table = {}
    for c_source, c_num in zip(carbon_source_list, carbon_numbers):
        if c_source not in base_reactions:
            print('%s not in base model'%c_source)
            continue
        tmp_table, multi_max_yield = saveYieldCalculationResult(target_models, model_name, c_source, c_num, prev_lb, ngam)
        for multi_chem in multi_max_yield:
            for i, data in enumerate(tmp_table['target ID']):
                if data == multi_chem:
                    tmp_table.iloc[i] = multi_max_yield[multi_chem]
                    
        tmp_table.to_csv('./outputs/%s/%s_%s_yield_resource_%s.csv'%(model_name, model_name, c_source[3:-2], ngam))
        tmp_table.sort_values('molar yield (mol/mol)')
        c_yield_table[c_source] = tmp_table

    if model_name == 'yeast850':
        for model_id, model in target_models.items():
            model.reactions.get_by_id('GLYCDy').upper_bound = 1000

    # Maximum achievable yield, microaerobic
    initial_setting(target_models, {prev_carbon_source:(0,0), oxygen_uptake_reaction:(-0.5,0)})
    ngam = 'achievable_microaerobic'

    c_yield_table_ngam = {}
    for c_source, c_num in zip(carbon_source_list, carbon_numbers):
        if c_source not in base_reactions:
            print('%s not in base model'%c_source)
            continue
        tmp_table, multi_max_yield = saveYieldCalculationResult(target_models, model_name, c_source, c_num, prev_lb, ngam, achievable=True)
        for multi_chem in multi_max_yield:
            for i, data in enumerate(tmp_table['target ID']):
                if data == multi_chem:
                    tmp_table.iloc[i] = multi_max_yield[multi_chem]
        
        tmp_table.to_csv('./outputs/%s/%s_%s_yield_resource_%s.csv'%(model_name, model_name, c_source[3:-2], ngam))
        tmp_table.sort_values('molar yield (mol/mol)')
        c_yield_table_ngam[c_source] = tmp_table



    # Maximum theoretical yield, microaerobic
    initial_setting(target_models, {prev_carbon_source:(0,0), oxygen_uptake_reaction:(-0.5,0)})
    ngam = 'theoretical_microaerobic'

    c_yield_table = {}
    for c_source, c_num in zip(carbon_source_list, carbon_numbers):
        if c_source not in base_reactions:
            print('%s not in base model'%c_source)
            continue
        tmp_table, multi_max_yield = saveYieldCalculationResult(target_models, model_name, c_source, c_num, prev_lb, ngam)
        for multi_chem in multi_max_yield:
            for i, data in enumerate(tmp_table['target ID']):
                if data == multi_chem:
                    tmp_table.iloc[i] = multi_max_yield[multi_chem]
                    
        tmp_table.to_csv('./outputs/%s/%s_%s_yield_resource_%s.csv'%(model_name, model_name, c_source[3:-2], ngam))
        tmp_table.sort_values('molar yield (mol/mol)')
        c_yield_table[c_source] = tmp_table