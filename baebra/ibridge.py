import os
from math import isnan

import numpy as np
import pandas as pd

from cobra.flux_analysis import pfba, moma

os.environ['NUMEXPR_NUM_THREADS'] = '4'

def _extract_met_direction(model):
    dict_reactant = {}
    dict_product = {}
    for met in model.metabolites:
        dict_reactant[met.id] = []
        dict_product[met.id] = []
        for rxn in met.reactions:
            coeff = rxn.get_coefficient(met.id)
            if coeff < 0:
                dict_reactant[met.id].append(rxn.id)
            elif coeff > 0:
                dict_product[met.id].append(rxn.id)
            else:
                raise 

    return dict_reactant, dict_product


def match_met_direction(cobra_model, universal_model):
    dict_reactant_model, dict_product_model = _extract_met_direction(cobra_model)
    dict_reactant_univ, dict_product_univ = _extract_met_direction(universal_model)

    dict_reactant = dict_reactant_model.copy()
    dict_product = dict_product_model.copy()

    for key, val in dict_reactant_univ.items():
        if key not in dict_reactant:
            dict_reactant[key] = val
        else:
            prev_lst = dict_reactant[key]
            dict_reactant[key] = list(set(prev_lst) | set(val))

    for key, val in dict_product_univ.items():
        if key not in dict_product:
            dict_product[key] = val
        else:
            prev_lst = dict_product[key]
            dict_product[key] = list(set(prev_lst) | set(val)) 

    return dict_reactant, dict_product



def match_met_direction_innate(cobra_model):
    dict_reactant_model, dict_product_model = _extract_met_direction(cobra_model)
    return dict_reactant_model, dict_product_model



def run_MOMAs(model, biomass_rxn, target_rxn, num_pts=10):
    with model as m:
        wt_dist = pfba(m)
        wt_bio = wt_dist.fluxes[biomass_rxn]

    
    with model as m:
        m.objective = target_rxn
        m.objective_direction = 'max'
        max_const = m.slim_optimize()
        m.objective_direction = 'min'
        min_const = m.slim_optimize()

    results = {}
    for f_target in np.linspace(min_const, max_const, num_pts):
        with model as m:
            # m.reactions.get_by_id(target_rxn).bounds = (f_target * 0.95, f_target * 1.05)
            m.reactions.get_by_id(target_rxn).bounds = (f_target, f_target)
            f_bio = m.slim_optimize()
        with model as m:
            # m.reactions.get_by_id(target_rxn).bounds = (f_target * 0.95, f_target * 1.05)
            m.reactions.get_by_id(target_rxn).bounds = (f_target, f_target)
            # m.reactions.get_by_id(biomass_rxn).bounds = (f_bio * 0.95, f_bio * 1.05)
            m.reactions.get_by_id(biomass_rxn).lower_bound = f_bio * 0.95
            sol_moma = moma(m, solution=wt_dist, linear=True)
            results[f_target] = sol_moma.fluxes
    df = pd.DataFrame.from_dict(results)

    df_corr = df.abs().T.corr()
    df_cov = df.abs().T.cov()
    
    return df_corr, df_cov



def select_candidates(model, df_corr, df_cov, corr_threshold=0, cov_threshold=0.1):
    # df_pcorr = df_corr[df_corr > corr_threshold].dropna()
    df_pcov = df_cov[df_cov > cov_threshold].dropna()
    # pos_rxns = list(set(df_pcorr.index)&set(df_pcov.index))
    pos_rxns = list(df_pcov.index)

    # df_ncorr = df_corr[df_corr < -corr_threshold].dropna()
    df_ncov = df_cov[df_cov < -cov_threshold].dropna()
    # neg_rxns = list(set(df_ncorr.index)&set(df_ncov.index))
    neg_rxns = list(df_ncov.index)
    
    print('Number of positive candidate reactions: %d'%(len(pos_rxns)))
    print('Number of negative candidate reactions: %d'%(len(neg_rxns)))

    dict_met_scores = {}
    for met in model.metabolites:
        rxn_candidates = [rxn.id for rxn in met.reactions]
        rxn_pos_candidates = list(set(pos_rxns) & set(rxn_candidates))
        rxn_neg_candidates = list(set(neg_rxns) & set(rxn_candidates))


        score_pos = df_pcov[rxn_pos_candidates].abs().sum()
        score_neg = df_ncov[rxn_neg_candidates].abs().sum()
        score_norm = (score_pos - score_neg)/np.sqrt(1 + len(rxn_candidates))

        dict_met_scores[met.id] = (score_pos, score_neg, score_norm)

    return dict_met_scores




def _make_candidate_reaction_sets(dict_met_scores, exclude_met_dir): 
    ex_metabolites = []
    with open(exclude_met_dir, 'r') as fp:
        for line in fp:
            ex_metabolites.append(line.strip())
            
    met_pos = []
    met_neg = []
    met_neut = []

    for met_id, (score_pos, score_neg, score_norm) in dict_met_scores.items():
        if met_id[-1] != 'c':
            continue
        elif met_id in ex_metabolites:
            continue

        if score_norm > 0:
            met_pos.append(met_id)
        elif score_norm < 0:
            met_neg.append(met_id)
        else:
            met_neut.append(met_id)

    return met_pos, met_neg, met_neut



def find_reactions(dict_reactant, dict_product, dict_met_scores, exclude_met_dir):

    met_pos, met_neg, _ = _make_candidate_reaction_sets(
        dict_met_scores, exclude_met_dir
        )

    target_over = {}
    target_down = {}
    for met1 in met_neg:
        SoC_neg = abs(dict_met_scores[met1][-1])
        Un_consumes = dict_reactant[met1]
        Un_produces = dict_product[met1]

        for met2 in met_pos:
            result_over = {}
            result_down = {}

            SoC_pos = abs(dict_met_scores[met2][-1])
            Up_produces = dict_product[met2]

            R_np = set(Up_produces) & set(Un_consumes)
            if len(R_np) != 0:
                for rxn_id in R_np:
                    result_over[rxn_id] = SoC_pos - SoC_neg

            Up_consumes = dict_reactant[met2]
            R_pn = set(Up_consumes) & set(Un_produces)
            if len(R_pn) != 0:
                for rxn_id in R_pn:
                    result_down[rxn_id] = SoC_neg - SoC_pos

            target_over[(met1, met2)] = result_over
            target_down[(met1, met2)] = result_down

    return target_over, target_down



def save_result(target_dict, model_reactions, over=True):
    if over:
        target_type = 'Up'
    else:
        target_type = 'Down'
    df_list = []
    for key, val in target_dict.items():
        met1, met2 = key
        for rxn_id, soc in val.items():
            df_list.append([rxn_id, met1, met2, soc, target_type])

    if over:
        columns_names = ['Target_reactions', 'Negative_metabolite', 'Positive_metabolite', 'SoC', 'Regulation_type']
    else:
        columns_names = ['Target_reactions', 'Positive_metabolite', 'Negative_metabolite', 'SoC', 'Regulation_type']

    if len(df_list) == 0:
        df_list = pd.DataFrame(['', '', '', '', '']).T
        df_list.columns = columns_names
        df_list.drop(0)
        return df_list

    else:
        df_list = pd.DataFrame(df_list)
        df_list.columns = columns_names


    df_list = df_list.sort_values(['SoC'], ascending=False)
    df_list.index = np.arange(len(df_list))
    soc_max = df_list['SoC'].max()
    soc_min = df_list['SoC'].min()
    df_list['Normalized_SoC'] = (df_list['SoC'] - soc_min) / (soc_max-soc_min)
    df_list = df_list[df_list['Normalized_SoC'] >= 0.5]
    return df_list



def predict_ibridge_targets(model, exclude_met_dir, num_step=10):
    target_chem_id = '_'.join(model.id.split('_')[1:])
    target_chem_id = target_chem_id.split('_path_')[0]
    target_rxn = f'EXT_{target_chem_id}_c'

    for rxn in model.reactions:
        if rxn.objective_coefficient == 1:
            biomass_rxn = rxn.id

    with model as m:
        max_bio = m.slim_optimize()

        if isnan(max_bio):
            return []
        elif max_bio < 1e-6:
            return []
        
        m.objective = target_rxn
        max_prod = m.slim_optimize()

        if isnan(max_prod):
            return []
        elif max_prod < 1e-6:
            return []
    
    model_reactions = [rxn.id for rxn in model.reactions]
    dict_reactant, dict_product = match_met_direction_innate(model)

    df_corr, df_cov = run_MOMAs(model, biomass_rxn, target_rxn, num_step)
    df_corr = df_corr[target_rxn].loc[model_reactions]
    df_cov = df_cov[target_rxn].loc[model_reactions]

    dict_met_scores = select_candidates(model, df_corr, df_cov, 0.1, 0.1)
    target_over, target_down = find_reactions(dict_reactant, dict_product, dict_met_scores, exclude_met_dir)

    up_targets = save_result(target_over, model_reactions, over=True)
    down_targets = save_result(target_down, model_reactions, over=False)

    rxn_duplicated = set(up_targets['Target_reactions']) & set(down_targets['Target_reactions'])
    rxn_duplicated = list(rxn_duplicated)

    df1 = up_targets[~up_targets['Target_reactions'].isin(rxn_duplicated)]
    df2 = down_targets[~down_targets['Target_reactions'].isin(rxn_duplicated)]
    df = pd.concat([df1, df2], sort=False).reset_index(drop=True)

    return df