import os
import re
import copy
import logging
import warnings
from math import isnan 

import numpy as np
import pandas as pd

from gurobipy import multidict, tuplelist, Model, quicksum, GRB
from cobra.flux_analysis.loopless import add_loopless

from tqdm import tqdm

from baebra.calculator import _make_anaerobic_condition
from baebra.calculator import calculate_yield


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


def set_model_condition(model, c_source, prev_c_source='glc__D', air='aerobic', maintenance=None):
    if air == 'aerobic':
        model.reactions.EX_o2_e.lower_bound = -1000
    elif air == 'anaerobic':
        if 'yeast' in model.id:
            model = _make_anaerobic_condition(model)
        model.reactions.EX_o2_e.lower_bound = 0
    elif air == 'microaerobic':
        if 'yeast' in model.id:
            model = _make_anaerobic_condition(model)
        model.reactions.EX_o2_e.lower_bound = -0.5
    else:
        raise Exception("Wrong air condition")
    
    c_source_rxn_id = f'EX_{c_source}_e'
    prev_c_source_rxn_id = f'EX_{prev_c_source}_e'

    p_C_num = re.compile('C(\d*)')

    c_source_chem = model.metabolites.get_by_id(c_source+'_e')
    prev_c_source_chem = model.metabolites.get_by_id(prev_c_source+'_e')

    c_source_C_num = float(p_C_num.match(c_source_chem.formula).group(1))
    prev_c_source_C_num = float(p_C_num.match(prev_c_source_chem.formula).group(1))

    carbon_ratio = c_source_C_num / prev_c_source_C_num

    prev_c_source_rxn = model.reactions.get_by_id(prev_c_source_rxn_id)
    prev_lb = prev_c_source_rxn.lower_bound
    prev_c_source_rxn.lower_bound = 0
    model.reactions.get_by_id(c_source_rxn_id).lower_bound = prev_lb / carbon_ratio

    if maintenance:
        model.reactions.get_by_id(maintenance).bounds = (0, 0)

    return model


def predict_heterologous_reactions(target_model, universal_model, 
                                   c_source='EX_glc__D_e',
                                   limit_reaction_num=1, num_cpu=1,
                                   num_step=51):
    
    target_chem_id = '_'.join(target_model.id.split('_')[1:])
    target_chem_id = target_chem_id.split('_path_')[0]
    target_reaction = f'EXT_{target_chem_id}_c'

    
    expanded_model = copy.deepcopy(target_model)
    original_rxns = [rxn.id for rxn in target_model.reactions]

    with target_model as m:
        m.objective = target_reaction
        max_f = m.slim_optimize()
        if isnan(max_f) or abs(max_f) < 1e-3:
            max_f = 0
            return False
        m.reactions.get_by_id(target_reaction).bounds = (max_f, max_f)
        m.objective = c_source
        m.objective_direction = 'max'
        min_c = m.slim_optimize()

        if isnan(min_c) or abs(min_c) < 1e-3:
            min_c = 0
            print('Unable to uptake carbon source')
            return False

    prev_theoretical_yield = abs(max_f/min_c)
    #
    original_max_yield = max_f
        

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
            if each_reaction.id in original_rxns:
                rm_reactions.append(each_reaction)
            pass
        else:
            rm_reactions.append(each_reaction)

    print('No. of removed reactions in universal model: %d'%(len(rm_reactions)))
    universal_model.remove_reactions(rm_reactions)

    added_reactions, remove_reactions = _filter_same_reactions(expanded_model, universal_model)
    print('No. of added reactions: %d'%(len(added_reactions)))
    universal_model.remove_reactions(remove_reactions)
    expanded_model.add_reactions(added_reactions)

    with expanded_model as m:
        m.objective = target_reaction
        max_f = m.slim_optimize()
        if isnan(max_f) or abs(max_f) < 1e-3:
            max_f = 0
            return False
        m.reactions.get_by_id(target_reaction).bounds = (max_f, max_f)
        m.objective = c_source
        m.objective_direction = 'max'
        min_c = m.slim_optimize()

        if isnan(min_c) or abs(min_c) < 1e-3:
            min_c = 0
            print('Unable to uptake carbon source in the expanded model')
            return False
        
    theoretical_yield = abs(max_f/min_c)
    if theoretical_yield - prev_theoretical_yield < 1e-3:
        print('no yield improvement')
        return False
    elif theoretical_yield < 1e-3:
        print('zero yield')
        return False
    
    print('Improved max yield: %0.6f'%theoretical_yield)

    start_f = original_max_yield
    limit_f = original_max_yield * 1.5
    print('From %s to %s'%(start_f, limit_f))
    yield_list = np.linspace(start_f, limit_f, num_step)
    target_identified = False
    identified_targets = []

    obj = GapFilling()
    obj.set_limit_number(limit_reaction_num)
    obj.set_cpu_number(num_cpu)
    obj.set_univ_model(universal_model)


    for each_constraint in yield_list[1:]:
        print('Yield constraint : %0.3fX\t%0.6f'%(each_constraint/original_max_yield, each_constraint))
        obj.set_threshold(each_constraint)
        obj.set_additional_constr(identified_targets)
        result_info = obj.run_gap_filling(target_model, target_reaction)

        if len(result_info) == 0:
            break
        for each_solution in result_info:
            target_reactions = result_info[each_solution]
            if len(target_reactions) > 0:
                target_identified = True
            target_reactions.sort()
            identified_targets.append(target_reactions)

    if target_identified:
        tmp = []
        for target_reactions in identified_targets:
            tmp.append(';'.join(target_reactions))
        tmp = set(tmp)
        target_deduplicated = []
        for targets in tmp:
            targets = targets.split(';')
            targets.sort()
            target_deduplicated.append(targets)
        return target_deduplicated
    else:
        return target_identified
    


class GapFilling(object):
    def __init__(self):
        self.threshold = 0.0001
        self.limit_num = 5
        self.num_cpu = 1
        self.identified_targets = []
        self.universal_model = None


    def run_GapFill(self, target_reaction, universal_reactions=[], inf_flag=False, solution_number=100):
        result_info = {}
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        lower_boundary_constraints = self.lower_boundary_constraints
        upper_boundary_constraints = self.upper_boundary_constraints

        if not inf_flag:
            for key in lower_boundary_constraints:
                if lower_boundary_constraints[key] == float("-inf"):
                    lower_boundary_constraints[key] = -1000.0

            for key in upper_boundary_constraints:
                if upper_boundary_constraints[key] == float("inf"):
                    upper_boundary_constraints[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model('Gap filling')
        m.setParam('OutputFlag', 0)
        # m.setParam('DualReductions', 0)
        m.reset()

        m.params.Threads = self.num_cpu
        m.update()

        # create variables
        epsilon = 0.001

        v = {}
        fplus = {}
        fminus = {}
        b_bool = {}

        for each_reaction in model_reactions:
            v[each_reaction] = m.addVar(lb=lower_boundary_constraints[each_reaction],
                                        ub=upper_boundary_constraints[each_reaction], 
                                        name=each_reaction
                                )
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)

        for each_reaction in universal_reactions:
            b_bool[each_reaction] = m.addVar(vtype=GRB.BINARY, name=each_reaction)

        m.update()

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] == (fplus[each_reaction] - fminus[each_reaction]))

        for each_reaction in model_reactions:
            m.addConstr(
                (fplus[each_reaction] - fminus[each_reaction]) \
                      >= lower_boundary_constraints[each_reaction]
            )
            m.addConstr(
                (fplus[each_reaction] - fminus[each_reaction]) \
                    <= upper_boundary_constraints[each_reaction]
            )

        for each_reaction in universal_reactions:
            m.addConstr(
                (fplus[each_reaction] - fminus[each_reaction]) \
                    <= 1000.0 * b_bool[each_reaction]
            )
            m.addConstr(
                (fplus[each_reaction] - fminus[each_reaction]) \
                    >= -1000.0 * b_bool[each_reaction]
            )
            m.addConstr(
                (fplus[each_reaction] - fminus[each_reaction]) \
                    + 1000.0 * (1 - b_bool[each_reaction]) >= epsilon
            )

        m.addConstr(
            quicksum((b_bool[each_reaction]) for each_reaction in universal_reactions)\
                  <= self.limit_number
        )

        m.addConstr((fplus[target_reaction] - fminus[target_reaction]) 
                    >= self.threshold
        )

        m.update()

        # Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(
                (fplus[each_reaction] - fminus[each_reaction]) * coffvalue[metabolite, each_reaction] for
                metabolite, each_reaction in
                pairs.select(each_metabolite, '*')) == 0)

        m.update()

        # neglect previously identified targets
        for prev_targets in self.identified_targets:
            m.addConstr(
                quicksum(b_bool[rxn] for rxn in prev_targets) \
                    <= len(prev_targets)-1
            )
        m.update()

        m.setObjective(
            quicksum((b_bool[each_reaction]) for each_reaction in universal_reactions),
            GRB.MINIMIZE
        )

        attr_list = []
        for each_solution in range(solution_number):
            m.optimize()
            if m.status == GRB.Status.OPTIMAL:
                result_info[each_solution] = []
                ReactionFlux = {}
                for reaction in universal_reactions:
                    if b_bool[reaction].x > 0.5:
                        result_info[each_solution].append(reaction)
                attr_list.append(m.addConstr(
                    quicksum(
                        b_bool[rxn] for rxn in result_info[each_solution]) \
                            <= len(result_info[each_solution]) - 1
                ))
                m.update()

            else:
                m.reset()
                return result_info

        for attr in attr_list:
            removeConstraintIndex = m.getConstrs().index(attr)
            m.remove(m.getConstrs()[removeConstraintIndex])

        m.reset()
        return result_info
        


    def run_gap_filling(self, metabolic_gap_model, target_reaction):
        self.load_cobra_model(copy.deepcopy(metabolic_gap_model))

        universal_model = self.universal_model

        cobra_reactions = [rxn.id for rxn in self.cobra_model.reactions]
        added_reactions = []
        universal_reactions = []

        for each_reaction in universal_model.reactions:
            if each_reaction.id not in cobra_reactions:
                added_reactions.append(each_reaction)
                universal_reactions.append(each_reaction.id)
        self.cobra_model.add_reactions(added_reactions)

        self.load_cobra_model(self.cobra_model)

        result_info = self.run_GapFill(
            target_reaction=target_reaction, 
            universal_reactions=universal_reactions
        )

        return result_info



    def load_cobra_model(self, cobra_model):
        self.cobra_model = cobra_model
        model = cobra_model
        model_metabolites = []
        model_reactions = []
        model_genes = []
        lower_boundary_constraints = {}
        upper_boundary_constraints = {}
        objective_reaction = ''
        for each_metabolite in model.metabolites:
            model_metabolites.append(each_metabolite.id)
        
        model_genes = [each_gene.id for each_gene in model.genes]
        
        Smatrix = {}

        for each_reaction in model.reactions:
            if each_reaction.objective_coefficient == 1.0:
                objective_reaction = each_reaction.id

            reactant_list = list(each_reaction.reactants)
            reactant_coff_list = list(each_reaction.get_coefficients(reactant_list))
            product_list = list(each_reaction.products)
            product_coff_list = list(each_reaction.get_coefficients(product_list))

            for i in range(len(reactant_list)):
                Smatrix[(reactant_list[i].id, each_reaction.id)] = reactant_coff_list[i]

            for i in range(len(product_list)):
                Smatrix[(product_list[i].id, each_reaction.id)] = product_coff_list[i]

            model_reactions.append(each_reaction.id)
            lb = each_reaction.lower_bound
            ub = each_reaction.upper_bound
            if lb < -1000.0:
                lb = float('-inf')
            if ub > 1000.0:
                ub = float('inf')
            lower_boundary_constraints[each_reaction.id] = lb
            upper_boundary_constraints[each_reaction.id] = ub

        self.model_metabolites = model_metabolites
        self.model_reactions = model_reactions
        self.model_genes = model_genes
        self.Smatrix = Smatrix
        self.lower_boundary_constraints = lower_boundary_constraints
        self.upper_boundary_constraints = upper_boundary_constraints
        self.objective = objective_reaction

        return
    

    def set_threshold(self, threshold):
        self.threshold = threshold

    def set_limit_number(self, limit_number):
        self.limit_number = limit_number

    def set_cpu_number(self, num_cpu):
        self.num_cpu = num_cpu

    def set_additional_constr(self, identified_targets):
        self.identified_targets = identified_targets
        return
    
    def set_univ_model(self, universal_model):
        self.universal_model = universal_model
    




def validate_target(target_model, universal_model, target_identified, c_source, loopless=True):
    valid_targets = {}
    target_chem_id = '_'.join(target_model.id.split('_')[1:])
    target_chem_id = target_chem_id.split('_path_')[0]
    target_reaction = f'EXT_{target_chem_id}_c'
    for reactions in tqdm(target_identified):

        added_reactions = []
        for each_reaction in reactions:
            obj = universal_model.reactions.get_by_id(each_reaction)
            added_reactions.append(obj)

        temp_model = copy.deepcopy(target_model)
        temp_model.add_reactions(added_reactions)
        add_loopless(temp_model)

        temp_model.objective = target_reaction
        max_f = temp_model.slim_optimize()
        if isnan(max_f) or abs(max_f) < 1e-3:
            max_f = 0.0

        temp_model.reactions.get_by_id(target_reaction).bounds = (max_f, max_f)
        temp_model.objective = c_source
        temp_model.objective_direction = 'max'
        min_c = temp_model.slim_optimize()

        if isnan(min_c) or abs(min_c) < 1e-3:
            new_yield = 0.0
        else:
            new_yield = abs(max_f/min_c)

        valid_targets[';'.join(reactions)] = new_yield

    return valid_targets