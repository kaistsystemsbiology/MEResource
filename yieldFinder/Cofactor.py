import copy

from gurobipy import *
from gurobipy import multidict, Model, quicksum, tuplelist, GRB
from yieldFinder.Simulator import Simulator


class SwapCofactor(Simulator):
    def __init__(self):
        self.threshold = 0.0001
        self.limit_num = 5
        self.num_cpu = 1

    def run_SwapCofactor(self, target_reaction, candidate_reactions={}, solution_number=100):
        added_reactions = []
        result_info = {}
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        lower_boundary_constraints = self.lower_boundary_constraints
        upper_boundary_constraints = self.upper_boundary_constraints

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
        m.reset()

        m.params.Threads = self.num_cpu
        m.update()

        # create variables
        epsilon = 0.001

        v = {}
        fplus = {}
        fminus = {}
        b_bool = {}
        target_bool = {}

        for each_reaction in model_reactions:
            v[each_reaction] = m.addVar(lb=lower_boundary_constraints[each_reaction],
                                        ub=upper_boundary_constraints[each_reaction], name=each_reaction)
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)

        for each_reaction in candidate_reactions:
            target_bool[each_reaction] = m.addVar(vtype=GRB.BINARY, name=each_reaction+'_target')
            b_bool[each_reaction] = m.addVar(vtype=GRB.BINARY, name=each_reaction)
            swap_reactions = candidate_reactions[each_reaction]
            for swap_rxn in swap_reactions:
                b_bool[swap_rxn] = m.addVar(vtype=GRB.BINARY, name=swap_rxn)

        m.update()

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] == (fplus[each_reaction] - fminus[each_reaction]))

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] >= lower_boundary_constraints[each_reaction])
            m.addConstr(v[each_reaction]  <= upper_boundary_constraints[each_reaction])

        for each_reaction in candidate_reactions:
            m.addConstr(v[each_reaction]  <= 1000.0 * b_bool[each_reaction])
            m.addConstr(v[each_reaction]  >= -1000.0 * b_bool[each_reaction])
            swap_reactions = candidate_reactions[each_reaction]
            for swap_rxn in swap_reactions:
                m.addConstr(v[swap_rxn]  <= 1000.0 * b_bool[swap_rxn])
                m.addConstr(v[swap_rxn]  >= -1000.0 * b_bool[swap_rxn])
            same_reactions = [each_reaction]
            same_reactions += swap_reactions
            m.addConstr(quicksum(b_bool[each_rxn] for each_rxn in same_reactions) == 1)
            m.addConstr(quicksum(b_bool[each_rxn] for each_rxn in swap_reactions) <= target_bool[each_reaction])

        m.addConstr(quicksum(target_bool[each_reaction] for each_reaction in target_bool) <= self.limit_number)

        m.addConstr((fplus[target_reaction] - fminus[target_reaction]) >= self.threshold)

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
            m.addConstr(quicksum(b_bool[rxn] for rxn in prev_targets) <= len(prev_targets)-1)
        m.update()

        m.setObjective(quicksum((target_bool[each_reaction]) for each_reaction in target_bool), GRB.MINIMIZE)

        attr_list = []
        for each_solution in range(solution_number):
            m.optimize()
            if m.status == GRB.Status.OPTIMAL:
                result_info[each_solution] = []
                ReactionFlux = {}
                for reaction in target_bool:
                    if target_bool[reaction].x > 0.5:
                        for candiddate_swap in candidate_reactions[reaction]:
                            if b_bool[candiddate_swap].x > 0.5:
                                result_info[each_solution].append(candiddate_swap)
                for each_reaction in result_info[each_solution]:
                    attr_list.append(m.addConstr(b_bool[each_reaction] == 0))
                    attr_list.append(m.addConstr((fplus[each_reaction] - fminus[each_reaction]) == 0))
                m.update()

            else:
                m.reset()

                return result_info

        for attr in attr_list:
            removeConstraintIndex = m.getConstrs().index(attr)
            m.remove(m.getConstrs()[removeConstraintIndex])

        m.reset()
        return result_info
        

    def run_swap_cofactor(self, metabolic_gap_model, target_reaction, candidate_reactions):
        cobra_model = self.cobra_model

        cobra_reactions = [each_rxn.id for each_rxn in cobra_model.reactions]
        added_reactions = []

        for each_reaction in candidate_reactions:
            generated_rxns = candidate_reactions[each_reaction]
            added_reactions += generated_rxns


        result_info = self.run_SwapCofactor(
            target_reaction=target_reaction, 
            candidate_reactions=candidate_reactions
        )
        
        return result_info

    def set_threshold(self, threshold):
        self.threshold = threshold

    def set_limit_number(self, limit_number):
        self.limit_number = limit_number

    def set_cpu_number(self, num_cpu):
        self.num_cpu = num_cpu

    def set_additional_constr(self, identified_targets):
        self.identified_targets = identified_targets
        return





class CoCofactor(Simulator):
    def __init__(self):
        self.threshold = 0.0001
        self.limit_num = 5
        self.num_cpu = 1

    def run_CoCofactor(self, target_reaction, candidate_reactions={}, solution_number=100):
        result_info = {}
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        lower_boundary_constraints = self.lower_boundary_constraints
        upper_boundary_constraints = self.upper_boundary_constraints

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
        m.reset()

        m.params.Threads = self.num_cpu
        m.update()

        # create variables

        v = {}
        fplus = {}
        fminus = {}
        b_bool = {}
        target_bool = {}
        direction_bool = {}

        for each_reaction in model_reactions:
            v[each_reaction] = m.addVar(lb=lower_boundary_constraints[each_reaction],
                                        ub=upper_boundary_constraints[each_reaction], name=each_reaction)
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)

        for original_rxn, swapped_rxns in candidate_reactions.items():
            target_bool[original_rxn] = m.addVar(vtype=GRB.BINARY, name=original_rxn+'_target')
            # b_bool[original_rxn] = m.addVar(vtype=GRB.BINARY, name=original_rxn)
            for swap_rxn in swapped_rxns:
                b_bool[swap_rxn] = m.addVar(vtype=GRB.BINARY, name=swap_rxn)

        for original_rxn, swapped_rxns in candidate_reactions.items():
            direction_bool[original_rxn] = m.addVar(vtype=GRB.BINARY, name=original_rxn+'_target_direction')
            for swap_rxn in swapped_rxns:
                direction_bool[swap_rxn] = m.addVar(vtype=GRB.BINARY, name=swap_rxn+'_direction')

        m.update()

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] == (fplus[each_reaction] - fminus[each_reaction]))

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] >= lower_boundary_constraints[each_reaction])
            m.addConstr(v[each_reaction] <= upper_boundary_constraints[each_reaction])


        for original_rxn, swapped_rxns in candidate_reactions.items():
            m.addConstr(v[original_rxn] >= -1000 * (1-direction_bool[original_rxn]))
            m.addConstr(v[original_rxn] <= 1000 * direction_bool[original_rxn])
            for swap_rxn in swapped_rxns:
                m.addConstr(v[swap_rxn] >= -1000 * (1-direction_bool[swap_rxn]))
                m.addConstr(v[swap_rxn] <= 1000 * direction_bool[swap_rxn])


        for original_rxn, swapped_rxns in candidate_reactions.items():
            # m.addConstr(v[original_rxn]  <= 1000.0 * b_bool[original_rxn])
            # m.addConstr(v[original_rxn]  >= -1000.0 * b_bool[original_rxn])
            for swap_rxn in swapped_rxns:
                m.addConstr(v[swap_rxn] <= 1000.0 * b_bool[swap_rxn])
                m.addConstr(v[swap_rxn] >= -1000.0 * b_bool[swap_rxn])

                # Let the swapped reaction and the original reaction have same direction (to prevent cycle)
                '''
                m.addConstr(v[swap_rxn] * v[original_rxn] >= 0)
                Linearize the above equation as below.
                '''
                m.addConstr(direction_bool[swap_rxn] - direction_bool[original_rxn] == 0)

            # m.addConstr(quicksum(b_bool[each_rxn] for each_rxn in swapped_rxns) == 1)
            m.addConstr(quicksum(b_bool[each_rxn] for each_rxn in swapped_rxns) <= target_bool[original_rxn])

        m.addConstr(quicksum(target_bool[each_reaction] for each_reaction in target_bool) <= self.limit_number)

        m.addConstr((fplus[target_reaction] - fminus[target_reaction]) >= self.threshold)

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
            m.addConstr(quicksum(b_bool[rxn] for rxn in prev_targets) <= len(prev_targets)-1)
        m.update()

        m.setObjective(quicksum((target_bool[each_reaction]) for each_reaction in target_bool), GRB.MINIMIZE)

        attr_list = []
        for each_solution in range(solution_number):
            m.optimize()
            if m.status == GRB.Status.OPTIMAL:
                result_info[each_solution] = []
                for reaction in target_bool:
                    if target_bool[reaction].x > 0.5:
                        for candiddate_swap in candidate_reactions[reaction]:
                            if b_bool[candiddate_swap].x > 0.5:
                                result_info[each_solution].append(candiddate_swap)
                for each_reaction in result_info[each_solution]:
                    attr_list.append(m.addConstr(b_bool[each_reaction] == 0))
                    attr_list.append(m.addConstr((fplus[each_reaction] - fminus[each_reaction]) == 0))
                m.update()

            else:
                m.reset()

                return result_info

        for attr in attr_list:
            removeConstraintIndex = m.getConstrs().index(attr)
            m.remove(m.getConstrs()[removeConstraintIndex])

        m.reset()
        return result_info
        

    def run_coutilize_cofactor(self, target_reaction, candidate_reactions):
        added_reactions = []

        for each_reaction in candidate_reactions:
            generated_rxns = candidate_reactions[each_reaction]
            added_reactions += generated_rxns

        result_info = self.run_CoCofactor(
            target_reaction=target_reaction, 
            candidate_reactions=candidate_reactions
        )
        
        return result_info

    def set_threshold(self, threshold):
        self.threshold = threshold

    def set_limit_number(self, limit_number):
        self.limit_number = limit_number

    def set_cpu_number(self, num_cpu):
        self.num_cpu = num_cpu

    def set_additional_constr(self, identified_targets):
        self.identified_targets = identified_targets
        return
