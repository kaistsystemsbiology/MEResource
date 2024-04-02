import copy

from gurobipy import multidict, tuplelist, Model, quicksum, GRB
from yieldFinder.Simulator import Simulator


class GapFilling(Simulator):
    def __init__(self):
        self.threshold = 0.0001
        self.limit_num = 5
        self.num_cpu = 1
        self.identified_targets = []

    def run_GapFill(self, target_reaction, flux_constraints={}, universal_reactions=[], inf_flag=False, solution_number=100):
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
                                        ub=upper_boundary_constraints[each_reaction], name=each_reaction)
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)

        for each_reaction in universal_reactions:
            b_bool[each_reaction] = m.addVar(vtype=GRB.BINARY, name=each_reaction)

        m.update()

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] == (fplus[each_reaction] - fminus[each_reaction]))

        for each_reaction in model_reactions:
            m.addConstr((fplus[each_reaction] - fminus[each_reaction]) >= lower_boundary_constraints[each_reaction])
            m.addConstr((fplus[each_reaction] - fminus[each_reaction]) <= upper_boundary_constraints[each_reaction])

        for each_reaction in universal_reactions:
            m.addConstr((fplus[each_reaction] - fminus[each_reaction]) <= 1000.0 * b_bool[each_reaction])
            m.addConstr((fplus[each_reaction] - fminus[each_reaction]) >= -1000.0 * b_bool[each_reaction])
            m.addConstr(
                (fplus[each_reaction] - fminus[each_reaction]) + 1000.0 * (1 - b_bool[each_reaction]) >= epsilon)

        m.addConstr(quicksum((b_bool[each_reaction]) for each_reaction in universal_reactions) <= self.limit_number)
        # m.addConstr(quicksum((b_bool[each_reaction]) for each_reaction in universal_reactions) >= 1)

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

        m.setObjective(quicksum((b_bool[each_reaction]) for each_reaction in universal_reactions), GRB.MINIMIZE)

        attr_list = []
        for each_solution in range(solution_number):
            m.optimize()
            if m.status == GRB.Status.OPTIMAL:
                result_info[each_solution] = []
                ReactionFlux = {}
                for reaction in universal_reactions:
                    if b_bool[reaction].x > 0.5:
                        result_info[each_solution].append(reaction)
                # for each_reaction in result_info[each_solution]:
                #     attr_list.append(m.addConstr(b_bool[each_reaction] == 0))
                #     attr_list.append(m.addConstr((fplus[each_reaction] - fminus[each_reaction]) == 0))
                attr_list.append(m.addConstr(
                    quicksum(b_bool[rxn] for rxn in result_info[each_solution]) <= len(result_info[each_solution])-1
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
        

    def load_universal_model(self, universal_model):
        self.universal_model = universal_model

    def fill_gap(self, target_reaction):
        universal_model = self.universal_model
        # cobra_model = self.cobra_model

        cobra_reactions = [each_reaction.id for each_reaction in self.cobra_model.reactions]
        added_reactions = []
        universal_reactions = []

        for each_reaction in universal_model.reactions:
            if each_reaction.id not in cobra_reactions:
                added_reactions.append(each_reaction)
                universal_reactions.append(each_reaction.id)

        self.cobra_model.add_reactions(added_reactions)
        self.load_cobra_model(self.cobra_model)

        result_info = self.run_GapFill(target_reaction=target_reaction, universal_reactions=universal_reactions)

        return result_info

    def run_gap_filling(self, universal_model, metabolic_gap_model, target_reaction):
        self.load_universal_model(universal_model)
        self.load_cobra_model(copy.deepcopy(metabolic_gap_model))
        # self.load_cobra_model(metabolic_gap_model)

        result_info = self.fill_gap(target_reaction)
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