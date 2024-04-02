import xml.etree.ElementTree as elemTree
from libchebipy import ChebiEntity
import pandas as pd


class Metabolite():
    def __init__(self, met_id):
        self.id = met_id
        self.charge = None
        self.formula = None
        self.name = None
        self.mass = None
        self.inchikey = None
        self.extractMetInfo()
        self.reactions = []
        
    def extractMetInfo(self):
        chebientity = ChebiEntity(self.id)
        self.charge = chebientity.get_charge()
        self.formula = chebientity.get_formula()
        self.name = chebientity.get_name()
        self.mass = chebientity.get_mass()
        self.inchikey = chebientity.get_inchi_key()
        self.inchi = chebientity.get_inchi()
        self.smiles = chebientity.get_smiles()
        
    def get_reactions(self):
        for rxn in self.reactions:
            metabolits = rxn.metabolites
            rxn_eqn = rxn.reaction
            print(f'{rxn.id}\t{rxn_eqn}')
            for met in rxn.metabolites:
                rxn_eqn = rxn_eqn.replace(f'({met.id})', met.name)
            print(f'\t{rxn_eqn}\n')
        

class Reaction():
    def __init__(self, rxn_id, eqn):
        self.id = rxn_id
        self.reactants = None
        self.products = None
        self.metabolites = None
        self.coefficients = None
        self.ec = None
        self.bounds = None

        self.reaction = eqn
        self.getReactionInfo(eqn)
        
    def get_coefficient(self, met_id):
        return self.coefficients(metabolite)
    
    def knock_out(self):
        self.bounds = (0.0, 0.0)
        
    def seeStoichimetry(self):
        reactant_list = [met.id for met in self.reactants]
        product_list = [met.id for met in self.products]

        met_list = reactant_list + product_list
        coeff_list = [str(self.coefficient[met_id]) for met_id in met_list]
        tmp_dict = {'RXN_ID':[self.id], 
                    'Metabolite':met_list, 
                    'Coefficient':coeff_list}
        table = pd.DataFrame.from_dict(tmp_dict, orient='index')
        return table

        
    def getCoefficient(self, eqn, direction='reactant'):
        coeff_map = {}
        metabolites = eqn.split(' + ')
        for met in metabolites:
            tmp = met.split(' ')
            if len(tmp) == 1:
                coeff = 1
                met = tmp[0]
            elif len(tmp) == 2:
                coeff = float(tmp[0])
                met = tmp[1]
            else:
                print(eqn)
                raise ImplementationError

            if direction == 'reactant':
                coeff_map[met] = -coeff
            elif direction == 'product':
                coeff_map[met] = coeff
        return coeff_map


    def getReactionInfo(self, eqn):
        if ' <=> ' in eqn:
            self.bounds = (-1000, 1000)
            self.lower_bound = -1000
            self.upper_bound = 1000
            sep = ' <=> '
        elif ' => ' in eqn:
            self.bounds = (0, 1000)
            self.lower_bound = 0
            self.upper_bound = 1000
            sep = ' => '
        else:
            raise
        reactants, products = eqn.split(sep)
        reactants_coeff_map = self.getCoefficient(reactants, direction='reactant')
        products_coeff_map = self.getCoefficient(products, direction='product')
        reactants = list(reactants_coeff_map.keys())
        products = list(products_coeff_map.keys())
        metabolites = reactants + products
        coefficients = {}
        for met in reactants:
            coefficients[met[1:-1]] = reactants_coeff_map[met]
        for met in products:
            coefficients[met[1:-1]] = products_coeff_map[met]
            
        self.reactants = reactants
        self.products = products
        self.metabolites = metabolites
        self.coefficient = coefficients
        return
    
    
    

def replaceReaction(reaction, met_dict):
    reaction.reactants = [met_dict[met[1:-1]] for met in reaction.reactants]
    reaction.products = [met_dict[met[1:-1]] for met in reaction.products]
    reaction.metabolites = reaction.reactants + reaction.products
    
    
    
def saveExtractedData(metabolite_dict, reaction_dict):
    
    met_name = {}
    met_formula = {}
    met_charge = {}
    met_inchikey = {}
    met_mass = {}
    met_inchi = {}
    met_smiles = {}

    for met in metabolite_dict.values():
        met_name[met.id] = met.name
        met_formula[met.id] = met.formula
        met_charge[met.id] = met.charge
        met_inchikey[met.id] = met.inchikey
        met_inchi[met.id] = met.inchi
        met_smiles[met.id] = met.smiles
        met_mass[met.id] = met.mass

    met_table = pd.DataFrame({'Name':met_name, 
                              'Formula':met_formula,
                              'Charge':met_charge, 
                              'Mass':met_mass, 
                              'InChi':met_inchi,
                              'Smiles':met_smiles, 
                              'InChiKey':met_inchikey})
    met_table = met_table.sort_index()
    met_table.to_csv('./Metabolite.csv')

    rxn_direction = {}
    rxn_ec = {}
    for rxn in reaction_dict.values():
        if rxn.lower_bound < 0:
            rxn_direction[rxn.id] = '<=>'
        else:
            rxn_direction[rxn.id] = '>'
        rxn_ec[rxn.id] = rxn.ec
    table_rxn = pd.DataFrame({'Direction':rxn_direction, 'EC_number':rxn_ec})
    table_rxn = table_rxn.sort_index()
    table_rxn.to_csv('./Reaction.csv')

    rxn_ids = [rxn for rxn in reaction_dict]
    rxn_ids.sort()

    with open('./Stoichiometry.csv', 'w') as fp:
        for rxn_id in rxn_ids:
            rxn = reaction_dict[rxn_id]
            reactant_list = [met.id for met in rxn.reactants]
            product_list = [met.id for met in rxn.products]

            met_list = reactant_list + product_list
            coeff_list = [str(rxn.coefficient[met_id]) for met_id in met_list]
            met_str = ','.join(met_list)
            coeff_str = ','.join(coeff_list)
            fp.write(f'RXN_ID,{rxn_id}\n')
            fp.write(f'Metabolite,{met_str}\n')
            fp.write(f'Coefficient,{coeff_str}\n')

    
    
def extractData(rheaDB_dir='./baebra/rhea-ebeye.xml'):
    tree = elemTree.parse(rheaDB_dir)
    root = tree.getroot()
    child = root.getchildren()[-1]

    extractedRxn = {}
    for entry in child:
        rxn_id = entry.attrib['id']
        ec = None
        eqn = None
        for ref in entry[1]:
            if ref.attrib['dbname'] == 'IntEnz':
                ec = ref.attrib['dbkey']
        for additional_item in entry[2]:
            if additional_item.attrib['name'] == 'chebiIdEquation':
                eqn = additional_item.text
        extractedRxn[rxn_id] = [ec, eqn]


    for rxn in extractedRxn:
        eqn = extractedRxn[rxn][1]
        if ' = ' in eqn:
            extractedRxn[rxn][1] = eqn.replace(' = ', ' <=> ')

    reaction_list = []
    noChEBI_list = []
    for rxn in extractedRxn:
        eqn =  extractedRxn[rxn][1]
        if 'no ChEBI ID' in eqn:
            noChEBI_list.append(rxn)
            continue
        reaction = Reaction(rxn, eqn)
        reaction.ec = extractedRxn[rxn][0]
        reaction_list.append(reaction)


    metabolite_str_list = []
    for rxn in reaction_list:
        for met in rxn.metabolites:
            metabolite_str_list.append(met[1:-1])
    metabolite_str_list = list(set(metabolite_str_list))

    metabolite_list = []
    for met in metabolite_str_list:
        metabolite_list.append(Metabolite(met))
#     print(f'Number of metasbolite extracted: {len(metabolite_list)}')

    met_dict = {met.id: met for met in metabolite_list}
    for rxn in reaction_list:
        replaceReaction(rxn, met_dict)


    for rxn in reaction_list:
        for met in rxn.metabolites:
            met.reactions.append(rxn)
            
    rxn_dict = {rxn.id: rxn for rxn in reaction_list}
    met_dict = {met.id: met for met in metabolite_list}
            
            
    Info = {'metabolite_dict': met_dict, 'reaction_dict':rxn_dict}
    return Info
    
    
    
if __name__=='__main__':
    rheaDB_dir = './rhea-ebeye.xml'
    Info = extractData(rheaDB_dir)
    metabolite_dict = Info['metabolite_dict']
    reaction_dict = Info['reaction_dict']
    saveExtractedData(metabolite_dict, reaction_dict)