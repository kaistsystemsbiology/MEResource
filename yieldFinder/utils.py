from argparse import ArgumentParser


def argument_parser(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = ArgumentParser()
    parser.add_argument('-o', '--output_dir', 
                        required=True, help='Output directory')
    parser.add_argument('-i1', '--input_model', 
                        required=True, help='Input model directory')
    parser.add_argument('-i2', '--model_name', 
                        required=True, help='Input model name')
    parser.add_argument('-l', '--log_dir', required=False, 
                        default='yieldImproving.log', help='Log file directory')
    parser.add_argument('-n', '--num_reaction', required=False, type=int,
                        default=3, help='Limit of the number of reactions')
    parser.add_argument('-c', '--carbon_source', required=False, 
                        default='glc__D', 
                        choices=['glc__D', 'xyl__D', 'arab__L', 'fru', 'gal', 'glyc', 'sucr'],
                        help='Carbon source (glc__D/xyl__D/arab__L/fru/gal/glyc/sucr)')
    parser.add_argument('-o2', '--oxygen_availability', choices=['aerobic', 'anaerobic', 'microaerobic'],
                        required=False, default='aerobic',
                        help='Oxygen availability')
    parser.add_argument('-u', '--universal_model', required=False, 
                        default='./inputs/annotated_universal_model.json', 
                        help='Universal model directory')
    parser.add_argument('-m', '--maintenance', required=False, 
                        default='ATPM', 
                        help='Non growth associated energy requirement reaction')
    parser.add_argument('-cpu', '--cpu_number', required=False, type=int,
                        default=1, 
                        help='Number of cpu thread to use')
    return parser