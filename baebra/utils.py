from argparse import ArgumentParser


def argument_parser(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = ArgumentParser()
    parser.add_argument('-i', '--input_model', 
                        required=True, help='Input model name')

    parser.add_argument('-l', '--log_dir', required=False, 
                        default='yieldImproving.log', help='Log file directory')
    parser.add_argument('-n', '--num_reaction', required=False, type=int,
                        default=3, help='Limit of the number of reactions')
    parser.add_argument('-u', '--universal_model', required=False, 
                        default='./input/annotated_universal_model.json', 
                        help='Universal model directory')
    parser.add_argument('-m', '--maintenance', required=False, 
                        default='ATPM', 
                        help='Non growth associated energy requirement reaction')
    parser.add_argument('-cpu', '--cpu_number', required=False, type=int,
                        default=1, 
                        help='Number of cpu thread to use')
    return parser




def argument_parser_pipeline(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = ArgumentParser()
    parser.add_argument('-b', '--base_config',
                        required=False, default='./config/base.yaml')
    parser.add_argument('-c', '--config', 
                        required=False, type=str, default=None)
    parser.add_argument('-i', '--input_model', 
                        required=False, help='Input model name', default=None)
    parser.add_argument('-t', '--target_chem', required=False,
                        help='Target chemical CHEBI ID', type=str, default=None)
    parser.add_argument('-o', '--output_dir', required=False,
                        help='Output directory', type=str, default=None)
    return parser