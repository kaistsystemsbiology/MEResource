# Comprehensive evaluation of the capabilities of microbial cell factories
## Procedure

**Note**: 
This source code was developed in Linux, and has been tested in Ubuntu 16.04 with Python 3.6.

1. Clone the repository

        git clone https://github.com/kaistsystemsbiology/MEResource.git

2. Create and activate virtual environment

        conda env create -f environment.yml
        conda activate baebra

3. run resource calculation pipeline (for yeast850, use run_pipeline_with_fva.py) ( ~ 1~10 hours, depending on the model)

        python run_pipeilne.py -c eco -t CHEBI:16724
        python run_pipeilne_with_fva.py -c sce -t CHEBI:16724