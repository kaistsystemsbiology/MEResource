# Baebra

## Genome-scale metabolic model reconstruction and yield calculation for target chemical
## Procedure

**Note**: 
This source code was developed in Linux, and has been tested in Ubuntu 16.06 with Python 3.7.

1. Clone the repository

        git clone https://github.com/kaistsystemsbiology/MEResource.git

2. Create and activate virtual environment

        conda env create -f environment.yml
        conda activate baebra

3. run resource calculation pipeline

        python run_pipeilne.py -c eco -t CHEBI:16724