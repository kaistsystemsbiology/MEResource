# Baebra

## Genome-scale metabolic model reconstruction and yield calculation for target chemical
## Procedure

**Note**: 
This source code was developed in Linux, and has been tested in Ubuntu 16.06 with Python 3.7.

1. Clone the repository

        git clone https://anlito@bitbucket.org/anlito/baebra.git

2. Create and activate virtual environment

        conda env create -f environment.yml
        conda activate baebra

3. Generate pathway introduced models

        python ModelGenerator.py -i iML1515

4. Calculate yields for the models

        python Yield_calculation.py -i iML1515

5. Identify cofactor alterations for the improvement of the maximum theoretical yield

        python runCofactorFinder.py -o ./results/results_iML1515_aerobic_glc__D_Cofactor -i1 ./TargetChemicalModels/iML1515/iML1515_ac.xml -i2 iML1515 -cpu 3 -n 1 -c glc__D -o2 True

6. Identify heterologous reactions for the improvement of the maximum theoretical yield

        python runReactionFinder.py -o ./results/results_iML1515_aerobic_glc__D_Reaction -i1 ./TargetChemicalModels/iML1515/iML1515_3php.xml -i2 iML1515 -cpu 3 -n 3 -c glc__D -o2 True