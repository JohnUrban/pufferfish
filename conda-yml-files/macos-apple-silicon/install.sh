######################################################
source ~/software/conda/source-conda-base.sh 

CONDA_SUBDIR=osx-64 conda env create -f sandbox-pufferfish.clean.yml
conda activate pufferfish-sandbox-clean
CONDA_SUBDIR=osx-64 conda install pip numpy scipy matplotlib pandas requests scikit-learn biopython bedtools
