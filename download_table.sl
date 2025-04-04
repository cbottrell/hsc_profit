#!/bin/bash -l
#SBATCH --account=pawsey0119
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=slurm/table.out
#SBATCH --error=slurm/table.err
#SBATCH --job-name=table
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=connor.bottrell@uwa.edu.au
#SBATCH --export=NONE
   
# Set working directory
# Should be in a writable path with some space, like /scratch
dir="${MYSOFTWARE}/work/Simulations/Scripts/Morphologies"
 
# Set the name of the python enviroment we will use
python_env="setonix_cpu"
 
# Load dependencies
# Note that the specific version will change over time

module load py-pandas/2.1.2
module load py-numpy/1.26.1
module load py-cython/3.0.4
module load py-mpi4py/3.1.5-py3.11.6
module load py-scikit-learn/1.3.2
module load py-h5py/3.8.0
module load py-matplotlib/3.8.1
module load py-pip/23.1.2-py3.11.6

# Activate virtual environment with non-modular packages
source ${MYSOFTWARE}/python-env/${python_env}/bin/activate

python download_table.py
