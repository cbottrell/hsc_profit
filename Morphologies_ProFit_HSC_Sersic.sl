#!/bin/bash -l
#SBATCH --account=pawsey0119
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G 
#SBATCH --time=23:59:59
#SBATCH --array=0-15
#SBATCH --output=slurm/profit_%A_%a.out
#SBATCH --error=slurm/profit_%A_%a.err
#SBATCH --job-name=ProFit
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
 
# Modules 
module load r/4.3.0

# R library paths
export PATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/bin:$PATH
export LD_LIBRARY_PATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/lib64:$LD_LIBRARY_PATH
export CPATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/include:$CPATH
export PKG_CONFIG_PATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/lib64/pkgconfig:$PKG_CONFIG_PATH

# Job environment variables

export UNIVERSE='IllustrisTNG'
export SIMULATION='TNG100-1'
export DATABASE='TNG100_1'
export SNAPMIN=72
export SNAPMAX=91
export MSTARMIN=10

export TABLE='Morphologies_ProFit_HSC_Sersic'

export JOB_ARRAY_SIZE=$SLURM_ARRAY_TASK_COUNT
export JOB_ARRAY_INDEX=$SLURM_ARRAY_TASK_ID

echo $SLURM_JOB_NUM_NODES
echo $SLURM_NTASKS

srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c 1 -m block:block:block python -m mpi4py Morphologies_ProFit_HSC_Sersic.py
