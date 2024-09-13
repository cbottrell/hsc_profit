#!/bin/bash -l
#SBATCH --account=pawsey0119
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=12:00:00
#SBATCH --output=slurm/jupyter.out
# #SBATCH --error=slurm/jupyter.err
#SBATCH --job-name=jupyter
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=bottrell
#SBATCH --export=NONE
   
# Set our working directory
# Should be in a writable path with some space, like /scratch
dir="${MYSOFTWARE}/work/Simulations/Scripts/Morphologies"
 
# Set the name of the python enviroment we will use
conda_env="setonix_cpu"
 
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
source ${MYSOFTWARE}/python-env/${conda_env}/bin/activate
 
# Modules 
module load r/4.3.0
module load opencv/4.8.0

export PATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/bin:$PATH
export LD_LIBRARY_PATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/lib64:$LD_LIBRARY_PATH
export CPATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/include:$CPATH
export PKG_CONFIG_PATH=/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/libtiff-4.5.1-rbf2zukoijn4rswztoxpzdrabkzqwgi3/lib64/pkgconfig:$PKG_CONFIG_PATH
  
# You should not need to edit the lines below
   
# Prepare the working directory
mkdir -p ${dir}
cd ${dir}
  
# Get the hostname
# We'll set up an SSH tunnel to connect to the Juypter notebook server
host=$(hostname)
   
# Set the port for the SSH tunnel
# This part of the script uses a loop to search for available ports on the node;
# this will allow multiple instances of GUI servers to be run from the same host node
port="9100"
pfound="0"
while [ $port -lt 65535 ] ; do
  check=$( ss -tuna | awk '{print $4}' | grep ":$port *" )
  if [ "$check" == "" ] ; then
    pfound="1"
    break
  fi
  : $((++port))
done
if [ $pfound -eq 0 ] ; then
  echo "No available communication port found to establish the SSH tunnel."
  echo "Try again later. Exiting."
  exit
fi
  
   
echo "*****************************************************"
echo "Setup - from your laptop do:"
echo "ssh -N -f -L ${port}:${host}:${port} $USER@$PAWSEY_CLUSTER.pawsey.org.au"
echo "*****"
echo "The launch directory is: $dir"
echo "*****************************************************"
echo ""
echo "*****************************************************"
echo "Terminate - from your laptop do:"
echo "kill \$( ps x | grep 'ssh.*-L *${port}:${host}:${port}' | awk '{print \$1}' )"
echo "*****************************************************"
echo ""
    
#Launch the notebook
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
    jupyter lab \
  --no-browser \
  --port=${port} --ip=0.0.0.0 \
  --notebook-dir=${dir}
