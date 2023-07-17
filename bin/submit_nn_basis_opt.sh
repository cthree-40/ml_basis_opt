#!/bin/bash -l
#SBATCH --job-name=dev_qchem_job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=24
#SBATCH --mem=150GB

# Set up environmental variables
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULE=dynamic
export OMP_PROC_BIND=true
export OMP_PLACES=threads

# Load anaconda environment
module load miniconda

# Load QChem environment
source /home/clm96/.qcdevsetup
export QCTHREADS=$OMP_NUM_THREADS

echo "Submitting basis optimization job."
echo " "
echo "Current directory: $(pwd)"
echo " "
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "OMP_SCHEDULE = $OMP_SCHEDULE"
echo "OMP_PROC_BIND = $OMP_PROC_BIND"
echo "OMP_PLACES = $OMP_PLACES"
echo " "

# Load optimization environment
conda activate qi_opt

python3 /home/clm96/pi_project/software/basis_opt/malbon_optimizer/nn_basis_opt/source/basis_optimization.s_shell.py
