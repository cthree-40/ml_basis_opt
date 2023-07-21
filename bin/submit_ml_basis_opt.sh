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

export NN_DIR="/home/clm96/pi_project/software/basis_opt/malbon_optimizer/ml_basis_opt"
export NN_SDIR="$NN_DIR/source"

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

# Look for input file (This is required for proper code function.)
if [ -f "basis_opt.input" ]; then
    sed '/### GLOBAL VAR FROM INPUT ###/ r basis_opt.input' $NN_SDIR/basis_optimization.s_shell.py > ./basis_optimization.s_shell.py 
    python3 ./basis_optimization.s_shell.py
else
    echo "Input file not found!"
    echo "Please create an input file to continue."
fi
