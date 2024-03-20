#!/bin/bash -l
#SBATCH --job-name=dev_qchem_job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=24
#SBATCH --mem=150GB

# Set up environmental variables
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULE=dynamic
#export OMP_PROC_BIND=true
#export OMP_PLACES=threads

export ML_DIR=${MLBASOPT}
export ML_SDIR="$ML_DIR/source"

# Load modules
source $ML_DIR/bin/ml_basis.module_load.della.sh

# Load QChem environment
source ~/.qcdevsetup
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

# Look for input file (This is required for proper code function.)
if [ -f "basis_opt.input" ]; then
    sed '/### GLOBAL VAR FROM INPUT ###/ r basis_opt.input' $ML_SDIR/basis_optimization.py > ./basis_optimization.py 
    python3 ./basis_optimization.py
else
    echo "Input file not found!"
    echo "Please create an input file to continue."
fi
