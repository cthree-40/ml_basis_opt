#!/bin/bash -l
#SBATCH --job-name=dev_qchem_job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=24
#SBATCH --mem=150GB

# Set necessary environmental variables
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULE=dynamic
#export OMP_PROC_BIND=true
#export OMP_PLACES=threads

# Set up development qchem environment
source ~/.qcdevsetup
export QCTHREADS=$OMP_NUM_THREADS

source $QCHEMDEVDIR/qchem_module_load.sh

echo "Submitting BATCH dev-QCHEM job."
echo "Current directory: $(pwd)"
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "OMP_SCHEDULE = $OMP_SCHEDULE"

echo "Start directory: ${1}"
echo "Final directory: ${2}"

for i in $( seq ${1} ${2} ); do
    cd $i
    
    if [ -e "hehhe/hehhe.input" ]; then
        cd hehhe
        cp hehhe.nbox_npts.txt nbox_npts.txt
        cp hehhe.nbox_data.txt nbox_data.txt
        # Run QChem
        qchem -nt $OMP_NUM_THREADS hehhe.input > hehhe.output
    cd ../
    fi
    
    cd ../
done
