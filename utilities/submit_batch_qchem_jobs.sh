#!/bin/bash -l
#SBATCH --job-name=dev_qchem_job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=24
#SBATCH --mem=150GB
#SBATCH --partition=compute
#SBATCH --account=uic317

# Set necessary environmental variables
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULE=dynamic
#export OMP_PROC_BIND=true
#export OMP_PLACES=threads

# Set up development qchem environment
source ~/.qcdevsetup
export QCTHREADS=$OMP_NUM_THREADS

echo "Submitting BATCH dev-QCHEM job."
echo "Current directory: $(pwd)"
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "OMP_SCHEDULE = $OMP_SCHEDULE"

echo "Start directory: ${1}"
echo "Final directory: ${2}"

molec=("fhf_2.58" "fhf_2.62" "fhf_2.66" "fhf_2.70" "fhf_2.74" "h2_0.74" "h_anion" "hehhe" "hcn" "hnc" "fhf")

for i in $( seq ${1} ${2} ); do
    cd $i

    # Loop over molecular systems
    for sys in "${molec[@]}"; do

	# If the input was created for this system, run Qchem
	if [ -e "${sys}/${sys}.input" ]; then

	    # Copy nbox files to directory
	    if [ -e "../${sys}.nbox_data" ]; then
		cp ${sys}.nbox_data ${sys}/
		cp ${sys}.nbox_npts ${sys}/
	    fi

	    cd ${sys}

            # Run QCHem
            qchem -nt $OMP_NUM_THREADS ${sys}.input > ${sys}.output
            cd ../
	    
	fi
	
    done

    cd ../
    
done
