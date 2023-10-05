#!/bin/bash -l

# Build training set

if [ $# -eq 0 ]
then
    echo "Usage: build_data_set.sh [start] [final] [nstates]"
    exit 1
fi

if [ -z "$1" ]
then
    echo "Missing first agument."
    exit 1
fi
startdir=$1

if [ -z "$2" ]
then
    echo "Missing second argument."
    exit 1
fi
finaldir=$2

if [ -z "$3" ]
then
    echo "Missing third argument."
    exit 1
fi
nstates=$3



for i in $( seq $startdir $finaldir )
do
    cd $i
    
    # Get state information
    for j in $( ls ); do
        if [[ -d "${j}"  && -f "${j}/${j}.output" ]]; then
            cd $j
            python3 $MLBASOPT/utilities/get_states_from_file.py -n $nstates -m $j 
            cd ../
        fi
    done

    # Compute RMSEs
    $MLBASOPT/utilities/compute_rmse_density_from_batchdir.pl
    $MLBASOPT/utilities/compute_rmse_zeropoint_from_batchdir.pl
    $MLBASOPT/utilities/compute_rmse_excited_states_from_batchdir.pl
    $MLBASOPT/utilities/compute_rmse_all_states_from_batchdir.pl

    # [w_dens] [w_xst] [w_gst] [w_ast]
    # w_dens = Ground state density
    # w_xst  = Excited state energies (relative)
    # w_gst  = Ground state energy    (absolute)
    # w_ast  = All state energies     (absolute)
    $MLBASOPT/utilities/compute_penalty_function_from_batchdir.pl 0.0 0.0 0.0 0.0 > rmse.dat

    cd ../
done

rm -rf training.dat

for i in $( seq $startdir $finaldir )
do
    cat $i/rmse.dat >> training.dat
done

