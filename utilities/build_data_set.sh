#!/bin/bash -l

# Build training set

startdir=$1
finaldir=$2
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

    $MLBASOPT/utilities/compute_penalty_function_from_batchdir.pl 1.0 0.0 0.0 1.0 > rmse.dat

    cd ../
done

rm -rf training.dat

for i in $( seq $startdir $finaldir )
do
    cat $i/rmse.dat >> training.dat
done

