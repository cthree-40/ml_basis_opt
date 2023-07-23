#!/bin/bash

# Build training set

startdir=$1
finaldir=$2

for i in $( seq $startdir $finaldir )
do
    cd $i
    
    # Get parameters from hcn/hcn.input file
    $MLBASOPT/utilities/get_var_values_from_file.sh hcn/hcn.input > var.dat
    
    # Compute RMSE
    $MLBASOPT/utilities/compute_rmse_density_from_batchdir.pl > rmse.dat
    
    cd ../
done

for i in $( seq $startdir $finaldir )
do
    cat $i/rmse.dat >> training.dat
done

