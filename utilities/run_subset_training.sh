#!/bin/bash

if [ $# -eq 0 ]
then
    echo "Usage: run_subset_training.sh [number of subsets]"
    exit 1
fi

if [ -z "$1" ]
then
    echo "Missing first argument."
    exit 1
fi
nbatch=$1

# Split up training data
$MLBASOPT/utilities/setup_subset_training.pl training.dat ${nbatch}

for i in $( seq 1 $nbatch )
do
    # Make subset directories
    mkdir ${i}
    cp ./* ${i}/
    cp training.dat_${i} ${i}/training.dat
    rm ${i}/training.dat_*
    cp ${i}/training.dat ${i}/testing.dat

done

    
