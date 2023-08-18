#!/bin/bash

jname=$1

if [[ $jname == "" ]]; then
    echo "Please enter slurm job number."
    exit 1
fi

cat slurm-${jname}.out | grep -v "outside bounds" | grep -v "variances" | grep -v "warnings" > output.log-${jname}
