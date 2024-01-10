#!/bin/bash

cd $1

cd hcn

ln -sf hcn.output hcn.out
ln -sf hcn.input hcn.in

sbatch -J td.hcn.${1} -t 0:30:00 -n 1 -c 4 --mem=24GB $QCHEMDEVDIR/job_submit_script.sh hcn

cd ../

cd ../
