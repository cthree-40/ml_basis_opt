#!/bin/bash

cd $1

cd fhf

ln -sf fhf.output fhf.out
ln -sf fhf.input fhf.in

sbatch -J td.fhf.${1} -t 00:30:00 -n 1 -c 12 --mem=24GB $QCHEMDEVDIR/job_submit_script.sh fhf

cd ../

cd ../
