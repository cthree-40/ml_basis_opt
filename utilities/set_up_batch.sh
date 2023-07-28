#!/bin/bash

for i in $( seq ${1} ${2} ); do
    
    mkdir $i
    mkdir $i/hcn
    mkdir $i/hehhe
    
    if [ -e "hcn.input.${i}" ]; then
        mv hcn.input.$i $i/hcn/hcn.input
        cp hcn.nbox* $i/hcn/
    fi

    if [ -e "hehhe.input.${i}" ]; then
        mv hehhe.input.$i $i/hehhe/hehhe.input
        cp hehhe.nbox* $i/hehhe
    fi

    mv var.$i $i/var.dat

done
