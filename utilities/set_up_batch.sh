#!/bin/bash

for i in $( seq ${1} ${2} ); do
    
    mkdir $i
    mkdir $i/h_anion
    mkdir $i/hcn
    mkdir $i/hehhe
    mkdir $i/fhf
    
    if [ -e "hcn.input.${i}" ]; then
        mv hcn.input.$i $i/hcn/hcn.input
        cp hcn.nbox* $i/hcn/
        cp hcn*fgh* $i/hcn/
    fi

    if [ -e "hehhe.input.${i}" ]; then
        mv hehhe.input.$i $i/hehhe/hehhe.input
        cp hehhe.nbox* $i/hehhe
        cp hehhe*fgh* $i/hehhe
    fi

    if [ -e "fhf.input.${i}" ]; then
        mv fhf.input.$i $i/fhf/fhf.input
        cp fhf.nbox* $i/fhf/
        cp fhf*fgh* $i/fhf/
    fi

    if [ -e "h_anion.input.${i}" ]; then
        mv h_anion.input.$i $i/h_anion/h_anion.input
        cp h_anion*fgh* $i/h_anion/
    fi

    mv var.$i $i/var.dat

done
