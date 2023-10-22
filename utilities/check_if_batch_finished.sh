#!/bin/bash


startdir=$1
finaldir=$2

echo "Uncompleted jobs:"
for i in $( seq ${startdir} ${finaldir} ); do
    cd $i

    # Check FHF
    if [ -f fhf/fhf.output ]; then
        if ! grep -q "Thank you very much for using Q-Chem." fhf/fhf.output; then
            echo "${i}/fhf"
        fi
    else
        echo "${i}/fhf"
    fi

    # Check HeHHe+
    if [ -f hehhe/hehhe.output ]; then
        if ! grep -q "Thank you very much for using Q-Chem." hehhe/hehhe.output; then
            echo "${i}/hehhe"
        fi
    else
        echo "${i}/hehhe"
    fi

    # Check HCN
    if [ -f hcn/hcn.output ]; then
        if ! grep -q "Thank you very much for using Q-Chem." hcn/hcn.output; then
            echo "${i}/hcn"
        fi
    else
        echo "${i}/hcn"
    fi

    cd ../
done

    
