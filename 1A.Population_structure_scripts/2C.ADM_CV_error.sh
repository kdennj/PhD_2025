#!/bin/bash -e


for ((K=2; K<=10; K+=1)); do
    echo Starting run for K=${K}
    for ((i=1; i<=10; i++)); do
        cd K${K}_log_files
        awk -v K=$K '$1=="CV"{print K, $4}' K${K}_R${i}.log >> all_CV_K${K}.txt
        cd ..
    done
done
