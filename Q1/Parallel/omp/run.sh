#!/bin/bash
for i in 1 2 4 8
do
    export OMP_NUM_THREADS=$i
    ./mm | grep "^[0-9]\."
done