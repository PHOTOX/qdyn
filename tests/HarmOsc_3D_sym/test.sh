#!/bin/bash

qdyn > qdyn.out

echo "#n      #exact  #calculated" > results.dat
nstates=$(ls energies.*.dat | wc -l)

for i in $(seq 1 $nstates); do
  file=energies.$i.dat
  tail -1 $file | awk '{print $2}' >> o
done

paste exact_results.dat o >> results.dat

rm wf* ener* o

cat results.dat
