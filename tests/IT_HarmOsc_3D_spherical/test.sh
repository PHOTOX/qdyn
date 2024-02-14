#!/bin/bash

if [[ -f ERROR ]]
then
  rm ERROR
fi

echo "Time:"
time qdyn > qdyn.out

echo

echo "# n     exact         calculated  E_diff      correct?" > results.dat
echo "# ----------------------------------------------------" >> results.dat
nstates=$(ls energies.*.dat | wc -l)

for i in $(seq 1 $nstates); do
  file=energies.$i.dat
  tail -1 $file | awk '{print $2 "  "}' >> o
done

paste exact_results.dat o | awk ' {if (sqrt(($3-$2)^2) > 1e-5) printf "  %d %12.6f  %12.6f  %10.6f    no \n",$1,$2,$3,$3-$2; else printf "  %d %12.6f  %12.6f  %10.6f    yes \n",$1,$2,$3,$3-$2}' >> results.dat

rm wf* ener* o 

cat results.dat

if [[ ! -z $(grep no results.dat) ]]
then
  echo "ERROR: test failed!!"
  touch ERROR
  touch ../ERROR
fi
