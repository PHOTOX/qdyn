#!/bin/bash

if [[ -z $(which python3) ]]
then
  echo -e "\npython3 command not evailable! Skipping!\n"
  exit
fi

if [[ -f ERROR ]]
then
  rm ERROR
fi

echo "Time:"
time qdyn > qdyn.out

echo

omega=$(python3 oscillation.py)

echo "# omega_ref   omega_calc  correct?" > results.dat

awk -v "omega=$omega" ' {if (sqrt(($1-omega)^2) > 1e-5) printf "%10.6f  %10.6f    no \n",$1,omega; else printf "%10.6f  %10.6f    yes \n",$1,omega}' exact_results.dat >> results.dat

rm wf* ener*

cat results.dat

if [[ ! -z $(grep no results.dat) ]]
then
  echo "ERROR: test failed!!"
  touch ERROR
  touch ../ERROR
fi
