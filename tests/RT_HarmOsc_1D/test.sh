#!/bin/bash

# removing SKIP and ERROR files
if [[ -f ERROR ]]
then
  rm ERROR
fi

if [[ -f SKIP ]]
then
  rm SKIP
fi

# checking python3
if [[ -z $(which python3) ]]
then
  echo -e "\nERROR: python3 command not evailable! Skipping!\n"
  touch ../SKIP
  touch SKIP
  exit
fi

# checking scipy
if [[ -z $(pip3 list | grep scipy) ]]
then
  echo -e "\nERROR: scipy package in python 3 not found.\nAvailable packagaes are: \n$(pip3 list)"
  touch ../SKIP
  touch SKIP
  exit
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
