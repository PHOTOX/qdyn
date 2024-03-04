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

# checking numpy
if [[ -z $(pip3 list | grep numpy) ]]
then
  echo -e "\nERROR: numpy package in python 3 not found.\nAvailable packagaes are: \n$(pip3 list)"
  touch ../SKIP
  touch SKIP
  exit
fi

echo "Time:"
time qdyn > qdyn.out

echo

test_result=$(python3 check_double_well.py)

if [[ $test_result == "True" ]]
then
  echo "Test passed - adiabatic populations match the reference." > results.dat
elif [[ $test_result == "False" ]]
then
  echo "Test failed - adiabatic populations don't match the reference and energy is conserved." > results.dat
else
  echo "Test failed - error in python script, it prints '$test_result' instead of True or False." > results.dat
fi

rm wf* ener*

cat results.dat

if [[ ! -z $(grep "Test failed" results.dat) ]]
then
  echo "ERROR: test failed!!"
  touch ERROR
  touch ../ERROR
fi
