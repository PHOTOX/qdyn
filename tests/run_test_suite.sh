#!/bin/bash

# tests to be run
tests="IT_HarmOsc_1D IT_HarmOsc_1D_project_rot IT_Morse_1D IT_HarmOsc_2D_sym IT_HarmOsc_2D_asym IT_HarmOsc_3D_sym IT_HarmOsc_3D_asym RT_HarmOsc_1D"

echo -e "\nRunning test suite"

# checking if qdyn is in PATH and trying to add it if necessary
if [[ -z $(which qdyn) ]];
then
    echo -e "\nWARNING: qdyn is not set in your path. Trying to export qdyn to path from ../src."
    # checking if test are in the same folder as src with qdyn
    if [[ -z $(ls ../src/qdyn) ]]
    then
        echo "ERROR: ../src/qdyn not found. Exiting."
        exit
    else
        pwd=$(pwd | sed 's/tests/src/g')
        export PATH=$pwd/:$PATH
        echo "* qdyn exported to path: $(which qdyn)."
    fi
fi

# tests
for i in $tests;
do
  cd $i
  echo -e "\n**********************\n  $i\n**********************\n"
  ./test.sh
  cd ../
done

# looking for skipped tests due to lacking python or others
if [[ -e SKIP ]];
then
  echo
  for i in $(find */SKIP) 
  do 
    echo -e "\033[0;33mSKIP found: $i"
  done
  rm SKIP
fi

# looking for errors
if [[ -e ERROR ]];
then
  echo 
  for i in $(find */ERROR) 
  do 
    echo -e "\033[0;31mERROR found: $i"
  done
  rm ERROR
else
  echo -e "\n\033[0;32mAll tests passed!\n"
fi
