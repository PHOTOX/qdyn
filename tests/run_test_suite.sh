#!/bin/bash

tests="IT_HarmOsc_1D IT_HarmOsc_1D_project_rot IT_Morse_1D IT_HarmOsc_2D_sym IT_HarmOsc_2D_asym IT_HarmOsc_3D_sym IT_HarmOsc_3D_asym RT_HarmOsc_1D"

echo -e "\nRunning test suite"

for i in $tests;
do
  cd $i
  echo -e "\n**********************\n  $i\n**********************\n"
  ./test.sh
  cd ../
done

if [[ -e ERROR ]];
then
  echo 
  for i in $(find */ERROR) 
  do 
    echo -e "\033[0;31mERROR found: $i"
  done
  rm ERROR
else
  echo -e "\n\033[0;32mAll tests passed!"
fi

echo
