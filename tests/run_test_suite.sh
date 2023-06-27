#!/bin/bash

tests="HarmOsc_1D Morse_1D HarmOsc_2D_sym HarmOsc_2D_asym HarmOsc_3D_sym HarmOsc_3D_asym"

echo -e "\nRunning test suite"

for i in $tests;
do
  cd $i
  echo -e "\n*********************\n  $i\n*********************\n"
  ./test.sh
  cd ../
done

if [[ ! -e */ERROR ]];
then
  echo -e "\n\033[0;32mAll tests passed!"
else
  echo 
  for i in $(find */ERROR) 
  do 
    echo -e "\033[0;31mERROR found: $i"
  done
fi

echo
