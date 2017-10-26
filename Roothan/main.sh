#! /bin/bash

echo -e "\n--------------Welcome to the program--------------------\nChoose the basis state to start with: \n [1] Harmonic Basis \n [2] Square Well Basis States"

read -p "Enter your choice (1,2): " choice

case $choice in
1) g++ -std=c++14 -o main main-routine.cpp lib_harmonic.cpp
   ./main harmonic
;;

2) g++ -std=c++14 -o main main-routine.cpp lib_squarewell.cpp
   ./main squarewell
;;

*) echo "Please enter any of the numbers 1,2,3 as your choice"
;;

esac

echo -e "\n"
