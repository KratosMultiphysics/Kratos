#!/bin/bash
xterm  -hold -e mpirun -np 1 Emperor emperorInput.xml & EPID=$!
sleep 0.4s
xterm  -hold -e mpirun -np 1 dummyCSM CSM/dummyCSMInput.xml &
sleep 0.4s 
xterm  -hold -e mpirun -np 1 python CFD/KratosOpenMP.py  &
wait $EPID
