#!/bin/bash
xterm  -e mpirun -np 1 Emperor emperorInput.xml & EPID=$!
sleep 0.4s
xterm  -e mpirun -np 1 python subSystem1/subSystem1.py &
sleep 0.4s 
xterm  -e mpirun -np 1 python subSystem2/subSystem2.py &
wait $EPID
