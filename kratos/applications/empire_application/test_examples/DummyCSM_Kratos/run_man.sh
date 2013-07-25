#!/bin/bash
xterm  -hold -e mpirun -np 1 Emperor emperorInput.xml & EPID=$!
sleep 0.4s
xterm  -hold -e mpirun -np 1 dummyCSM dummyCSMInput.xml &
sleep 0.4s 
xterm  -hold -e mpirun -np 1 Kratos empireKratos.xml  &
wait $EPID
