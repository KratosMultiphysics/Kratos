#!/bin/sh

rm -f CMakeCache.txt
rm -f *.cmake
rm -rf CMakeFiles\

cmake ../applications/CoSimulationApplication                                                \
-DCO_SIMULATION_APPLICATION=ON                                                               \
-DCO_SIMULATION_APPLICATION_PYTHON=ON                                                       
