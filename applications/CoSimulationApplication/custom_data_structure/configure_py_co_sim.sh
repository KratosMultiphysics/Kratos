#!/bin/sh

# this file is supposed to be called from "cmake_build"
# it configures and installs the python-only version of the CoSimulationApplication

## Instructions:
# Copy this file to "cmake_build" and run "sh configure_py_co_sim.sh"
# => no compilation is required
# => no paths have to be specified

rm -f CMakeCache.txt
rm -f *.cmake
rm -rf CMakeFiles\

cmake ../applications/CoSimulationApplication  \
-DCO_SIMULATION_APPLICATION_PYTHON=ON

make install