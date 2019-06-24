#!/bin/sh

# this file is supposed to be called from "cmake_build"
# it configures and installs the python-only version of the CoSimulationApplication

rm -f CMakeCache.txt
rm -f *.cmake
rm -rf CMakeFiles\

cmake ../applications/CoSimulationApplication  \
-DCO_SIMULATION_APPLICATION_PYTHON=ON

make install