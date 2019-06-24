#!/bin/sh

# this file is supposed to be called from "cmake_build"

rm -f CMakeCache.txt
rm -f *.cmake
rm -rf CMakeFiles\

cmake ../applications/CoSimulationApplication                                                \
-DCO_SIMULATION_APPLICATION_PYTHON=ON

make install