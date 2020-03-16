#!/bin/bash

# README: run this script to setup case

# clean working directory
if [ -d ./CFD ]
then
    rm -rf ./CFD
fi
mkdir ./CFD

# copy setup files
cp setup_tube_flow/solver_parameters.json CFD/