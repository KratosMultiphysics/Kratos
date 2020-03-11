#!/bin/bash

# README: run this script to setup case

# clean working directory
if [ -d ./CSM ]
then
    rm -rf ./CSM
fi
mkdir ./CSM

# copy setup files
cp setup_pipe_structure/solver_parameters.json CSM/