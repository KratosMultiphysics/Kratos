#!/bin/bash

# README: run this script to reset and setup case

# clean working directory
if [ -d ./CFD ] 
then
    rm -rf ./CFD 
fi

mkdir ./CFD
cp mesh.jou CFD/
cp case.jou CFD/
cp udf_inlet.c CFD/
cd CFD

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh.jou > setup_gambit.log 2>&1
ml -GAMBIT

# make fluent case
ml ANSYS_CFD/2019R1
fluent 3ddp -g -i case.jou > setup_fluent.log 2>&1
ml -ANSYS_CFD

cd ..
