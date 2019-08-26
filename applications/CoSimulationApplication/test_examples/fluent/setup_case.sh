#!/bin/bash

# README: run this script to setup case

# clean working directory
if [ -d ./CFD ] 
then
    rm -rf ./CFD 
fi

mkdir ./CFD
cp mesh.jou CFD/
cp case.jou CFD/
cd CFD

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh.jou
ml -GAMBIT

# make fluent case
ml ANSYS_CFD/2019R1
fluent 2ddp -g -i case.jou
ml -ANSYS_CFD

cd ..
