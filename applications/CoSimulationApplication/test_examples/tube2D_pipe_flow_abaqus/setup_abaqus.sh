#!/bin/bash

# README: run this script to remove old data and setup case

export INTEL_LICENSE_FILE=28518@157.193.126.6
source /apps/SL6.3/Intel/compiler/2015.3.187/bin/compilervars.sh intel64
module load ABAQUS/6.14

# clean working directory
if [ -d ./CSM ]
then
    rm -rf ./CSM
fi

# create new folder
cp -r setup_abaqus CSM
cd CSM

# cp ../setup_abaqus/{Base.inp,makeHostFile.sh} /.
source makeHostFile.sh

cd ..