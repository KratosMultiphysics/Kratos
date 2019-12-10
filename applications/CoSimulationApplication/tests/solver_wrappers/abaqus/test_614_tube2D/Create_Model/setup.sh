#!/bin/sh 
module load ABAQUS/6.14

abaqus cae noGUI=makeInp.py

mv CSM_Time0.inp ../

rm *.rpy*
