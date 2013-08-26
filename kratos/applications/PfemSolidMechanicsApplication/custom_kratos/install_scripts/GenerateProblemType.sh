#!/bin/bash

echo -n "Create a new problem type (y [yes], n [no] ) [default= y]: "
read -e PROBLEM_CREATION

if [ -z $PROBLEM_CREATION ]; then
    PROBLEM_CREATION="y"
fi

PROGRAM_DIR=$HOME
PROGRAM_FOLDER=kratos


echo -n "Copy the new Script form python_scripts (y [yes], n [no] ) [default= y]: "
read -e PYTHON_SCRIPTS

if [ -z $PYTHON_SCRIPTS ]; then
    PYTHON_SCRIPTS="y"
fi


if [ $PYTHON_SCRIPTS == "y" ]; then

cp -u $PROGRAM_DIR/$PROGRAM_FOLDER/applications/PfemSolidApplication/python_scripts/script.py  $PROGRAM_DIR/$PROGRAM_FOLDER/applications/PfemSolidApplication/custom_problemtype/pfem_solid_application_scripts/

fi


if [ $PROBLEM_CREATION == "y" ]; then

cd $PROGRAM_DIR/$PROGRAM_FOLDER/applications/PfemSolidApplication/custom_problemtype

rm -R Kratos_Pfem_Solid_Application.gid

python $PROGRAM_DIR/$PROGRAM_FOLDER/problemtype_generator/problemtype.py Pfem_Solid_Application

echo ProblemType CREATED

else

echo ProblemType NOT CREATED

fi