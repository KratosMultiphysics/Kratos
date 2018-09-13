#!/bin/bash

### THE USER SHOULD CHANGE THE FOLLOWING THREE LINES
export NEW_NAME_LOWERCASE=convection_diffusionr 
export NEW_NAME_MIXED=ConvectionDiffusionr
export NEW_NAME_UPPERCASE=CONVECTIONDIFFUSIONR 
### DO NOT CHANGE ANYTHING AFTER THIS POINT !!!

cd test_application
mv test_application.h ${NEW_NAME_LOWERCASE}_application.h
mv test_application.cpp ${NEW_NAME_LOWERCASE}_application.cpp
mv custom_python/test_python_application.cpp custom_python/${NEW_NAME_LOWERCASE}_python_application.cpp
mv TestApplication.py  ${NEW_NAME_MIXED}Application.py

##KRATOS_TEST_APPLICATION_H_INCLUDED    ----> KRATOS_NEWNAME_APPLICATION_H_INCLUDED
export MACRO_TEXT=KRATOS_${NEW_NAME_UPPERCASE}_APPLICATION_H_INCLUDED
sed -i "s/KRATOS_TEST_APPLICATION_H_INCLUDED/${MACRO_TEXT}/g" ${NEW_NAME_LOWERCASE}_application.h 
sed -i "s/KRATOS_TEST_APPLICATION_H_INCLUDED/${MACRO_TEXT}/g" ${NEW_NAME_LOWERCASE}_application.cpp 
sed -i "s/KRATOS_TEST_APPLICATION_H_INCLUDED/${MACRO_TEXT}/g" custom_python/${NEW_NAME_LOWERCASE}_python_application.cpp 

##KratosTestApplication    ----> KratosNewNameApplication
sed -i "s/KratosTestApplication/Kratos${NEW_NAME_MIXED}Application/g" ${NEW_NAME_LOWERCASE}_application.h 
sed -i "s/KratosTestApplication/Kratos${NEW_NAME_MIXED}Application/g" ${NEW_NAME_LOWERCASE}_application.cpp 
sed -i "s/KratosTestApplication/Kratos${NEW_NAME_MIXED}Application/g" custom_python/${NEW_NAME_LOWERCASE}_python_application.cpp 
sed -i "s/KratosTestApplication/Kratos${NEW_NAME_MIXED}Application/g" CMakeLists.txt
sed -i "s/KratosTestApplication/Kratos${NEW_NAME_MIXED}Application/g" ${NEW_NAME_MIXED}Application.py

##test_application    ----> new_name_application
sed -i "s/test_application/${NEW_NAME_LOWERCASE}_application/g" ${NEW_NAME_MIXED}Application.py


#include test_application.h    ----> #include newname_application.h
sed -i "s/test_application.h/${NEW_NAME_LOWERCASE}_application.h/g" ${NEW_NAME_LOWERCASE}_application.h 
sed -i "s/test_application.h/${NEW_NAME_LOWERCASE}_application.h/g" ${NEW_NAME_LOWERCASE}_application.cpp 
sed -i "s/test_application.h/${NEW_NAME_LOWERCASE}_application.h/g" custom_python/${NEW_NAME_LOWERCASE}_python_application.cpp 

#KRATOS_TEST_APPLICATION_SOURCES   ----> #KRATOS_NEWNAME_APPLICATION_SOURCES
sed -i "s/KRATOS_TEST_APPLICATION_SOURCES/KRATOS_${NEW_NAME_UPPERCASE}_APPLICATION_SOURCES/g" CMakeLists.txt

# # test_application.cpp    ----> new_name_application.cpp
# # test_python_application.cpp    ----> newname_python_application.cpp
sed -i "s/test_application.cpp/${NEW_NAME_LOWERCASE}_application.cpp/g" CMakeLists.txt 
sed -i "s/test_python_application.cpp/${NEW_NAME_LOWERCASE}_python_application.cpp/g" CMakeLists.txt

# # TestApplication.py    ----> NewNameApplication.py
sed -i "s/TestApplication.py/${NEW_NAME_MIXED}Application.py/g" CMakeLists.txt 


cd ..
mv test_application ${NEW_NAME_LOWERCASE}_application
