#!/bin/sh

cp ~/Kratos/cmake_build/*.sh ~/ 
cd
rm -rf Kratos
git clone https://github.com/KratosMultiphysics/Kratos/
cd Kratos
git checkout -b CoSimApp-development origin/CoSimApp-development
