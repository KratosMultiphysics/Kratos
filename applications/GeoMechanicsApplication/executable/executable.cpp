// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#include "executable.h"
#include <iostream>

int main(int argc, char** argv) {

	//string workingDirectory = "C:/Development/Kratos/KratosStandaloneTest/test1";
    //string projectName = "ProjectParameters_stage1.json";

    string workingDirectory = argv[1];
    string projectName = argv[2];

    auto execute = Kratos::KratosExecute();
    execute.cpp_geomechanics(workingDirectory, projectName);
}
