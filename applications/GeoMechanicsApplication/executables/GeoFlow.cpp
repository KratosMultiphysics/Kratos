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

#include "GeoFlow.h"
#include <iostream>

int main(int argc, char** argv) {

	//string workingDirectory = "C:/Development/Kratos/KratosStandaloneTest/test1";
    //string projectName = "ProjectParameters_stage1.json";

    string workingDirectory = argv[1];
    string projectName = argv[2];
    string hasPiping = argv[3];
    bool hasPipingFlag;
    hasPipingFlag = hasPiping=="1";
    cout <<"Test :"<< argv[3] << std::endl;
    auto execute = Kratos::KratosExecute();
    execute.geoflow(workingDirectory, projectName, hasPipingFlag);
}
