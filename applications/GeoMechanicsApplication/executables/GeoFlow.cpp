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

    try
    {
        string workingDirectory = argv[1];
        string projectName = argv[2];
        double minCriticalHead = stod(argv[3]);
        double maxCriticalHead = stod(argv[4]);
        double stepCriticalHead = stod(argv[5]);

        auto execute = Kratos::KratosExecute();
        execute.geoflow(workingDirectory, projectName, minCriticalHead, maxCriticalHead, stepCriticalHead);
    }
    catch (runtime_error e)
    {
        cout << "Runtime error: " << e.what();
    }

	
    

}
