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

int main(int argc, char **argv)
{

    try
    {
        if (argc != 6)
        {
            std::cerr << "Invalid arguments detected, usage: KratosGeoFlow.exe <working directory> <project parameters file> <minimum critical head> <maximum critical head> <critical head step size>" << std::endl;
            return -1;
        }

        string workingDirectory = argv[1];
        string projectName = argv[2];

        try
        {
            double minCriticalHead = stod(argv[3]);
            double maxCriticalHead = stod(argv[4]);
            double stepCriticalHead = stod(argv[5]);

            auto execute = Kratos::KratosExecute();
            execute.geoflow(workingDirectory, projectName, minCriticalHead, maxCriticalHead, stepCriticalHead);
        }
        catch (...)
        {
            std::cerr << "Could not parse critical head parameters to double values. Please check the input and try again.";
        }
    }
    catch (runtime_error e)
    {
        cout << "Runtime error: " << e.what();
    }
}
