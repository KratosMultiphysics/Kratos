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

void emptyLog(char *log) {}

bool emptyCancel()
{
    return false;
}

int main(int argc, char **argv)
{
    try
    {
        if (argc < 6)
        {
            std::cerr << "Invalid arguments detected, usage: KratosGeoFlow.exe <working directory> <project parameters file> <minimum critical head> <maximum critical head> <critical head step size> <optional: critical head boundary model part name>" << std::endl;
            return -1;
        }

        std::string workingDirectory = argv[1];
        std::string projectName = argv[2];

        std::string criticalHeadBoundaryModelPartName;
        if (argc == 7)
        {
            criticalHeadBoundaryModelPartName = argv[6];
        }

        try
        {
            double minCriticalHead = std::stod(argv[3]);
            double maxCriticalHead = std::stod(argv[4]);
            double stepCriticalHead = std::stod(argv[5]);

            auto execute = Kratos::KratosExecute();
            execute.execute_flow_analysis(workingDirectory,
                                          projectName,
                                          minCriticalHead,
                                          maxCriticalHead,
                                          stepCriticalHead,
                                          criticalHeadBoundaryModelPartName,
                                          &emptyLog,
                                          &emptyLog,
                                          &emptyCancel);
        }
        catch (...)
        {
            std::cerr << "Could not parse critical head parameters to double values. Please check the input and try again.";
        }
    }
    catch (std::runtime_error &e)
    {
        std::cout << "Runtime error: " << e.what();
    }
}
