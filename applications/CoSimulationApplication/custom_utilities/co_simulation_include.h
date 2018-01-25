//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

#ifndef CO_SIMULATION_UTILITIES
#define CO_SIMULATION_UTILITIES

// system includes
#include "stdio.h"
#include <iostream>

std::string slash = "/";
std::string dot = ".";
std::string availExtension = "coav";

enum DataLocationOnMesh {
    ON_NODES,
    ON_ELEMENTS
};

bool CoSimulation_FileExists(const char *iFileName)
{
    return (access(iFileName, F_OK) != -1);
}

void CoSimulation_Wait(int iSeconds)
{
    sleep(iSeconds);
}

bool CoSimulation_IsDataFieldAvailable(const char *iFieldName)
{
    return CoSimulation_FileExists(iFieldName);
}

bool CoSimulation_MakeDataFieldAvailable(const char *iFieldName)
{
    std::ofstream outputFile(iFieldName);
    return outputFile.is_open();
}

bool CoSimulation_MakeDataFieldNotAvailable(const char *iFieldName)
{
    if (access( iFieldName, F_OK ) != -1) 
        return (remove(iFieldName) != 0);
    else 
        return 0;
}


#endif /* CO_SIMULATION_UTILITIES */