//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/kratos_application.h"
#include "includes/registry.h"

// Registering processes
#include "processes/process.h"
#include "processes/output_process.h"

namespace Kratos
{

void KratosApplication::RegisterProcesses()
{
    KRATOS_REGISTER_PROCESS_WITH_PROTOTYPE("KratosMultiphysics", "Process", Process())
    KRATOS_REGISTER_PROCESS_WITH_PROTOTYPE("KratosMultiphysics", "OutputProcess", OutputProcess())
}

}