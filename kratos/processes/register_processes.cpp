//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "process.h"
#include "includes/kernel.h"

// registering processes
#include "output_process.h"

namespace Kratos
{
void KratosApplication::RegisterProcesses() {
    KRATOS_REGISTER_PROCESS("KratosMultiphysics", "OutputProcess", OutputProcess)
}
}