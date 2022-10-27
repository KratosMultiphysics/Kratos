//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Ruben Zorrilla
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/kernel.h"
#include "includes/registry.h"
#include "operations/operation.h"

// Registering operations
#include "operations/foo_operation.h"

namespace Kratos
{

void KratosApplication::RegisterOperations()
{
    KRATOS_REGISTER_OPERATION("KratosMultiphysics", "Operation", Operation())
}

}