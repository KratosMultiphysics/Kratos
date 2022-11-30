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
#include "includes/kratos_application.h"
#include "includes/registry_auxiliaries.h"

// Registering operations
#include "operations/operation.h"

namespace Kratos
{

void KratosApplication::RegisterOperations()
{
    RegistryAuxiliaries::RegisterOperationWithPrototype("KratosMultiphysics", "Operation", Operation());
}

}