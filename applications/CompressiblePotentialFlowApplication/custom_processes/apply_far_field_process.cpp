//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez
//

// Project includes
#include "apply_far_field_process.h"

namespace Kratos {

// Constructor for ApplyFarFieldProcess Process
ApplyFarFieldProcess::ApplyFarFieldProcess(ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart)
{}

void ApplyFarFieldProcess::Execute()
{
    std::cout << "Hello World" << std::endl;
}
} // namespace Kratos.
