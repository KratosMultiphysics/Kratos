//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Andreas Winterstein
//


// System includes


// External includes


// Project includes
#include "calculate_mesh_velocity_utility.h"


namespace Kratos
{

CalculateMeshVelocityUtility::CalculateMeshVelocityUtility(ModelPart& rModelPart,
                                                           Parameters Settings)
    : mrModelPart(rModelPart)
{
    // Check if ModelPart has MESH_DISPLACEMENT, MESH_VELOCITY & MESH_ACCELERATION

}

void CalculateMeshVelocityUtility::CalculateMeshVelocities()
{
    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
    KRATOS_ERROR_IF(delta_time < 1.0e-24) << "Detected delta_time = 0! "
        << "check if the time step is created correctly "
        << "for the current time step" << std::endl;

    // if method == "bdf"
    //     ....

}



}  // namespace Kratos.


