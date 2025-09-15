//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED)
#define KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "utilities/time_discretization.h"

namespace Kratos {
class ModelPart; // forward-declaring to not having to include it here
namespace MeshVelocityCalculation {

void KRATOS_API(MESH_MOVING_APPLICATION) CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF1& rBDF);

void KRATOS_API(MESH_MOVING_APPLICATION) CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF2& rBDF);

void KRATOS_API(MESH_MOVING_APPLICATION) CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::Newmark& rGenAlpha);

void KRATOS_API(MESH_MOVING_APPLICATION) CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::Bossak& rGenAlpha);

void KRATOS_API(MESH_MOVING_APPLICATION) CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::GeneralizedAlpha& rGenAlpha);

} // namespace MeshVelocityCalculation
} // namespace Kratos.

#endif // KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED  defined
