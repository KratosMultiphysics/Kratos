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
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/mesh_moving_variables.h"
#include "utilities/time_discretization.h"

namespace Kratos {
namespace MeshVelocityCalculation {

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF1& rBDF);

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF2& rBDF);

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::Newmark& rGenAlpha);

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::Bossak& rGenAlpha);

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::GeneralizedAlpha& rGenAlpha);

} // namespace MeshVelocityCalculation
} // namespace Kratos.

#endif // KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED  defined
