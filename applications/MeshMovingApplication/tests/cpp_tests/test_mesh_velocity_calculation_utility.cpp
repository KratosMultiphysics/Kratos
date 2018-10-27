//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/checks.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/variable_utils.h"

// Application includes
#include "mesh_moving_application.h"
#include "custom_utilities/calculate_mesh_velocity_utility.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesBDF, MeshMovingApplicationFastSuite)
{
    Model current_model;
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesGeneralizedAlpha, MeshMovingApplicationFastSuite)
{
    Model current_model;
}

}
}  // namespace Kratos.
