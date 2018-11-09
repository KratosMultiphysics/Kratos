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
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/variable_utils.h"
#include "includes/mesh_moving_variables.h"

// Application includes
#include "custom_utilities/calculate_mesh_velocity_utility.h"

namespace Kratos {
namespace Testing {

void CreateModelPartForMeshVelComputation(
    Model& rModel,
    const std::size_t BufferSize,
    const bool RequiresMeshAcceleration)
{
    auto& r_model_part = rModel.CreateModelPart("MeshVelMP");

    Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.0, 0.0, 0.0);
    Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.0, 1.0, 0.0);
    Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 1.0, 1.0, 0.0);
    Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node<3> > geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 3,
        "element_name": "Element2D3N"
    })");

    r_model_part.SetBufferSize(BufferSize);

    r_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    if (RequiresMeshAcceleration) {
        r_model_part.AddNodalSolutionStepVariable(MESH_ACCELERATION);
    }

    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesBDF1, MeshMovingApplicationFastSuite)
{
    Model current_model;
    const std::size_t buffer_size = CalculateMeshVelocityUtility::GetMinimumBufferSize("bdf1");

    CreateModelPartForMeshVelComputation(current_model, buffer_size, false);

    auto& r_model_part = current_model.GetModelPart("MeshVelMP");
    Parameters utility_params(R"(
        {
            "time_scheme": "bdf1
        })");

    CalculateMeshVelocityUtility utility(r_model_part, utility_params);

}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesBDF2, MeshMovingApplicationFastSuite)
{
    Model current_model;
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesGeneralizedAlpha, MeshMovingApplicationFastSuite)
{
    Model current_model;
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesBossak, MeshMovingApplicationFastSuite)
{
    Model current_model;
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesNewmark, MeshMovingApplicationFastSuite)
{
    Model current_model;
}

}
}  // namespace Kratos.
