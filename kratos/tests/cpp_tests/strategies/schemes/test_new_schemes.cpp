//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"
#include "solving_strategies/schemes/new_scheme.h"

namespace Kratos::Testing
{

// typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
// typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

// typedef ModelPart::DofsArrayType DofsArrayType;

namespace
{

static void SetUpTestSchemesModelPart(ModelPart& rModelPart)
{
    const int domain_size = 3;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 0.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 0.0, 0.0, 0.0);

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
    }
}

}

KRATOS_TEST_CASE_IN_SUITE(NewScheme, KratosCoreFastSuite)
{
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    SetUpTestSchemesModelPart(r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_type" : "block"
    })");
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    auto p_scheme = Kratos::make_unique<NewScheme<SparseSpaceType, LocalSpaceType>>(scheme_settings);

    // Create the DOF set
    p_scheme->SetUpDofArray(r_test_model_part);

}

}  // namespace Kratos::Testing.

