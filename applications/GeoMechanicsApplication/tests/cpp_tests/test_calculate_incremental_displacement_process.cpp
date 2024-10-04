// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include "custom_processes/calculate_incremental_displacement_process.h"
#include "geo_mechanics_application_variables.h"
#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(CalculateIncrementalDisplacementProcess, KratosGeoMechanicsFastSuite)
{
    // set up model part
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(INCREMENTAL_DISPLACEMENT);

    auto p_node_1 = r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    std::vector<intrusive_ptr<Node>> nodes = {p_node_1, p_node_2};

    for (auto p_node : nodes) {
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
    }

    // set up displacement values
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT, 1) = Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};

    p_node_2->FastGetSolutionStepValue(DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{7.0, 12.0, 9.0};
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT, 1) = Kratos::array_1d<double, 3>{11.0, 8.0, 14.0};

    // initialize process to be tested
    auto process = CalculateIncrementalDisplacementProcess(r_model_part, {});
    process.Execute();

    const auto expected_incremental_displacement_1 = Kratos::array_1d<double, 3>{-3.0, -3.0, -3.0};
    const auto expected_incremental_displacement_2 = Kratos::array_1d<double, 3>{-4.0, 4.0, -5.0};

    KRATOS_CHECK_VECTOR_NEAR(p_node_1->FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT),
                             expected_incremental_displacement_1, 1e-12)
    KRATOS_CHECK_VECTOR_NEAR(p_node_2->FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT),
                             expected_incremental_displacement_2, 1e-12)
}

} // namespace Kratos::Testing
