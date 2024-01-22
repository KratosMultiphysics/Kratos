// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Wijtze Pieter Kikstra
//

#include "testing/testing.h"
#include "custom_processes/extrapolate_and_smooth_results_process.hpp"
#include "containers/model.h"
#include "geo_mechanics_application_variables.h"

using namespace Kratos;
namespace Kratos::Testing {

namespace {

ModelPart& CreateValidModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy", 2);
    result.AddNodalSolutionStepVariable(NODAL_AREA);
    result.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    auto pNode = result.CreateNewNode(1, 0., 0., 0.);
    result.GetProcessInfo()[DOMAIN_SIZE] = 3;
    pNode->FastGetSolutionStepValue(NODAL_AREA) = 20.;
    Matrix unit_nodal_stress_tensor = identity_matrix<double>(result.GetProcessInfo()[DOMAIN_SIZE]);
    pNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = unit_nodal_stress_tensor;

    return result;
}

} // namespace


KRATOS_TEST_CASE_IN_SUITE(CreateExtrapolateAndSmoothProcess, KratosGeoMechanicsFastSuite) {
    Model model;
    auto& model_part = CreateValidModelPart(model);
    const Parameters settings;

    bool has_thrown = false;
    try{
        ExtrapolateAndSmoothResultsProcess my_process(model_part, settings);
    }
    catch(...){
        has_thrown = true;
    }
    KRATOS_EXPECT_FALSE(has_thrown)
}

KRATOS_TEST_CASE_IN_SUITE(CheckNodalAreaAndStressSetToZero, KratosGeoMechanicsFastSuite) {
    Model model;
    auto& model_part = CreateValidModelPart(model);
    model_part.GetProcessInfo().SetValue(NODAL_SMOOTHING, true);
    const Parameters settings;

    ExtrapolateAndSmoothResultsProcess my_process(model_part, settings);

    my_process.ExecuteBeforeOutputStep();
    const auto node1_area = model_part.GetNode(1).FastGetSolutionStepValue(NODAL_AREA);
    KRATOS_EXPECT_DOUBLE_EQ(node1_area, 0.);
    const auto node1_stress = model_part.GetNode(1).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
    KRATOS_EXPECT_DOUBLE_EQ(node1_stress(2,2),0.);
}

KRATOS_TEST_CASE_IN_SUITE(CheckProcesInfoNodalSmoothingNotSet, KratosGeoMechanicsFastSuite) {
    Model model;
    auto& model_part = CreateValidModelPart(model);
    Parameters settings;

    ExtrapolateAndSmoothResultsProcess my_process(model_part, settings);
    my_process.ExecuteBeforeOutputStep();

    // The wrong value of 20 set in the model_part initialization should remain.
    const auto node1_area = model_part.GetNode(1).FastGetSolutionStepValue(NODAL_AREA);
    KRATOS_EXPECT_DOUBLE_EQ(node1_area, 20.);
    // The wrong value of unity matrix in the model_part initialization should remain.
    const auto node1_stress = model_part.GetNode(1).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
    KRATOS_EXPECT_DOUBLE_EQ(node1_stress(0,0), 1.);

}

}