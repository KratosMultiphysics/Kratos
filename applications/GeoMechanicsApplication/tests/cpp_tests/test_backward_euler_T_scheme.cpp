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
//

#include "custom_strategies/schemes/backward_euler_T_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

namespace Kratos::Testing {

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

KRATOS_TEST_CASE_IN_SUITE(BackwardEulerScheme_UpdatesVariablesDerivatives_WhenPredictIsCalled,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;
    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);

    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);

    const double current_temperature = 1.0;
    const double previous_temperature = 0.0;
    const double delta_time = 4.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->FastGetSolutionStepValue(TEMPERATURE) = current_temperature;
    p_node->FastGetSolutionStepValue(TEMPERATURE, 1) = previous_temperature;

    KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(DT_TEMPERATURE), 0.0);

    scheme.Initialize(model_part);
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;
    scheme.Predict(model_part, dof_set, A, Dx, b);

    const double expected_dt_temperature = 0.25;
    KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(DT_TEMPERATURE),
                            expected_dt_temperature);
}

} // namespace Kratos::Testing