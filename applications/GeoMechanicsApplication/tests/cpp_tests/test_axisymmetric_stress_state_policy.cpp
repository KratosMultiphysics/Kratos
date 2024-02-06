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
//                   Marjan Fathian
//

#include "containers/model.h"
#include "custom_elements/axisymmetric_stress_state_policy.h"
#include "custom_elements/stress_state_policy.h"
#include "geo_mechanics_application_constants.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateBMatrixWithValidGeometryReturnsCorrectResults, KratosGeoMechanicsFastSuite)
{
    std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressStatePolicy>();

    Model      model;
    ModelPart& model_part = model.CreateModelPart("Main");

    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    model_part.CreateNewElement("UPwSmallStrainElement2D3N", 1, node_ids, model_part.CreateNewProperties(0));

    Vector Np(3);
    Np <<= 1.0, 2.0, 3.0;

    // clang-format off
    Matrix GradNpT(3, 2);
    GradNpT <<= 1.0, 2.0,
                3.0, 4.0,
                5.0, 6.0;
    // clang-format on

    const Matrix calculated_matrix =
        p_stress_state_policy->CalculateBMatrix(GradNpT, Np, model_part.GetElement(1).GetGeometry());

    // clang-format off
    Matrix expected_matrix(4, 6);
    expected_matrix <<= 1   ,0  ,3   ,0 ,5   ,0,
                        0   ,2  ,0   ,4 ,0   ,6,
                        0.2 ,0  ,0.4 ,0 ,0.6 ,0,
                        2   ,1  ,4   ,3 ,6   ,5;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(calculated_matrix, expected_matrix, 1e-12)
}

} // namespace Kratos::Testing
