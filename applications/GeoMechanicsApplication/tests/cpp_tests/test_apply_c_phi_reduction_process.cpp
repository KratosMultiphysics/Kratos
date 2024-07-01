// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//
#include "containers/model.h"
#include "custom_processes/apply_c_phi_reduction_process.hpp"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "stub_linear_elastic_law.h"
#include "testing/testing.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace
{
ModelPart& PrepareCPhiTestModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy");

    // Set up the test model part mesh
    const auto             p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto             p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    const auto             p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    const auto             p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions": 2,
        "element_name": "Element2D3N",
        "condition_name": "LineCondition",
        "create_skin_sub_model_part": true
    })");

    StructuredMeshGeneratorProcess(domain_geometry, result, mesher_parameters).Execute();

    auto  p_dummy_law             = std::make_shared<Testing::StubLinearElasticLaw>();
    auto& r_model_part_properties = result.GetProperties(0);
    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_dummy_law);
    Vector umat_parameters(6);
    umat_parameters <<= 10000000, 0.2, 10.0, 25.0, 25.0, 1000;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    r_model_part_properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 3);
    r_model_part_properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 4);
    r_model_part_properties.SetValue(NUMBER_OF_UMAT_PARAMETERS, 6);

    return result;
}

void CheckReducedCPhi(const ModelPart& rModelPart, double COrig, double PhiOrig, double ReductionFactor)
{
    for (Element& rElement : rModelPart.Elements()) {
        const auto element_properties     = rElement.GetProperties();
        const auto umat_properties_vector = element_properties.GetValue(UMAT_PARAMETERS);
        const auto c_index   = element_properties.GetValue(INDEX_OF_UMAT_C_PARAMETER) - 1;
        const auto phi_index = element_properties.GetValue(INDEX_OF_UMAT_PHI_PARAMETER) - 1;

        KRATOS_EXPECT_DOUBLE_EQ(umat_properties_vector(c_index), ReductionFactor * COrig);

        const double phi_rad = MathUtils<>::DegreesToRadians(PhiOrig);
        const double tan_phi = std::tan(phi_rad);
        KRATOS_EXPECT_DOUBLE_EQ(std::tan(MathUtils<>::DegreesToRadians(umat_properties_vector(phi_index))),
                                ReductionFactor * tan_phi);
    };
}
} // namespace

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(CheckCAndPhiReducedAfterCallingApplyCPhiReductionProcess, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareCPhiTestModelPart(model);

    ApplyCPhiReductionProcess process{r_model_part, {}};
    process.ExecuteInitializeSolutionStep();

    CheckReducedCPhi(r_model_part, 10.0, 25.0, 0.9);
}

KRATOS_TEST_CASE_IN_SUITE(CheckCAndPhiTwiceReducedAfterCallingApplyCPhiReductionProcessTwice, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareCPhiTestModelPart(model);

    ApplyCPhiReductionProcess process{r_model_part, {}};
    process.ExecuteInitializeSolutionStep();
    process.ExecuteFinalizeSolutionStep();
    process.ExecuteInitializeSolutionStep();

    CheckReducedCPhi(r_model_part, 10.0, 25.0, 0.8);
}

} // namespace Kratos::Testing