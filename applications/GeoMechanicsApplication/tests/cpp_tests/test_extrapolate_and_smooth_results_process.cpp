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
#include "processes/structured_mesh_generator_process.h"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"

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

ModelPart& CreateValidMeshedModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy_mesh", 2);

    // add variables to model part before meshing, these need to be there for the FinalizeSolutionStep of the UPw element
    result.AddNodalSolutionStepVariable(NODAL_AREA);
    result.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);

    result.AddNodalSolutionStepVariable(NODAL_DAMAGE_VARIABLE);


    result.AddNodalSolutionStepVariable(DISPLACEMENT);
    result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    result.AddNodalSolutionStepVariable(VELOCITY);
    result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    result.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);

    // Set up the test model part mesh
    // For some strange reason counterclockwise gives a negative area, so I used clockwise here
    auto p_point_1 = Kratos::make_intrusive<Node>(1, 1.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, 1.0, 2.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 3.0, 2.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 3.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({"number_of_divisions": 2,
                                              "element_name": "UPwSmallStrainElement2D3N",
                                              "condition_name": "LineCondition",
                                              "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(domain_geometry, result, mesher_parameters).Execute();

    auto p_linear_law            = std::make_shared<GeoLinearElasticPlaneStrain2DLaw>();
    auto& r_model_part_properties = result.GetProperties(0);
    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_linear_law);

    // Provide every element with a CAUCHY_STRESS_VECTOR for its integration points
    Vector one_ip_stress_vector(VOIGT_SIZE_2D_PLANE_STRAIN, 2.0);
    for (auto& rElement : result.Elements()) {
        rElement.Initialize(result.GetProcessInfo());
        std::vector<Vector> all_ip_stress_vectors;
        for (auto i = 0; i < rElement.GetGeometry().IntegrationPointsNumber(rElement.GetIntegrationMethod()); ++i) {
            all_ip_stress_vectors.push_back(one_ip_stress_vector);
        }
        rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, all_ip_stress_vectors, result.GetProcessInfo());
    }

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

    my_process.ExecuteInitializeSolutionStep();
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
    my_process.ExecuteInitializeSolutionStep();

    // The wrong value of 20 set in the model_part initialization should remain.
    const auto node1_area = model_part.GetNode(1).FastGetSolutionStepValue(NODAL_AREA);
    KRATOS_EXPECT_DOUBLE_EQ(node1_area, 20.);
    // The wrong value of unity matrix in the model_part initialization should remain.
    const auto node1_stress = model_part.GetNode(1).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
    KRATOS_EXPECT_DOUBLE_EQ(node1_stress(0,0), 1.);

}

KRATOS_TEST_CASE_IN_SUITE(CheckCauchyStressVectorPresentInIntegrationPoints, KratosGeoMechanicsFastSuite) {
    Model model;
    auto& model_part = CreateValidMeshedModelPart(model);
    Parameters settings;
    std::vector<Vector> all_ip_stress_vectors;

    // Get Cauchy stresses ( the Calculate does not calculate, only retrieve )
    model_part.GetElement(4).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, all_ip_stress_vectors, model_part.GetProcessInfo());

    KRATOS_EXPECT_DOUBLE_EQ(2.0, all_ip_stress_vectors[2][3]);
}

KRATOS_TEST_CASE_IN_SUITE(CheckSummedNodalAreaAndAveragedStressTensor, KratosGeoMechanicsFastSuite) {
    Model model;
    auto& model_part = CreateValidMeshedModelPart(model);

    model_part.GetProcessInfo().SetValue(NODAL_SMOOTHING, true);
    const Parameters settings;
    ExtrapolateAndSmoothResultsProcess my_process(model_part, settings);
    my_process.ExecuteInitializeSolutionStep();

    // loop over the elements calling the FinalizeSolutionStep
    for (auto& rElement : model_part.Elements()) {
        rElement.FinalizeSolutionStep(model_part.GetProcessInfo());
    }

    my_process.ExecuteBeforeOutputStep();

    KRATOS_EXPECT_DOUBLE_EQ(1.0, model_part.GetNode(1).FastGetSolutionStepValue(NODAL_AREA));
    KRATOS_EXPECT_DOUBLE_EQ(1.5, model_part.GetNode(2).FastGetSolutionStepValue(NODAL_AREA));
    KRATOS_EXPECT_DOUBLE_EQ(0.5, model_part.GetNode(3).FastGetSolutionStepValue(NODAL_AREA));
    KRATOS_EXPECT_DOUBLE_EQ(3.0, model_part.GetNode(5).FastGetSolutionStepValue(NODAL_AREA));

    KRATOS_EXPECT_NEAR(2.0, model_part.GetNode(1).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(0,0), 1.E-6);
    KRATOS_EXPECT_NEAR(2.0, model_part.GetNode(2).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(0,0), 1.E-6);
    KRATOS_EXPECT_NEAR(2.0, model_part.GetNode(3).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(0,0), 1.E-6);
    KRATOS_EXPECT_NEAR(2.0, model_part.GetNode(5).FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(0,0), 1.E-6);
}

}