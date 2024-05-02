// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "testing/testing.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix3D10NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(model);

    Properties properties(0);
    // Please note these are not representative values, it just ensures the values are set
    properties.SetValue(DENSITY_WATER, 1000.0);
    properties.SetValue(POROSITY, 0.3);
    properties.SetValue(DENSITY_SOLID, 2500.0);

    auto        p_stress_state_policy = std::make_unique<ThreeDimensionalStressState>();
    const auto  integration_method    = r_model_part.GetElement(1).GetIntegrationMethod();
    const auto& r_geom                = r_model_part.GetElement(1).GetGeometry();

    const Geometry<Node>::IntegrationPointsArrayType& integration_points =
        r_geom.IntegrationPoints(integration_method);
    std::vector<RetentionLaw::Pointer> p_retention_law(integration_points.size());
    for (unsigned int i = 0; i < p_retention_law.size(); ++i) {
        p_retention_law[i] = RetentionLawFactory::Clone(properties);
        p_retention_law[i]->InitializeMaterial(
            properties, r_geom, row(r_geom.ShapeFunctionsValues(integration_method), i));
    }

    ProcessInfo process_info;

    // execution
    const SizeType number_of_integration_points = r_geom.IntegrationPoints(integration_method).size();
    const Geometry<Node>::Pointer p_pressure_geometry =
        make_shared<Tetrahedra3D4<Node>>(r_geom(0), r_geom(1), r_geom(2), r_geom(3));
    const auto& Np_container = p_pressure_geometry->ShapeFunctionsValues(integration_method);

    const auto solid_densities = GeoTransportEquationUtilities::CalculateSoilDensityVector(
        r_geom, number_of_integration_points, Np_container, p_retention_law, properties, process_info);
    Vector expected_solid_densities(r_geom.PointsNumber() + p_pressure_geometry->PointsNumber());
    expected_solid_densities <<= 2050, 2050, 2050, 2050, 2050, 2050, 2050, 2050, 2050, 2050, 2050,
        2050, 2050, 2050;

    KRATOS_CHECK_VECTOR_NEAR(solid_densities, expected_solid_densities, 1e-4)

    const auto integration_coefficients = GeoEquationOfMotionUtilities::CalculateIntegrationCoefficientInitialConfiguration(
        r_geom, r_geom.IntegrationPoints(integration_method), integration_method, *p_stress_state_policy);
    Vector expected_integration_coefficients(r_geom.PointsNumber() + p_pressure_geometry->PointsNumber());
    expected_integration_coefficients <<= 0.0122488, 0.0122488, 0.0122488, 0.0122488, 0.0187813,
        0.0187813, 0.0187813, 0.0187813, 0.007091, 0.007091, 0.007091, 0.007091, 0.007091, 0.007091;

    KRATOS_CHECK_VECTOR_NEAR(integration_coefficients, expected_integration_coefficients, 1e-4)

    const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
        r_geom, integration_method, solid_densities, integration_coefficients);

    Matrix expected_mass_matrix(r_geom.WorkingSpaceDimension() * r_geom.PointsNumber(),
                                r_geom.WorkingSpaceDimension() * r_geom.PointsNumber());
    // clang-format off
    expected_mass_matrix <<=
4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,
4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,
4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,
0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,
0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,
0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,
0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,
0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,
0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,0.813492,0.813492,0.813492,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,
0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,
0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,
0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,0.813492,4.88095,4.88095,4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,
-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,
-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,
-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,
-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,
-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,
-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,
-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,
-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,
-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,
-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,
-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,
-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,
-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,
-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,
-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,13.0159,13.0159,13.0159,
-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,
-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317,
-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-4.88095,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,-3.25397,6.50794,6.50794,6.50794,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,13.0159,26.0317,26.0317,26.0317;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(mass_matrix_u, expected_mass_matrix, 1e-4)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Model      model;
    const auto nodal_variables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                               std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model, nodal_variables);

    Properties properties(0);
    // Please note these are not representative values, it just ensures the values are set
    properties.SetValue(DENSITY_WATER, 1000.0);
    properties.SetValue(POROSITY, 0.3);
    properties.SetValue(DENSITY_SOLID, 2500.0);

    auto        p_stress_state_policy = std::make_unique<ThreeDimensionalStressState>();
    const auto  integration_method    = r_model_part.GetElement(1).GetIntegrationMethod();
    const auto& r_geom                = r_model_part.GetElement(1).GetGeometry();

    const Geometry<Node>::IntegrationPointsArrayType& integration_points =
        r_geom.IntegrationPoints(integration_method);
    std::vector<RetentionLaw::Pointer> p_retention_law(integration_points.size());
    for (unsigned int i = 0; i < p_retention_law.size(); ++i) {
        p_retention_law[i] = RetentionLawFactory::Clone(properties);
        p_retention_law[i]->InitializeMaterial(
            properties, r_geom, row(r_geom.ShapeFunctionsValues(integration_method), i));
    }

    ProcessInfo process_info;

    // execution
    const SizeType number_of_integration_points = r_geom.IntegrationPoints(integration_method).size();
    const auto& N_container = r_geom.ShapeFunctionsValues(integration_method);

    const auto solid_densities = GeoTransportEquationUtilities::CalculateSoilDensityVector(
        r_geom, number_of_integration_points, N_container, p_retention_law, properties, process_info);
    Vector expected_solid_densities(r_geom.PointsNumber());
    expected_solid_densities <<= 2050, 2050, 2050, 2050;

    KRATOS_CHECK_VECTOR_NEAR(solid_densities, expected_solid_densities, 1e-4)

    const auto integration_coefficients = GeoEquationOfMotionUtilities::CalculateIntegrationCoefficientInitialConfiguration(
        r_geom, r_geom.IntegrationPoints(integration_method), integration_method, *p_stress_state_policy);
    Vector expected_integration_coefficients(r_geom.PointsNumber());
    expected_integration_coefficients <<= 0.0416667, 0.0416667, 0.0416667, 0.0416667;

    KRATOS_CHECK_VECTOR_NEAR(integration_coefficients, expected_integration_coefficients, 1e-4)

    const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
        r_geom, integration_method, solid_densities, integration_coefficients);

    Matrix expected_mass_matrix(r_geom.WorkingSpaceDimension() * r_geom.PointsNumber(),
                                r_geom.WorkingSpaceDimension() * r_geom.PointsNumber());
    // clang-format off
    expected_mass_matrix <<=
34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,
34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,
34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,
17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,
17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,
17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(mass_matrix_u, expected_mass_matrix, 1e-4)
}

} // namespace Kratos::Testing