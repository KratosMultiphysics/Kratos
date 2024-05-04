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

#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "testing/testing.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix2D6NDiffOrderGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(model);

    Properties properties(0);
    // Please note these are not representative values, it just ensures the values are set
    properties.SetValue(DENSITY_WATER, 1000.0);
    properties.SetValue(POROSITY, 0.0);
    properties.SetValue(DENSITY_SOLID, 1700.0);

    auto        p_stress_state_policy = std::make_unique<PlaneStrainStressState>();
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
        make_shared<Triangle2D3<Node>>(r_geom(0), r_geom(1), r_geom(2));
    const auto& Np_container = p_pressure_geometry->ShapeFunctionsValues(integration_method);

    const auto solid_densities = GeoTransportEquationUtilities::CalculateSoilDensities(
        r_geom, number_of_integration_points, p_pressure_geometry->PointsNumber(), Np_container,
        p_retention_law, properties, process_info);
    Vector expected_solid_densities(number_of_integration_points);
    expected_solid_densities <<= 1700, 1700, 1700;

    KRATOS_CHECK_VECTOR_NEAR(solid_densities, expected_solid_densities, 1e-4)

    const auto integration_coefficients =
        GeoEquationOfMotionUtilities::CalculateIntegrationCoefficientInitialConfiguration(
            r_geom, integration_method, *p_stress_state_policy);
    Vector expected_integration_coefficients(number_of_integration_points);
    expected_integration_coefficients <<= 0.000416667, 0.000416667, 0.000416667;

    KRATOS_CHECK_VECTOR_NEAR(integration_coefficients, expected_integration_coefficients, 1e-4)

    const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
        r_geom, integration_method, solid_densities, integration_coefficients);

    Matrix expected_mass_matrix(r_geom.WorkingSpaceDimension() * r_geom.PointsNumber(),
                                r_geom.WorkingSpaceDimension() * r_geom.PointsNumber());
    // clang-format off
       expected_mass_matrix <<=
    0.0524691,0.0524691,-0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,
    0.0524691,0.0524691,-0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,
    -0.0262346,-0.0262346,0.0524691,0.0524691,-0.0262346,-0.0262346,0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,
    -0.0262346,-0.0262346,0.0524691,0.0524691,-0.0262346,-0.0262346,0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,
    -0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0524691,0.0524691,-0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,
    -0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0524691,0.0524691,-0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,
    0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.28858,0.28858,0.209877,0.209877,0.209877,0.209877,
    0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.28858,0.28858,0.209877,0.209877,0.209877,0.209877,
    -0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,0.209877,0.209877,0.28858,0.28858,0.209877,0.209877,
    -0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,0.209877,0.209877,0.28858,0.28858,0.209877,0.209877,
    0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,0.209877,0.209877,0.209877,0.209877,0.28858,0.28858,
    0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,0.209877,0.209877,0.209877,0.209877,0.28858,0.28858;
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

    const auto solid_densities = GeoTransportEquationUtilities::CalculateSoilDensities(
        r_geom, number_of_integration_points, r_geom.PointsNumber(), N_container, p_retention_law,
        properties, process_info);
    Vector expected_solid_densities(number_of_integration_points);
    expected_solid_densities <<= 2050, 2050, 2050, 2050;

    KRATOS_CHECK_VECTOR_NEAR(solid_densities, expected_solid_densities, 1e-4)

    const auto integration_coefficients =
        GeoEquationOfMotionUtilities::CalculateIntegrationCoefficientInitialConfiguration(
            r_geom, integration_method, *p_stress_state_policy);
    Vector expected_integration_coefficients(number_of_integration_points);
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