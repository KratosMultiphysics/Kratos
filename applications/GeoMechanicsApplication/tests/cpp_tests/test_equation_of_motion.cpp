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

#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
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

    auto&       r_element   = r_model_part.GetElement(1);
    const auto& r_geom      = r_element.GetGeometry();
    auto        p_elem_prop = r_model_part.pGetProperties(0);
    p_elem_prop->SetValue(DENSITY_WATER, 1000.0);
    p_elem_prop->SetValue(POROSITY, 0.0);
    p_elem_prop->SetValue(DENSITY_SOLID, 1700.0);
    // set arbitrary constitutive law
    const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic2DInterfaceLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    ProcessInfo process_info;

    r_element.Initialize(process_info);

    Matrix mass_matrix;
    r_element.CalculateMassMatrix(mass_matrix, process_info);

    Matrix expected_mass_matrix(r_geom.WorkingSpaceDimension() * r_geom.PointsNumber() + 3,
                                r_geom.WorkingSpaceDimension() * r_geom.PointsNumber() + 3);
    // clang-format off
       expected_mass_matrix <<=
    0.0524691,0.0524691,-0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,0,0,0,
    0.0524691,0.0524691,-0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,0,0,0,
    -0.0262346,-0.0262346,0.0524691,0.0524691,-0.0262346,-0.0262346,0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0,0,0,
    -0.0262346,-0.0262346,0.0524691,0.0524691,-0.0262346,-0.0262346,0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0,0,0,
    -0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0524691,0.0524691,-0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,0,0,0,
    -0.0262346,-0.0262346,-0.0262346,-0.0262346,0.0524691,0.0524691,-0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,0,0,0,
    0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.28858,0.28858,0.209877,0.209877,0.209877,0.209877,0,0,0,
    0.0262346,0.0262346,0.0262346,0.0262346,-0.0524691,-0.0524691,0.28858,0.28858,0.209877,0.209877,0.209877,0.209877,0,0,0,
    -0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,0.209877,0.209877,0.28858,0.28858,0.209877,0.209877,0,0,0,
    -0.0524691,-0.0524691,0.0262346,0.0262346,0.0262346,0.0262346,0.209877,0.209877,0.28858,0.28858,0.209877,0.209877,0,0,0,
    0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,0.209877,0.209877,0.209877,0.209877,0.28858,0.28858,0,0,0,
    0.0262346,0.0262346,-0.0524691,-0.0524691,0.0262346,0.0262346,0.209877,0.209877,0.209877,0.209877,0.28858,0.28858,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(mass_matrix, expected_mass_matrix, 1e-4)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Model      model;
    const auto nodal_variables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                               std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model, nodal_variables);

    // Set the element properties
    auto p_elem_prop = r_model_part.pGetProperties(0);
    // set arbitrary constitutive law
    const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic2DInterfaceLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());
    // Please note these are not representative values, it just ensures the values are set
    p_elem_prop->SetValue(DENSITY_WATER, 1000.0);
    p_elem_prop->SetValue(POROSITY, 0.3);
    p_elem_prop->SetValue(DENSITY_SOLID, 2500.0);

    ProcessInfo process_info;

    auto&       r_element = r_model_part.GetElement(1);
    const auto& r_geom    = r_element.GetGeometry();
    r_element.Initialize(process_info);

    Matrix mass_matrix;
    r_element.CalculateMassMatrix(mass_matrix, process_info);

    Matrix expected_mass_matrix((r_geom.WorkingSpaceDimension() + 1) * r_geom.PointsNumber(),
                                (r_geom.WorkingSpaceDimension() + 1) * r_geom.PointsNumber());
    // clang-format off
    expected_mass_matrix <<=
    34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,0,0,0,0,
    34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,0,0,0,0,
    34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,17.0833,17.0833,17.0833,0,0,0,0,
    17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,0,0,0,0,
    17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,0,0,0,0,
    17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,17.0833,34.1667,34.1667,34.1667,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(mass_matrix, expected_mass_matrix, 1e-4)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDampingMatrixGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    constexpr std::size_t n = 10;

    constexpr double mass_matrix_value = 10;
    const auto       mass_matrix       = scalar_matrix(n, n, mass_matrix_value);

    constexpr double stiffness_matrix_value = 20;
    const auto       stiffness_matrix       = scalar_matrix(n, n, stiffness_matrix_value);

    double rayleigh_alpha = 0.0;
    double rayleigh_beta  = 1.0;
    auto   damping_matrix = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        rayleigh_alpha, rayleigh_beta, mass_matrix, stiffness_matrix);

    auto expected_damping_matrix = scalar_matrix(n, n, stiffness_matrix_value);

    KRATOS_CHECK_MATRIX_NEAR(damping_matrix, expected_damping_matrix, 1e-4)

    rayleigh_alpha = 1.0;
    rayleigh_beta  = 0.0;
    damping_matrix = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        rayleigh_alpha, rayleigh_beta, mass_matrix, stiffness_matrix);

    expected_damping_matrix = scalar_matrix(n, n, mass_matrix_value);

    KRATOS_CHECK_MATRIX_NEAR(damping_matrix, expected_damping_matrix, 1e-4)

    rayleigh_alpha = 0.5;
    rayleigh_beta  = 0.5;
    damping_matrix = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        rayleigh_alpha, rayleigh_beta, mass_matrix, stiffness_matrix);

    const double expected_matrix_value = rayleigh_alpha * mass_matrix_value + rayleigh_beta * stiffness_matrix_value;
    expected_damping_matrix = scalar_matrix(n, n, expected_matrix_value);

    KRATOS_CHECK_MATRIX_NEAR(damping_matrix, expected_damping_matrix, 1e-4)
}

} // namespace Kratos::Testing