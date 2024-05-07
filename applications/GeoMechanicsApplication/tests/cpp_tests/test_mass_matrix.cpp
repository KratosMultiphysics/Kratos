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

#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix2D6NDiffOrderGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Properties::Pointer p_properties(new Properties(0));
    // Please note these are not representative values, it just ensures the values are set
    p_properties->SetValue(DENSITY_WATER, 1000.0);
    p_properties->SetValue(POROSITY, 0.0);
    p_properties->SetValue(DENSITY_SOLID, 1700.0);
    // set arbitrary constitutive law
    const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic2DInterfaceLaw");
    p_properties->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    Model model;
    auto& r_model_part =
        ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(model, p_properties);

    auto&       r_element = r_model_part.GetElement(1);
    const auto& r_geom    = r_element.GetGeometry();

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

} // namespace Kratos::Testing