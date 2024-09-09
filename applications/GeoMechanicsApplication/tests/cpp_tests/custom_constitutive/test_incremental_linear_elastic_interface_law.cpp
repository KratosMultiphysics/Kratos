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
//                   Anne van de Graaf
//

#include "custom_constitutive/incremental_linear_elastic_interface_law.h"
#include "custom_geometries/line_interface_geometry.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace
{

constexpr auto absolute_tolerance = 1.0e-6;
constexpr auto relative_tolerance = 1.0e-6;

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesHas2DWorkingSpace, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesHasStrainSizeOfTwo, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_EQ(law.GetStrainSize(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesUsesCauchyStressMeasure, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_EQ(law.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesIsIncremental, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_TRUE(law.IsIncremental())
}

KRATOS_TEST_CASE_IN_SUITE(CloneOfLinearElasticLawForInterfacesIsIndependentOfOriginal, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto original_law = GeoIncrementalLinearElasticInterfaceLaw{};
    auto       p_cloned_law = original_law.Clone();

    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const GeoIncrementalLinearElasticInterfaceLaw*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesDoesNotRequireInitializationOfMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_FALSE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesChecksForCorrectMaterialProperties,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = GeoIncrementalLinearElasticInterfaceLaw{};
    auto       properties   = Properties{};
    const auto geometry     = LineInterfaceGeometry<Line2D2<Node>>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "No interface normal stiffness is defined")

    properties[INTERFACE_NORMAL_STIFFNESS] = -5.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface normal stiffness must be positive, but got -5")

    properties[INTERFACE_NORMAL_STIFFNESS] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface normal stiffness must be positive, but got 0")

    properties[INTERFACE_NORMAL_STIFFNESS] = 5.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "No interface shear stiffness is defined")

    properties[INTERFACE_SHEAR_STIFFNESS] = -2.5;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface shear stiffness must be positive, but got -2.5")

    properties[INTERFACE_SHEAR_STIFFNESS] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface shear stiffness must be positive, but got 0")

    properties[INTERFACE_SHEAR_STIFFNESS] = 2.5;
    KRATOS_EXPECT_EQ(law.Check(properties, geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(WhenNoInitialStateIsGivenStartWithZeroRelativeDisplacementAndZeroTraction,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    const auto dummy_properties            = Properties{};
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(dummy_properties, dummy_geometry, dummy_shape_function_values);

    auto value = Vector{};
    law.GetValue(STRAIN, value);
    const auto zero_vector = Vector{ZeroVector{2}};
    KRATOS_EXPECT_VECTOR_NEAR(value, zero_vector, absolute_tolerance)
    law.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_NEAR(value, zero_vector, absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(WhenAnInitialStateIsGivenStartFromThereAfterMaterialInitialization,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law                           = GeoIncrementalLinearElasticInterfaceLaw{};
    const auto initial_relative_displacement = Vector{ScalarVector{2, 0.5}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    law.SetInitialState(p_initial_state);

    const auto dummy_properties            = Properties{};
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(dummy_properties, dummy_geometry, dummy_shape_function_values);

    auto value = Vector{};
    law.GetValue(STRAIN, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, initial_relative_displacement, relative_tolerance)
    law.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, initial_traction, relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ComputedIncrementalTractionIsProductOfIncrementalRelativeDisplacementAndStiffness,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_parameters        = ConstitutiveLaw::Parameters{};
    auto relative_displacement = Vector{2};
    relative_displacement <<= 0.1, 0.3;
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{ZeroVector{2}};
    law_parameters.SetStressVector(traction);
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;
    law_parameters.SetMaterialProperties(properties);
    auto       law                         = GeoIncrementalLinearElasticInterfaceLaw{};
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_geometry, dummy_shape_function_values);

    law.CalculateMaterialResponseCauchy(law_parameters);

    auto expected_traction = Vector{2};
    expected_traction <<= 2.0, 3.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ComputedTractionIsSumOfPreviousTractionAndTractionIncrement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_parameters        = ConstitutiveLaw::Parameters{};
    auto relative_displacement = Vector{ZeroVector{2}};
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{ZeroVector{2}};
    law_parameters.SetStressVector(traction);
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;
    law_parameters.SetMaterialProperties(properties);
    auto       law                           = GeoIncrementalLinearElasticInterfaceLaw{};
    const auto initial_relative_displacement = Vector{ScalarVector{2, 5.0}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    law.SetInitialState(p_initial_state);
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_geometry, dummy_shape_function_values);

    // First step
    relative_displacement <<= 5.1, 5.3;
    law.CalculateMaterialResponseCauchy(law_parameters);
    law.FinalizeMaterialResponseCauchy(law_parameters);

    // Second step
    relative_displacement <<= 5.2, 5.6;
    law.CalculateMaterialResponseCauchy(law_parameters);
    law.FinalizeMaterialResponseCauchy(law_parameters);

    auto expected_traction = Vector{2};
    expected_traction <<= 30.0 + (5.1 - 5.0) * 20.0 + (5.2 - 5.1) * 20.0,
        30.0 + (5.3 - 5.0) * 10.0 + (5.6 - 5.3) * 10.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, relative_tolerance)
    auto expected_relative_displacement = Vector{2};
    expected_relative_displacement <<= 5.0 + (5.1 - 5.0) + (5.2 - 5.1), 5.0 + (5.3 - 5.0) + (5.6 - 5.3);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStrainVector(), expected_relative_displacement, relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesCanBeSavedToAndLoadedFromASerializer,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law                           = GeoIncrementalLinearElasticInterfaceLaw{};
    const auto initial_relative_displacement = Vector{ScalarVector{2, 0.5}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    law.SetInitialState(p_initial_state);

    const auto dummy_properties            = Properties{};
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(dummy_properties, dummy_geometry, dummy_shape_function_values);

    auto law_parameters        = ConstitutiveLaw::Parameters{};
    auto relative_displacement = Vector{2};
    relative_displacement <<= 0.1, 0.3;
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{2};
    traction <<= 20.0, 45.0;
    law_parameters.SetStressVector(traction);
    law.FinalizeMaterialResponseCauchy(law_parameters);

    auto       serializer = Serializer{new std::stringstream{}};
    const auto tag        = "test_tag"s;
    serializer.save(tag, law);

    auto restored_law = GeoIncrementalLinearElasticInterfaceLaw{};
    serializer.load(tag, restored_law);

    auto value = Vector{};
    restored_law.GetValue(STRAIN, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, relative_displacement, relative_tolerance)
    restored_law.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, traction, relative_tolerance)
    KRATOS_EXPECT_TRUE(restored_law.HasInitialState())
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(restored_law.GetInitialState().GetInitialStrainVector(),
                                       initial_relative_displacement, relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(restored_law.GetInitialState().GetInitialStressVector(),
                                       initial_traction, relative_tolerance)
}

} // namespace Kratos::Testing
