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
#include "custom_constitutive/interface_plane_strain.h"
#include "custom_constitutive/interface_three_dimensional_surface.h"
#include "custom_geometries/interface_geometry.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include "custom_utilities/registration_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

inline GeoIncrementalLinearElasticInterfaceLaw MakeLaw2D()
{
    return GeoIncrementalLinearElasticInterfaceLaw{std::make_unique<InterfacePlaneStrain>()};
}

inline GeoIncrementalLinearElasticInterfaceLaw MakeLaw3D()
{
    return GeoIncrementalLinearElasticInterfaceLaw{std::make_unique<InterfaceThreeDimensionalSurface>()};
}

inline void InitializeLawMaterial(GeoIncrementalLinearElasticInterfaceLaw& law,
                                  const Properties& properties = Properties{})
{
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_geometry, dummy_shape_function_values);
}

inline void SetInitialState(GeoIncrementalLinearElasticInterfaceLaw& law,
                            const Vector&                            initial_relative_displacement,
                            const Vector&                            initial_traction)
{
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    law.SetInitialState(p_initial_state);
}

inline Properties MakeProperties(double normal_stiffness, double shear_stiffness)
{
    Properties properties;
    properties[INTERFACE_NORMAL_STIFFNESS] = normal_stiffness;
    properties[INTERFACE_SHEAR_STIFFNESS]  = shear_stiffness;
    return properties;
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_HasCorrectWorkingSpace, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = MakeLaw2D();
    KRATOS_EXPECT_EQ(law_2D.WorkingSpaceDimension(), 2);

    auto law_3D = MakeLaw3D();
    KRATOS_EXPECT_EQ(law_3D.WorkingSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_HasCorrectStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law_2D = MakeLaw2D();
    KRATOS_EXPECT_EQ(law_2D.GetStrainSize(), 2);

    const auto law_3D = MakeLaw3D();
    KRATOS_EXPECT_EQ(law_3D.GetStrainSize(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_UsesCauchyStressMeasure, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = MakeLaw2D();
    KRATOS_EXPECT_EQ(law_2D.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);

    auto law_3D = MakeLaw3D();
    KRATOS_EXPECT_EQ(law_3D.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_IsIncremental, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = MakeLaw2D();
    KRATOS_EXPECT_TRUE(law_2D.IsIncremental())
}

KRATOS_TEST_CASE_IN_SUITE(CloneOfLinearElasticLawForInterfaces_IsIndependentOfOriginal,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto original_law = MakeLaw2D();
    auto       p_cloned_law = original_law.Clone();

    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const GeoIncrementalLinearElasticInterfaceLaw*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CloneOfLinearElasticLawForInterfaces_HasCorrectStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto original_law_2D = MakeLaw2D();
    const auto p_cloned_law_2D = original_law_2D.Clone();
    KRATOS_EXPECT_EQ(p_cloned_law_2D->GetStrainSize(), 2);

    // 3D  case
    const auto original_law_3D = MakeLaw3D();
    const auto p_cloned_law_3D = original_law_3D.Clone();

    KRATOS_EXPECT_EQ(p_cloned_law_3D->GetStrainSize(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_DoesNotRequireInitializationOfMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = MakeLaw2D();
    KRATOS_EXPECT_FALSE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_ChecksForCorrectMaterialProperties,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = MakeLaw2D();
    auto       properties   = Properties{};
    const auto geometry     = InterfaceGeometry<Line2D2<Node>>{};
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

KRATOS_TEST_CASE_IN_SUITE(TheCalculatedConstitutiveMatrixIsADiagonalMatrixContainingNormalAndShearStiffnessValues,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto material_properties                        = Properties{};
    material_properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    material_properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;

    auto law_parameters = ConstitutiveLaw::Parameters{};
    law_parameters.SetMaterialProperties(material_properties);

    auto law_2D = MakeLaw2D();

    // 2D  case
    // Act
    auto actual_constitutive_matrix = Matrix{};
    law_2D.CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, actual_constitutive_matrix);

    // Assert
    auto expected_constitutive_matrix = Matrix{2, 2};
    // clang-format off
    expected_constitutive_matrix <<= 20.0, 0.0,
                                     0.0, 10.0;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_constitutive_matrix, expected_constitutive_matrix,
                                       Defaults::relative_tolerance)

    // 3D  case
    // Arrange
    auto law_3D = MakeLaw3D();

    // Act
    actual_constitutive_matrix.clear();
    law_3D.CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, actual_constitutive_matrix);

    // Assert
    expected_constitutive_matrix = Matrix{3, 3};
    // clang-format off
    expected_constitutive_matrix <<= 20.0, 0.0,  0.0,
                                     0.0, 10.0,  0.0,
                                     0.0,  0.0, 10.0;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_constitutive_matrix, expected_constitutive_matrix,
                                       Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TryingToCalculateTheValueOfAnUnsupportedMatrixVariableRaisesAnError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto        law                                = MakeLaw2D();
    const auto& r_some_unsupported_matrix_variable = ENGINEERING_STRAIN_TENSOR;
    auto        dummy_parameters                   = ConstitutiveLaw::Parameters{};
    auto        value                              = Matrix{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.CalculateValue(dummy_parameters, r_some_unsupported_matrix_variable, value),
        "Can't calculate value of ENGINEERING_STRAIN_TENSOR: unsupported variable")
}

KRATOS_TEST_CASE_IN_SUITE(WhenNoInitialStateIsGivenStartWithZeroRelativeDisplacementAndZeroTraction,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto dummy_properties = Properties{};

    // 2D  case
    auto law_2D = MakeLaw2D();
    InitializeLawMaterial(law_2D, dummy_properties);

    auto value = Vector{};
    law_2D.GetValue(STRAIN, value);
    const auto zero_vector_2D = Vector{ZeroVector{2}};
    KRATOS_EXPECT_VECTOR_NEAR(value, zero_vector_2D, Defaults::absolute_tolerance)
    law_2D.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_NEAR(value, zero_vector_2D, Defaults::absolute_tolerance)

    // 3D  case
    auto law_3D = MakeLaw3D();
    InitializeLawMaterial(law_3D, dummy_properties);

    value.clear();
    law_3D.GetValue(STRAIN, value);
    const auto zero_vector_3D = Vector{ZeroVector{3}};
    KRATOS_EXPECT_VECTOR_NEAR(value, zero_vector_3D, Defaults::absolute_tolerance)
    law_3D.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_NEAR(value, zero_vector_3D, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(WhenAnInitialStateIsGivenStartFromThereAfterMaterialInitialization,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // 2D  case
    const auto initial_relative_displacement = Vector{ScalarVector{2, 0.5}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    auto law_2D = MakeLaw2D();
    law_2D.SetInitialState(p_initial_state);
    const auto dummy_properties = Properties{};
    InitializeLawMaterial(law_2D, dummy_properties);

    auto value = Vector{};
    law_2D.GetValue(STRAIN, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, initial_relative_displacement, Defaults::relative_tolerance)
    law_2D.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, initial_traction, Defaults::relative_tolerance)

    // 3D  case
    const auto initial_relative_displacement_3D = Vector{ScalarVector{3, 0.5}};
    const auto initial_traction_3D              = Vector{ScalarVector{3, 30.0}};
    p_initial_state = make_intrusive<InitialState>(initial_relative_displacement_3D, initial_traction_3D);
    auto law_3D = MakeLaw3D();
    law_3D.SetInitialState(p_initial_state);
    InitializeLawMaterial(law_3D, dummy_properties);

    value.clear();
    law_3D.GetValue(STRAIN, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, initial_relative_displacement_3D, Defaults::relative_tolerance)
    law_3D.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, initial_traction_3D, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TryingToGetTheValueOfAnUnsupportedVectorVariableRaisesAnError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto        law                                = MakeLaw2D();
    const auto& r_some_unsupported_vector_variable = GREEN_LAGRANGE_STRAIN_VECTOR;

    auto value = Vector{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.GetValue(r_some_unsupported_vector_variable, value),
        "Can't get value of GREEN_LAGRANGE_STRAIN_VECTOR: unsupported variable")
}

KRATOS_TEST_CASE_IN_SUITE(ComputedIncrementalTractionIsProductOfIncrementalRelativeDisplacementAndStiffness,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arange
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;
    auto law_parameters                    = ConstitutiveLaw::Parameters{};

    // 2D  case
    auto relative_displacement = Vector{2};
    relative_displacement <<= 0.1, 0.3;
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{ZeroVector{2}};
    law_parameters.SetStressVector(traction);

    law_parameters.SetMaterialProperties(properties);
    auto law_2d = MakeLaw2D();
    InitializeLawMaterial(law_2d, properties);

    law_2d.CalculateMaterialResponseCauchy(law_parameters);

    auto expected_traction = Vector{2};
    expected_traction <<= 2.0, 3.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, Defaults::relative_tolerance)

    // 3D  case
    relative_displacement = Vector{3};
    relative_displacement <<= 0.1, 0.3, 0.5;
    law_parameters.SetStrainVector(relative_displacement);
    traction = Vector{ZeroVector{3}};
    law_parameters.SetMaterialProperties(properties);
    auto law_3D = MakeLaw3D();
    InitializeLawMaterial(law_3D, properties);

    law_3D.CalculateMaterialResponseCauchy(law_parameters);

    expected_traction = Vector{3};
    expected_traction <<= 2.0, 3.0, 5.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ComputedTractionIsSumOfPreviousTractionAndTractionIncrement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;

    auto law_parameters = ConstitutiveLaw::Parameters{};

    // 2D  case
    auto relative_displacement = Vector{ZeroVector{2}};
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{ZeroVector{2}};
    law_parameters.SetStressVector(traction);

    law_parameters.SetMaterialProperties(properties);
    const auto initial_relative_displacement = Vector{ScalarVector{2, 5.0}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    auto law_2D = MakeLaw2D();
    law_2D.SetInitialState(p_initial_state);
    InitializeLawMaterial(law_2D, properties);

    // First step
    relative_displacement <<= 5.1, 5.3;
    law_2D.CalculateMaterialResponseCauchy(law_parameters);
    law_2D.FinalizeMaterialResponseCauchy(law_parameters);

    // Second step
    relative_displacement <<= 5.2, 5.6;
    law_2D.CalculateMaterialResponseCauchy(law_parameters);
    law_2D.FinalizeMaterialResponseCauchy(law_parameters);

    auto expected_traction = Vector{2};
    expected_traction <<= 30.0 + (5.1 - 5.0) * 20.0 + (5.2 - 5.1) * 20.0,
        30.0 + (5.3 - 5.0) * 10.0 + (5.6 - 5.3) * 10.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, Defaults::relative_tolerance)
    auto expected_relative_displacement = Vector{2};
    expected_relative_displacement <<= 5.0 + (5.1 - 5.0) + (5.2 - 5.1), 5.0 + (5.3 - 5.0) + (5.6 - 5.3);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStrainVector(),
                                       expected_relative_displacement, Defaults::relative_tolerance)

    // 3D  case
    relative_displacement = Vector{ZeroVector{3}};
    law_parameters.SetStrainVector(relative_displacement);
    traction = Vector{ZeroVector{3}};
    law_parameters.SetMaterialProperties(properties);
    const auto initial_relative_displacement_3D = Vector{ScalarVector{3, 5.0}};
    const auto initial_traction_3D              = Vector{ScalarVector{3, 30.0}};
    auto       p_initial_state_3D =
        make_intrusive<InitialState>(initial_relative_displacement_3D, initial_traction_3D);
    auto law_3D = MakeLaw3D();
    law_3D.SetInitialState(p_initial_state_3D);
    InitializeLawMaterial(law_3D, properties);

    // First step
    relative_displacement <<= 5.1, 5.3, 5.5;
    law_3D.CalculateMaterialResponseCauchy(law_parameters);
    law_3D.FinalizeMaterialResponseCauchy(law_parameters);

    // Second step
    relative_displacement <<= 5.2, 5.6, 6.0;
    law_3D.CalculateMaterialResponseCauchy(law_parameters);
    law_3D.FinalizeMaterialResponseCauchy(law_parameters);

    expected_traction = Vector{3};
    expected_traction <<= 30.0 + (5.1 - 5.0) * 20.0 + (5.2 - 5.1) * 20.0,
        30.0 + (5.3 - 5.0) * 10.0 + (5.6 - 5.3) * 10.0, 30.0 + (5.5 - 5.0) * 10.0 + (6.0 - 5.5) * 10.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, Defaults::relative_tolerance)
    expected_relative_displacement = Vector{3};
    expected_relative_displacement <<= 5.0 + (5.1 - 5.0) + (5.2 - 5.1),
        5.0 + (5.3 - 5.0) + (5.6 - 5.3), 5.0 + (5.5 - 5.0) + (6.0 - 5.5);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStrainVector(),
                                       expected_relative_displacement, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesCanBeSavedToAndLoadedFromASerializer,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // 2D  case
    auto       law_2D                        = MakeLaw2D();
    const auto initial_relative_displacement = Vector{ScalarVector{2, 0.5}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    auto p_initial_state = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    law_2D.SetInitialState(p_initial_state);

    const auto dummy_properties = Properties{};
    InitializeLawMaterial(law_2D, dummy_properties);

    auto law_parameters        = ConstitutiveLaw::Parameters{};
    auto relative_displacement = Vector{2};
    relative_displacement <<= 0.1, 0.3;
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{2};
    traction <<= 20.0, 45.0;
    law_parameters.SetStressVector(traction);
    law_2D.FinalizeMaterialResponseCauchy(law_parameters);

    const auto scoped_registration_law =
        ScopedSerializerRegistration{"InterfacePlaneStrain"s, InterfacePlaneStrain{}};
    auto       serializer = Serializer{new std::stringstream{}};
    const auto tag_2D     = "test_2D_tag"s;
    serializer.save(tag_2D, law_2D);

    auto restored_law_2D = GeoIncrementalLinearElasticInterfaceLaw{nullptr};
    serializer.load(tag_2D, restored_law_2D);

    auto value = Vector{};
    restored_law_2D.GetValue(STRAIN, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, relative_displacement, Defaults::relative_tolerance)
    restored_law_2D.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, traction, Defaults::relative_tolerance)
    KRATOS_EXPECT_TRUE(restored_law_2D.HasInitialState())
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(restored_law_2D.GetInitialState().GetInitialStrainVector(),
                                       initial_relative_displacement, Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(restored_law_2D.GetInitialState().GetInitialStressVector(),
                                       initial_traction, Defaults::relative_tolerance)

    // 3D case
    auto       law_3D                           = MakeLaw3D();
    const auto initial_relative_displacement_3D = Vector{ScalarVector{3, 0.5}};
    const auto initial_traction_3D              = Vector{ScalarVector{3, 30.0}};
    auto p_initial_state_3D = make_intrusive<InitialState>(initial_relative_displacement, initial_traction);
    law_3D.SetInitialState(p_initial_state);
    InitializeLawMaterial(law_3D, dummy_properties);

    relative_displacement = Vector{3};
    relative_displacement <<= 0.1, 0.3, 0.5;
    law_parameters.SetStrainVector(relative_displacement);
    traction = Vector{3};
    traction <<= 20.0, 45.0, 45.0;
    law_parameters.SetStressVector(traction);
    law_3D.FinalizeMaterialResponseCauchy(law_parameters);

    const auto scoped_registration_law_3D = ScopedSerializerRegistration{
        "InterfaceThreeDimensionalSurface"s, InterfaceThreeDimensionalSurface{}};
    const auto tag_3D = "test_tag"s;
    serializer.save(tag_3D, law_3D);

    auto restored_law_3D = GeoIncrementalLinearElasticInterfaceLaw{nullptr};
    serializer.load(tag_3D, restored_law_3D);

    value.clear();
    restored_law_3D.GetValue(STRAIN, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, relative_displacement, Defaults::relative_tolerance)
    restored_law_3D.GetValue(CAUCHY_STRESS_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, traction, Defaults::relative_tolerance)
    KRATOS_EXPECT_TRUE(restored_law_3D.HasInitialState())
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(restored_law_3D.GetInitialState().GetInitialStrainVector(),
                                       initial_relative_displacement, Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(restored_law_3D.GetInitialState().GetInitialStressVector(),
                                       initial_traction, Defaults::relative_tolerance)
}

} // namespace Kratos::Testing
