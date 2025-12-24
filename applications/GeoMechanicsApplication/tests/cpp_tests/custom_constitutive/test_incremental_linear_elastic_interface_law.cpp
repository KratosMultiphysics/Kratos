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
#include "custom_geometries/interface_geometry.hpp"
#include "custom_utilities/registration_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_2d_2.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace
{
GeoIncrementalLinearElasticInterfaceLaw CreateLaw2D()
{
    return GeoIncrementalLinearElasticInterfaceLaw{std::make_unique<InterfacePlaneStrain>()};
}

GeoIncrementalLinearElasticInterfaceLaw CreateLaw3D()
{
    return GeoIncrementalLinearElasticInterfaceLaw{std::make_unique<InterfaceThreeDimensionalSurface>()};
}

void InitializeLawMaterial(GeoIncrementalLinearElasticInterfaceLaw& rLaw,
                           const Properties&                        rProperties = Properties{})
{
    const auto dummy_geometry              = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    rLaw.InitializeMaterial(rProperties, dummy_geometry, dummy_shape_function_values);
}

void SetLawInitialState(GeoIncrementalLinearElasticInterfaceLaw& rLaw,
                        const Vector&                            rInitialRelativeDisplacement,
                        const Vector&                            rInitialTraction)
{
    const auto p_initial_state = make_intrusive<InitialState>(rInitialRelativeDisplacement, rInitialTraction);
    rLaw.SetInitialState(p_initial_state);
}

void TestInitialStates(GeoIncrementalLinearElasticInterfaceLaw& rLaw,
                       const Vector& rExpectedInitialRelativeDisplacement,
                       const Vector& rExpectedInitialTraction)
{
    KRATOS_EXPECT_TRUE(rLaw.HasInitialState())
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(rLaw.GetInitialState().GetInitialStrainVector(),
                                       rExpectedInitialRelativeDisplacement,
                                       Kratos::Testing::Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(rLaw.GetInitialState().GetInitialStressVector(),
                                       rExpectedInitialTraction, Kratos::Testing::Defaults::relative_tolerance)
}

void TestRelativeDisplacementVectorAndTractionVector(GeoIncrementalLinearElasticInterfaceLaw& rLaw,
                                                     const Vector& rExpectedRelativeDisplacementVector,
                                                     const Vector& rExpectedEffectiveTractionVector)
{
    auto value = Vector{};
    rLaw.GetValue(GEO_RELATIVE_DISPLACEMENT_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, rExpectedRelativeDisplacementVector,
                                       Kratos::Testing::Defaults::relative_tolerance)
    rLaw.GetValue(GEO_EFFECTIVE_TRACTION_VECTOR, value);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(value, rExpectedEffectiveTractionVector,
                                       Kratos::Testing::Defaults::relative_tolerance)
}
} // namespace

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_HasCorrectWorkingSpace, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = CreateLaw2D();
    KRATOS_EXPECT_EQ(law_2D.WorkingSpaceDimension(), 2);

    auto law_3D = CreateLaw3D();
    KRATOS_EXPECT_EQ(law_3D.WorkingSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_HasCorrectStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law_2D = CreateLaw2D();
    KRATOS_EXPECT_EQ(law_2D.GetStrainSize(), 2);

    const auto law_3D = CreateLaw3D();
    KRATOS_EXPECT_EQ(law_3D.GetStrainSize(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_UsesCauchyStressMeasure, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = CreateLaw2D();
    KRATOS_EXPECT_EQ(law_2D.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);

    auto law_3D = CreateLaw3D();
    KRATOS_EXPECT_EQ(law_3D.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_IsIncremental, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = CreateLaw2D();
    KRATOS_EXPECT_TRUE(law_2D.IsIncremental())

    auto law_3D = CreateLaw3D();
    KRATOS_EXPECT_TRUE(law_3D.IsIncremental())
}

KRATOS_TEST_CASE_IN_SUITE(CloneOfLinearElasticLawForInterfaces_IsIndependentOfOriginal,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto original_law = CreateLaw2D();
    auto       p_cloned_law = original_law.Clone();

    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const GeoIncrementalLinearElasticInterfaceLaw*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CloneOfLinearElasticLawForInterfaces_HasCorrectStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto original_law_2D = CreateLaw2D();
    const auto p_cloned_law_2D = original_law_2D.Clone();
    KRATOS_EXPECT_EQ(p_cloned_law_2D->GetStrainSize(), 2);

    // 3D  case
    const auto original_law_3D = CreateLaw3D();
    const auto p_cloned_law_3D = original_law_3D.Clone();

    KRATOS_EXPECT_EQ(p_cloned_law_3D->GetStrainSize(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_MovedLawHasCorrectStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       original_law_2D = CreateLaw2D();
    const auto moved_law_2D(std::move(original_law_2D));
    KRATOS_EXPECT_EQ(moved_law_2D.GetStrainSize(), 2);

    auto       original_law_3D = CreateLaw3D();
    const auto moved_law_3D    = std::move(original_law_3D);
    KRATOS_EXPECT_EQ(moved_law_3D.GetStrainSize(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_DoesNotRequireInitializationOfMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLaw2D();
    KRATOS_EXPECT_FALSE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfaces_ChecksForCorrectMaterialProperties,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = CreateLaw2D();
    auto       properties   = Properties{};
    const auto geometry     = InterfaceGeometry<Line2D2<Node>>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "INTERFACE_NORMAL_STIFFNESS does not exist in the material properties with Id 0.")

    properties[INTERFACE_NORMAL_STIFFNESS] = -5.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "INTERFACE_NORMAL_STIFFNESS in the material properties with Id 0 has "
        "an invalid value: -5 is out of the range (0, -).")

    properties[INTERFACE_NORMAL_STIFFNESS] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "INTERFACE_NORMAL_STIFFNESS in the material properties with Id 0 has "
        "an invalid value: 0 is out of the range (0, -).")

    properties[INTERFACE_NORMAL_STIFFNESS] = 5.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "INTERFACE_SHEAR_STIFFNESS does not exist in the material properties with Id 0.")

    properties[INTERFACE_SHEAR_STIFFNESS] = -2.5;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "INTERFACE_SHEAR_STIFFNESS in the material properties with Id 0 has "
        "an invalid value: -2.5 is out of the range (0, -).")

    properties[INTERFACE_SHEAR_STIFFNESS] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "INTERFACE_SHEAR_STIFFNESS in the material properties with Id 0 has "
        "an invalid value: 0 is out of the range (0, -).")

    properties[INTERFACE_SHEAR_STIFFNESS] = 2.5;
    KRATOS_EXPECT_EQ(law.Check(properties, geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(TheCalculatedConstitutiveMatrixIsADiagonalMatrixContainingNormalAndShearStiffnessValues,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;

    auto law_parameters = ConstitutiveLaw::Parameters{};
    law_parameters.SetMaterialProperties(properties);

    auto law_2D = CreateLaw2D();

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
    auto law_3D = CreateLaw3D();

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
    auto        law                                = CreateLaw2D();
    const auto& r_some_unsupported_matrix_variable = ENGINEERING_STRAIN_TENSOR;
    auto        dummy_parameters                   = ConstitutiveLaw::Parameters{};
    auto        value                              = Matrix{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.CalculateValue(dummy_parameters, r_some_unsupported_matrix_variable, value),
        "Can't calculate value of ENGINEERING_STRAIN_TENSOR: unsupported variable")
}

KRATOS_TEST_CASE_IN_SUITE(TryingToGetTheValueOfAnUnsupportedVectorVariableRaisesAnError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto        law                                = CreateLaw2D();
    const auto& r_some_unsupported_vector_variable = GREEN_LAGRANGE_STRAIN_VECTOR;

    auto value = Vector{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.GetValue(r_some_unsupported_vector_variable, value),
        "Can't get value of GREEN_LAGRANGE_STRAIN_VECTOR: unsupported variable")
}

KRATOS_TEST_CASE_IN_SUITE(WhenNoInitialStateIsGivenStartWithZeroRelativeDisplacementAndZeroTraction,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law_2D = CreateLaw2D();
    InitializeLawMaterial(law_2D);

    TestRelativeDisplacementVectorAndTractionVector(law_2D, Vector{ZeroVector{2}}, Vector{ZeroVector{2}});

    auto law_3D = CreateLaw3D();
    InitializeLawMaterial(law_3D);

    TestRelativeDisplacementVectorAndTractionVector(law_3D, Vector{ZeroVector{3}}, Vector{ZeroVector{3}});
}

KRATOS_TEST_CASE_IN_SUITE(WhenAnInitialStateIsGivenStartFromThereAfterMaterialInitialization,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law_2D                        = CreateLaw2D();
    const auto initial_relative_displacement = Vector{ScalarVector{2, 0.5}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    SetLawInitialState(law_2D, initial_relative_displacement, initial_traction);
    InitializeLawMaterial(law_2D);

    TestRelativeDisplacementVectorAndTractionVector(law_2D, initial_relative_displacement, initial_traction);

    auto       law_3D                           = CreateLaw3D();
    const auto initial_relative_displacement_3D = Vector{ScalarVector{3, 0.5}};
    const auto initial_traction_3D              = Vector{ScalarVector{3, 30.0}};
    SetLawInitialState(law_3D, initial_relative_displacement_3D, initial_traction_3D);
    InitializeLawMaterial(law_3D);

    TestRelativeDisplacementVectorAndTractionVector(law_3D, initial_relative_displacement_3D, initial_traction_3D);
}

KRATOS_TEST_CASE_IN_SUITE(ComputedIncrementalTractionIsProductOfIncrementalRelativeDisplacementAndStiffness,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arange
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 20.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 10.0;
    auto law_parameters                    = ConstitutiveLaw::Parameters{};

    auto relative_displacement = Vector{2};
    relative_displacement <<= 0.1, 0.3;
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{ZeroVector{2}};
    law_parameters.SetStressVector(traction);

    law_parameters.SetMaterialProperties(properties);
    auto law_2d = CreateLaw2D();
    InitializeLawMaterial(law_2d, properties);

    law_2d.CalculateMaterialResponseCauchy(law_parameters);

    auto expected_traction = Vector{2};
    expected_traction <<= 2.0, 3.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(law_parameters.GetStressVector(), expected_traction, Defaults::relative_tolerance)

    relative_displacement = Vector{3};
    relative_displacement <<= 0.1, 0.3, 0.5;
    law_parameters.SetStrainVector(relative_displacement);
    traction = Vector{ZeroVector{3}};
    law_parameters.SetMaterialProperties(properties);
    auto law_3D = CreateLaw3D();
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

    auto relative_displacement = Vector{ZeroVector{2}};
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{ZeroVector{2}};
    law_parameters.SetStressVector(traction);

    law_parameters.SetMaterialProperties(properties);
    auto law_2D = CreateLaw2D();
    SetLawInitialState(law_2D, Vector{ScalarVector{2, 5.0}}, Vector{ScalarVector{2, 30.0}});
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

    relative_displacement = Vector{ZeroVector{3}};
    law_parameters.SetStrainVector(relative_displacement);
    traction = Vector{ZeroVector{3}};
    law_parameters.SetMaterialProperties(properties);
    auto law_3D = CreateLaw3D();
    SetLawInitialState(law_3D, Vector{ScalarVector{3, 5.0}}, Vector{ScalarVector{3, 30.0}});
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
    auto       law_2D                        = CreateLaw2D();
    const auto initial_relative_displacement = Vector{ScalarVector{2, 0.5}};
    const auto initial_traction              = Vector{ScalarVector{2, 30.0}};
    SetLawInitialState(law_2D, initial_relative_displacement, initial_traction);
    InitializeLawMaterial(law_2D);

    auto law_parameters        = ConstitutiveLaw::Parameters{};
    auto relative_displacement = Vector{2};
    relative_displacement <<= 0.1, 0.3;
    law_parameters.SetStrainVector(relative_displacement);
    auto traction = Vector{2};
    traction <<= 20.0, 45.0;
    law_parameters.SetStressVector(traction);
    law_2D.FinalizeMaterialResponseCauchy(law_parameters);

    const auto scoped_registration_law =
        ScopedSerializerRegistration{std::make_pair("InterfacePlaneStrain"s, InterfacePlaneStrain{})};
    auto       serializer = Serializer{new std::stringstream{}};
    const auto tag_2D     = "test_2D_tag"s;
    serializer.save(tag_2D, law_2D);

    auto restored_law_2D = GeoIncrementalLinearElasticInterfaceLaw{nullptr};
    serializer.load(tag_2D, restored_law_2D);

    TestRelativeDisplacementVectorAndTractionVector(restored_law_2D, relative_displacement, traction);
    TestInitialStates(restored_law_2D, initial_relative_displacement, initial_traction);

    auto       law_3D                           = CreateLaw3D();
    const auto initial_relative_displacement_3D = Vector{ScalarVector{3, 0.5}};
    const auto initial_traction_3D              = Vector{ScalarVector{3, 30.0}};
    SetLawInitialState(law_3D, initial_relative_displacement_3D, initial_traction_3D);
    InitializeLawMaterial(law_3D);

    relative_displacement = Vector{3};
    relative_displacement <<= 0.1, 0.3, 0.5;
    law_parameters.SetStrainVector(relative_displacement);
    traction = Vector{3};
    traction <<= 20.0, 45.0, 45.0;
    law_parameters.SetStressVector(traction);
    law_3D.FinalizeMaterialResponseCauchy(law_parameters);

    const auto scoped_registration_law_3D = ScopedSerializerRegistration{
        std::make_pair("InterfaceThreeDimensionalSurface"s, InterfaceThreeDimensionalSurface{})};
    const auto tag_3D = "test_tag"s;
    serializer.save(tag_3D, law_3D);

    auto restored_law_3D = GeoIncrementalLinearElasticInterfaceLaw{nullptr};
    serializer.load(tag_3D, restored_law_3D);

    TestRelativeDisplacementVectorAndTractionVector(restored_law_3D, relative_displacement, traction);
    TestInitialStates(restored_law_3D, initial_relative_displacement_3D, initial_traction_3D);
}

} // namespace Kratos::Testing
