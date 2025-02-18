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

#include "custom_constitutive/mohr_coulomb_constitutive_law.hpp"
#include "custom_constitutive/plane_strain.h"

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_Clone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto original_law = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());

    // Act
    auto p_cloned_law = original_law.Clone();

    // Assert
    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const MohrCoulombConstitutiveLaw*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto                        law = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());
    ConstitutiveLaw::Parameters parameters;
    Properties                  properties(3);
    parameters.SetMaterialProperties(properties);
    const auto element_geometry = Geometry<Node>{};
    const auto process_info     = ProcessInfo{};

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_COHESION is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_COHESION, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_FRICTION_ANGLE is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_FRICTION_ANGLE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_DILATION_ANGLE is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_DILATION_ANGLE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_TENSION_CUTOFF is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_TENSION_CUTOFF, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property: 3")
    properties.SetValue(YOUNG_MODULUS, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: POISSON_RATIO is not defined or has an invalid value for property: 3")
    properties.SetValue(POISSON_RATIO, 1.0);
    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_SetValueAndGetValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto       law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());
    const auto process_info = ProcessInfo{};

    // Act
    constexpr auto set_value = 1.0;
    law.SetValue(STATE_VARIABLE, set_value, process_info);

    // Assert
    auto zero_value = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(law.GetValue(STATE_VARIABLE, zero_value), set_value);
    zero_value = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(law.GetValue(DAMAGE_VARIABLE, zero_value), set_value);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_CalculatePK2Stress, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());

    Vector stress_vector = ZeroVector(plane_strain->GetStrainSize());
    Vector strain_vector = ScalarVector(plane_strain->GetStrainSize(), 1.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(strain_vector);
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    parameters.SetMaterialProperties(properties);

    // Act
    law.CalculatePK2Stress(strain_vector, stress_vector, parameters);

    // Assert
    Vector expected_stress_vector(plane_strain->GetStrainSize());
    expected_stress_vector <<=2.5e+07, 2.5e+07, 2.5e+07, 3.846155e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(stress_vector, expected_stress_vector, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_FindRegionIndex, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());

    constexpr auto fme = 1.0;
    constexpr auto fte = 1.0;

    // Act&Assert
    constexpr auto expected_result = 0;
    KRATOS_EXPECT_EQ(law.FindRegionIndex(fme, fte), expected_result);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_CalculateElasticMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    parameters.SetMaterialProperties(properties);

    // Act
    // Matrix elastic_matrix;
    Matrix elastic_matrix = ZeroMatrix(plane_strain->GetStrainSize(), plane_strain->GetStrainSize()); // temporarely

    law.CalculateElasticMatrix(elastic_matrix, parameters);

    // Assert
    Matrix expected_elastic_matrix(plane_strain->GetStrainSize(), plane_strain->GetStrainSize());
    //clang-format off
    expected_elastic_matrix <<=13461538.461538462,5769230.769230769,5769230.769230769,0,
                               5769230.769230769,13461538.461538462,5769230.769230769,0,
                               5769230.769230769,5769230.769230769,13461538.461538462,0,
                               0,0,0,3846153.8461538465;
    //clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(elastic_matrix, expected_elastic_matrix, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_CalculateMohrCoulomb, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 10.0);
    properties.SetValue(GEO_COHESION, 0.5);
    properties.SetValue(GEO_DILATION_ANGLE, 5.0);
    properties.SetValue(GEO_TENSION_CUTOFF, 0.5);
    parameters.SetMaterialProperties(properties);

    // Act
    Vector cauchy_stress_vector = ZeroVector(plane_strain->GetStrainSize());
    law.CalculateMohrCoulomb(properties, cauchy_stress_vector);

    // Assert
    Vector expected_cauchy_stress_vector(plane_strain->GetStrainSize());
    expected_cauchy_stress_vector <<= 0, 0, 0, 0 ;
    KRATOS_EXPECT_VECTOR_EQ(cauchy_stress_vector, expected_cauchy_stress_vector);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_FinalizeMaterialResponseCauchy, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    parameters.SetMaterialProperties(properties);
    Vector strain_vector = ScalarVector(plane_strain->GetStrainSize(), 1.0);
    parameters.SetStrainVector(strain_vector);

    // Act
    law.FinalizeMaterialResponseCauchy(parameters); // set mStrainVectorFinalized as strain_vector
    auto stress_vector = Vector(0.0, static_cast<int>(plane_strain->GetStrainSize()));
    law.CalculatePK2Stress(strain_vector, stress_vector, parameters); // use mStrainVectorFinalized and update mStressVector
    law.FinalizeMaterialResponseCauchy(parameters); // update mStrainVectorFinalized and mStressVectorFinalized
    law.CalculatePK2Stress(strain_vector, stress_vector, parameters); // use mStrainVectorFinalized and mStressVectorFinalized

    // Assert
    auto expected_stress_vector = ZeroVector(plane_strain->GetStrainSize());
    KRATOS_EXPECT_VECTOR_EQ(stress_vector, expected_stress_vector);
}

} // namespace Kratos::Testing
