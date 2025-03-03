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

#include "custom_constitutive/coulomb_yield_surface.hpp"
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.hpp"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/stress_strain_utilities.h"

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include "utilities/math_utils.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_Clone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto original_law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    // Act
    auto p_cloned_law = original_law.Clone();

    // Assert
    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const MohrCoulombWithTensionCutOff*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto                        law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
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

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchy,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 10.0);
    properties.SetValue(GEO_COHESION, 0.5);
    properties.SetValue(GEO_DILATION_ANGLE, 5.0);
    properties.SetValue(GEO_TENSION_CUTOFF, 0.5);
    parameters.SetMaterialProperties(properties);

    // Act
    Vector cauchy_stress_vector = ZeroVector(plane_strain->GetStrainSize());
    cauchy_stress_vector <<= -10.0, -10.0, -10.0, 0.0;
    Vector strain_vector = ZeroVector(plane_strain->GetStrainSize());
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    Vector expected_cauchy_stress_vector(plane_strain->GetStrainSize());
    expected_cauchy_stress_vector <<= -10.0, -10.0, -10.0, 0.0;
    KRATOS_EXPECT_VECTOR_EQ(cauchy_stress_vector, expected_cauchy_stress_vector);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateRotationMatrix2D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<PlaneStrain>();
    auto law          = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    Vector stress_vector = ZeroVector(plane_strain->GetStrainSize());
    stress_vector <<= -6.0, -7.0, -14.0, 1.0;

    Vector eigenvalues_vector;
    Matrix eigenvectors_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(stress_vector, eigenvalues_vector, eigenvectors_matrix);

    Matrix rotation_matrix = law.CalculateRotationMatrix(eigenvectors_matrix);
    Matrix unit_matrix     = prod(rotation_matrix, trans(rotation_matrix));

    Matrix expected_unit_matrix = ZeroMatrix(3, 3);
    expected_unit_matrix <<= 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    KRATOS_EXPECT_MATRIX_NEAR(unit_matrix, expected_unit_matrix, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateRotationMatrix3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    auto plane_strain = std::make_unique<ThreeDimensional>();
    auto law          = MohrCoulombWithTensionCutOff(std::make_unique<ThreeDimensional>());

    Vector stress_vector = ZeroVector(plane_strain->GetStrainSize());
    stress_vector <<= 10.0, 50.0, 20.0, 40.0, 35.0, 45.0;

    Vector eigenvalues_vector;
    Matrix eigenvectors_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(stress_vector, eigenvalues_vector, eigenvectors_matrix);

    Matrix rotation_matrix = law.CalculateRotationMatrix(eigenvectors_matrix);
    Matrix unit_matrix     = prod(rotation_matrix, trans(rotation_matrix));

    Matrix expected_unit_matrix = ZeroMatrix(3, 3);
    expected_unit_matrix <<= 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    KRATOS_EXPECT_MATRIX_NEAR(unit_matrix, expected_unit_matrix, 1.0e-15);
}

} // namespace Kratos::Testing