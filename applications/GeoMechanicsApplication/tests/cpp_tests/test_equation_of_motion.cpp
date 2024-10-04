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

#include "custom_constitutive/linear_elastic_2D_interface_law.h"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix2D6NDiffOrderGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
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
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

    ProcessInfo process_info;

    r_element.Initialize(process_info);

    Matrix mass_matrix;
    r_element.CalculateMassMatrix(mass_matrix, process_info);

    Matrix expected_mass_matrix(r_geom.WorkingSpaceDimension() * r_geom.PointsNumber() + 3,
                                r_geom.WorkingSpaceDimension() * r_geom.PointsNumber() + 3);
    // clang-format off
       expected_mass_matrix <<=
    0.0524691,0.0,-0.0262346,0.0,-0.0262346,0.0,0.0262346,0.0,-0.0524691,0.0,0.0262346,0.0,0,0,0,
    0.0,0.0524691,0.0,-0.0262346,0.0,-0.0262346,0.0,0.0262346,0.0,-0.0524691,0.0,0.0262346,0,0,0,
    -0.0262346,0.0,0.0524691,0.0,-0.0262346,0.0,0.0262346,0.0,0.0262346,0.0,-0.0524691,0.0,0,0,0,
    0.0,-0.0262346,0.0,0.0524691,0.0,-0.0262346,0.0,0.0262346,0.0,0.0262346,0.0,-0.0524691,0,0,0,
    -0.0262346,0.0,-0.0262346,0.0,0.0524691,0.0,-0.0524691,0.0,0.0262346,0.0,0.0262346,0.0,0,0,0,
    0.0,-0.0262346,0.0,-0.0262346,0.0,0.0524691,0.0,-0.0524691,0.0,0.0262346,0.0,0.0262346,0,0,0,
    0.0262346,0.0,0.0262346,0.0,-0.0524691,0.0,0.28858,0.0,0.209877,0.0,0.209877,0.0,0,0,0,
    0.0,0.0262346,0.0,0.0262346,0.0,-0.0524691,0.0,0.28858,0.0,0.209877,0.0,0.209877,0,0,0,
    -0.0524691,0.0,0.0262346,0.0,0.0262346,0.0,0.209877,0.0,0.28858,0.0,0.209877,0.0,0,0,0,
    0.0,-0.0524691,0.0,0.0262346,0.0,0.0262346,0.0,0.209877,0.0,0.28858,0.0,0.209877,0,0,0,
    0.0262346,0.0,-0.0524691,0.0,0.0262346,0.0,0.209877,0.0,0.209877,0.0,0.28858,0.0,0,0,0,
    0.0,0.0262346,0.0,-0.0524691,0.0,0.0262346,0.0,0.209877,0.0,0.209877,0.0,0.28858,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(mass_matrix, expected_mass_matrix, 1e-4)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMassMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      model;
    const auto nodal_variables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                               std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model, nodal_variables);

    // Set the element properties
    auto p_elem_prop = r_model_part.pGetProperties(0);
    // set arbitrary constitutive law

    p_elem_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());
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
    34.1667,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,0,0,0,0,
    0.0,34.1667,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0,0,0,0,
    0.0,0.0,34.1667,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,17.0833,0,0,0,0,
    17.0833,0.0,0.0,34.1667,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,0,0,0,0,
    0.0,17.0833,0.0,0.0,34.1667,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0,0,0,0,
    0.0,0.0,17.0833,0.0,0.0,34.1667,0.0,0.0,17.0833,0.0,0.0,17.0833,0,0,0,0,
    17.0833,0.0,0.0,17.0833,0.0,0.0,34.1667,0.0,0.0,17.0833,0.0,0.0,0,0,0,0,
    0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,34.1667,0.0,0.0,17.0833,0.0,0,0,0,0,
    0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,34.1667,0.0,0.0,17.0833,0,0,0,0,
    17.0833,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,34.1667,0.0,0.0,0,0,0,0,
    0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,34.1667,0.0,0,0,0,0,
    0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,17.0833,0.0,0.0,34.1667,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(mass_matrix, expected_mass_matrix, 1e-4)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDampingMatrixGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
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

KRATOS_TEST_CASE_IN_SUITE(CalculateStiffnessMatrixGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr std::size_t voigt_size = 4;
    constexpr std::size_t n          = 3;

    const Matrix        b = ScalarMatrix(voigt_size, n, 0.5);
    std::vector<Matrix> b_matrices;
    b_matrices.push_back(b);
    b_matrices.push_back(b);

    const Matrix        constitutive = ScalarMatrix(voigt_size, voigt_size, 2.0);
    std::vector<Matrix> constitutive_matrices;
    constitutive_matrices.push_back(constitutive);
    constitutive_matrices.push_back(constitutive);

    std::vector<double> integration_coefficients{1.0, 2.0};

    auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        b_matrices, constitutive_matrices, integration_coefficients);

    const auto expected_stiffness_matrix = ScalarMatrix(n, n, 24.0);

    KRATOS_CHECK_MATRIX_NEAR(stiffness_matrix, expected_stiffness_matrix, 1e-4)
}

KRATOS_TEST_CASE_IN_SUITE(TheInternalForceVectorIsTheIntegralOfBTransposedTimesSigma, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto b_matrix                 = Matrix{ScalarMatrix{2, 8, 1.0}};
    const auto b_matrices               = std::vector<Matrix>{b_matrix, b_matrix};
    const auto stress_vector            = Vector{ScalarVector{2, 1.0}};
    const auto stress_vectors           = std::vector<Vector>{stress_vector, stress_vector};
    const auto integration_coefficients = std::vector<double>{0.25, 0.4};

    const auto expected_internal_force_vector = Vector{ScalarVector{8, 1.3}};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(GeoEquationOfMotionUtilities::CalculateInternalForceVector(
                                           b_matrices, stress_vectors, integration_coefficients),
                                       expected_internal_force_vector, Defaults::relative_tolerance)
}

// The following tests only raise errors when using debug builds
#ifdef KRATOS_DEBUG

KRATOS_TEST_CASE_IN_SUITE(CalculatingTheInternalForceVectorFailsWhenTheInputVectorsHaveDifferentSizes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto b_matrix       = Matrix{ScalarMatrix{2, 8, 1.0}};
    const auto b_matrices     = std::vector<Matrix>{b_matrix}; // Error: missing one matrix
    const auto stress_vector  = Vector{ScalarVector{2, 1.0}};
    const auto stress_vectors = std::vector<Vector>{stress_vector, stress_vector};
    const auto integration_coefficients = std::vector<double>{0.25, 0.4};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoEquationOfMotionUtilities::CalculateInternalForceVector(b_matrices, stress_vectors, integration_coefficients),
        "Cannot calculate the internal force vector: input vectors have different sizes")
}

KRATOS_TEST_CASE_IN_SUITE(CalculatingTheInternalForceVectorFailsWhenAllInputVectorsAreEmpty,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto b_matrices               = std::vector<Matrix>{};
    const auto stress_vectors           = std::vector<Vector>{};
    const auto integration_coefficients = std::vector<double>{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoEquationOfMotionUtilities::CalculateInternalForceVector(b_matrices, stress_vectors, integration_coefficients),
        "Cannot calculate the internal force vector: input vectors are empty")
}

KRATOS_TEST_CASE_IN_SUITE(CalculatingTheInternalForceVectorFailsWhenBMatricesHaveDifferentSizes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto b_matrices = std::vector<Matrix>{ScalarMatrix{2, 8, 1.0}, ScalarMatrix{1, 8, 1.0}}; // Error: matrices have different numbers of rows
    const auto stress_vector            = Vector{ScalarVector{2, 1.0}};
    const auto stress_vectors           = std::vector<Vector>{stress_vector, stress_vector};
    const auto integration_coefficients = std::vector<double>{0.25, 0.4};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoEquationOfMotionUtilities::CalculateInternalForceVector(b_matrices, stress_vectors, integration_coefficients),
        "Cannot calculate the internal force vector: B-matrices have different sizes")
}

KRATOS_TEST_CASE_IN_SUITE(CalculatingTheInternalForceVectorFailsWhenStressVectorsHaveDifferentSizes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto b_matrix   = Matrix{ScalarMatrix{2, 8, 1.0}};
    const auto b_matrices = std::vector<Matrix>{b_matrix, b_matrix};
    const auto stress_vectors = std::vector<Vector>{ScalarVector{2, 1.0}, ScalarVector{3, 1.0}}; // Error: vectors have different sizes
    const auto integration_coefficients = std::vector<double>{0.25, 0.4};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoEquationOfMotionUtilities::CalculateInternalForceVector(b_matrices, stress_vectors, integration_coefficients),
        "Cannot calculate the internal force vector: stress vectors have different sizes")
}

KRATOS_TEST_CASE_IN_SUITE(CalculatingTheInternalForceVectorFailsWhenTheMatrixVectorProductCantBeComputed,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Error: transpose of the B-matrix has more columns (3) than the number of stress components (2)
    const auto b_matrix                 = Matrix{ScalarMatrix{3, 8, 1.0}};
    const auto b_matrices               = std::vector<Matrix>{b_matrix, b_matrix};
    const auto stress_vector            = Vector{ScalarVector{2, 1.0}};
    const auto stress_vectors           = std::vector<Vector>{stress_vector, stress_vector};
    const auto integration_coefficients = std::vector<double>{0.25, 0.4};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoEquationOfMotionUtilities::CalculateInternalForceVector(b_matrices, stress_vectors, integration_coefficients), "Cannot calculate the internal force vector: matrix-vector product cannot be calculated due to size mismatch")
}

#endif

} // namespace Kratos::Testing