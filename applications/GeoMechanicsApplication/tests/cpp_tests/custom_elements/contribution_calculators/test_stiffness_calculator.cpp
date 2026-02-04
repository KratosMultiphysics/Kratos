// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_elements/contribution_calculators/stiffness_calculator.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/checks.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

namespace
{

using namespace Kratos;

class MockLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MockLaw);

    explicit MockLaw(const Matrix& rConstitutiveMatrix, const Vector& rStressVectors)
        : mConstitutiveMatrix(rConstitutiveMatrix), mStressVector(rStressVectors)
    {
    }

    Matrix& CalculateValue(Parameters&, const Variable<Matrix>&, Matrix& rValue) override
    {
        rValue = mConstitutiveMatrix;
        return mConstitutiveMatrix;
    }

    void CalculateMaterialResponseCauchy(Parameters& rValues) override
    {
        rValues.GetStressVector() = mStressVector;
    };

    SizeType GetStrainSize() const override { return mStressVector.size(); }

    Matrix mConstitutiveMatrix;
    Vector mStressVector;
};

template <unsigned int NumberOfUDof>
typename StiffnessCalculator<NumberOfUDof>::InputProvider CreateStiffnessInputProvider(
    const Matrix&                                rBMatrix,
    const Vector&                                rStrain,
    double                                       IntegrationCoefficient,
    const Properties&                            rProperties,
    const ProcessInfo&                           rProcessInfo,
    const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws,
    std::size_t                                  NumberOfIntegrationPoints = 1)
{
    auto get_b_matrices = [&rBMatrix, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rBMatrix);
    };
    auto get_integration_coefficients = [IntegrationCoefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, IntegrationCoefficient);
    };
    auto get_strains = [&rStrain, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rStrain);
    };
    auto get_element_properties = [&rProperties]() -> const auto& { return rProperties; };
    auto get_process_info       = [&rProcessInfo]() -> const auto& { return rProcessInfo; };
    auto get_constitutive_laws = [&rConstitutiveLaws]() -> const auto& { return rConstitutiveLaws; };

    return typename StiffnessCalculator<NumberOfUDof>::InputProvider(
        get_b_matrices, get_strains, get_integration_coefficients, get_element_properties,
        get_process_info, get_constitutive_laws);
}

} // namespace

namespace Kratos::Testing
{
TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestsStiffnessContribution)
{
    // Arrange
    constexpr std::size_t number_of_u_dof              = 4;
    constexpr auto        number_of_integration_points = std::size_t{2};

    const auto b_matrix = UblasUtilities::CreateMatrix({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    const auto strain = Vector{4, 0.5};
    const auto constitutive_matrix = UblasUtilities::CreateMatrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    const auto                            stress_vector = UblasUtilities::CreateVector({1, 2, 3});
    std::vector<ConstitutiveLaw::Pointer> mock_laws{
        std::make_shared<MockLaw>(constitutive_matrix, stress_vector),
        std::make_shared<MockLaw>(constitutive_matrix, stress_vector)};
    const auto properties              = Properties{};
    const auto process_info            = ProcessInfo{};
    const auto integration_coefficient = 0.5;

    StiffnessCalculator<number_of_u_dof>::InputProvider provider =
        CreateStiffnessInputProvider<number_of_u_dof>(b_matrix, strain, integration_coefficient, properties,
                                                      process_info, mock_laws, number_of_integration_points);
    StiffnessCalculator<number_of_u_dof> calculator(provider);

    // Act
    const auto actual_stiffness_matrix = calculator.LHSContribution().value();

    // Assert
    // The expected matrix is obtained by calculating B^T * C * B * weight for each
    // integration point and summing the results.
    auto expected_stiffness_matrix = UblasUtilities::CreateMatrix(
        {{1605, 1902, 2199, 2496}, {1854, 2196, 2538, 2880}, {2103, 2490, 2877, 3264}, {2352, 2784, 3216, 3648}});
    KRATOS_CHECK_MATRIX_NEAR(actual_stiffness_matrix, expected_stiffness_matrix, 1e-4)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestsStiffnessForceContribution)
{
    // Arrange
    constexpr std::size_t number_of_u_dof              = 4;
    constexpr auto        number_of_integration_points = std::size_t{2};

    const auto b_matrix = UblasUtilities::CreateMatrix({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    const auto strain = Vector{4, 0.5};
    const auto constitutive_matrix = UblasUtilities::CreateMatrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    const auto                            stress_vector = UblasUtilities::CreateVector({1, 2, 3});
    std::vector<ConstitutiveLaw::Pointer> mock_laws{
        std::make_shared<MockLaw>(constitutive_matrix, stress_vector),
        std::make_shared<MockLaw>(constitutive_matrix, stress_vector)};
    const auto properties              = Properties{};
    const auto process_info            = ProcessInfo{};
    const auto integration_coefficient = 0.5;

    StiffnessCalculator<number_of_u_dof>::InputProvider provider =
        CreateStiffnessInputProvider<number_of_u_dof>(b_matrix, strain, integration_coefficient, properties,
                                                      process_info, mock_laws, number_of_integration_points);
    StiffnessCalculator<number_of_u_dof> calculator(provider);

    // Act
    const auto actual_stiffness_force = calculator.RHSContribution();

    // Assert
    // The expected matrix is obtained by calculating - B^T * stress vector * weight for each
    // integration point and summing the results.
    auto expected_stiffness_force = UblasUtilities::CreateVector({-38, -44, -50, -56});
    KRATOS_CHECK_VECTOR_NEAR(actual_stiffness_force, expected_stiffness_force, 1e-4)
}

} // namespace Kratos::Testing