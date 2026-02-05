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

    [[nodiscard]] SizeType GetStrainSize() const override { return mStressVector.size(); }

    Matrix mConstitutiveMatrix;
    Vector mStressVector;
};

struct StiffnessInputs {
    const Matrix BMatrix = UblasUtilities::CreateMatrix({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    const Vector Strain = Vector{4, 0.5};
    const Matrix ConstitutiveMatrix = UblasUtilities::CreateMatrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    const Vector                          StressVector = UblasUtilities::CreateVector({1, 2, 3});
    std::vector<ConstitutiveLaw::Pointer> MockLaws{
        std::make_shared<MockLaw>(ConstitutiveMatrix, StressVector),
        std::make_shared<MockLaw>(ConstitutiveMatrix, StressVector)};
    const Properties  MyProperties           = Properties{};
    const ProcessInfo MyProcessInfo          = ProcessInfo{};
    const double      IntegrationCoefficient = 0.5;
};

template <unsigned int NumberOfUDof>
typename StiffnessCalculator<NumberOfUDof>::InputProvider CreateStiffnessInputProvider(
    const StiffnessInputs& rStiffnessInputs, std::size_t NumberOfIntegrationPoints = 1)
{
    auto get_b_matrices = [&rStiffnessInputs, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rStiffnessInputs.BMatrix);
    };

    auto get_integration_coefficients = [&rStiffnessInputs, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rStiffnessInputs.IntegrationCoefficient);
    };
    auto get_strains = [&rStiffnessInputs, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rStiffnessInputs.Strain);
    };
    auto get_element_properties = [&rStiffnessInputs]() -> const auto& {
        return rStiffnessInputs.MyProperties;
    };
    auto get_process_info = [&rStiffnessInputs]() -> const auto& {
        return rStiffnessInputs.MyProcessInfo;
    };
    auto get_constitutive_laws = [&rStiffnessInputs]() -> const auto& {
        return rStiffnessInputs.MockLaws;
    };

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

    const auto stiffness_inputs = StiffnessInputs{};
    const auto provider =
        CreateStiffnessInputProvider<number_of_u_dof>(stiffness_inputs, number_of_integration_points);
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

    const auto stiffness_inputs = StiffnessInputs{};
    const auto provider =
        CreateStiffnessInputProvider<number_of_u_dof>(stiffness_inputs, number_of_integration_points);
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