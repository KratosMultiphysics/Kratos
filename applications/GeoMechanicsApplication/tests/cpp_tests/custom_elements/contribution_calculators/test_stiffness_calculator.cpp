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

Kratos::ProcessInfo GetProcessInfo() { return Kratos::ProcessInfo(); }

using namespace Kratos;

class MockLaw : public Kratos::ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MockLaw);

    explicit MockLaw(const Matrix& rConstitutiveMatrix) : mConstitutiveMatrix(rConstitutiveMatrix)
    {
    }

    Kratos::Matrix& CalculateValue(Parameters& rParameterValues,
                                   const Kratos::Variable<boost::numeric::ublas::matrix<double>>& rThisVariable,
                                   Kratos::Matrix& rValue) override
    {
        rValue = mConstitutiveMatrix;
        return mConstitutiveMatrix;
    }

    void CalculateMaterialResponseCauchy(Parameters& rValues) override
    {
        rValues.GetStressVector() = UblasUtilities::CreateVector({1.0, 2.0, 3.0});
    };

    SizeType GetStrainSize() const override { return 3; }

    Matrix mConstitutiveMatrix;
};

} // namespace

namespace Kratos::Testing
{
TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestsStiffnessContribution)
{
    constexpr std::size_t number_of_u_dof = 4;

    Geo::IntegrationCoefficientsGetter integration_coefficients_getter = []() {
        return std::vector{1.0, 2.0};
    };
    const auto b = UblasUtilities::CreateMatrix({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    Geo::BMatricesGetter     b_matrices_getter = [b]() { return std::vector{b, b}; };
    Geo::StrainVectorsGetter strains_getter    = []() {
        return std::vector{Vector{4, 0.5}, Vector{4, 0.6}};
    };

    const auto constitutive_matrix = UblasUtilities::CreateMatrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    std::vector<ConstitutiveLaw::Pointer> mock_laws{std::make_shared<MockLaw>(constitutive_matrix),
                                                    std::make_shared<MockLaw>(constitutive_matrix)};
    Geo::ConstitutiveLawsGetter           constitutive_laws_getter =
        [mock_laws]() -> const std::vector<ConstitutiveLaw::Pointer>& { return mock_laws; };
    const auto            properties        = Properties{};
    Geo::PropertiesGetter properties_getter = [&properties]() -> const Properties& {
        return properties;
    };

    StiffnessCalculator<number_of_u_dof>::InputProvider provider(
        b_matrices_getter, strains_getter, integration_coefficients_getter, properties_getter,
        GetProcessInfo, constitutive_laws_getter);
    StiffnessCalculator<number_of_u_dof> calculator(provider);

    const auto actual_stiffness_matrix = calculator.LHSContribution().value();

    // The expected matrix is obtained by calculating B^T * C * B * weight for each
    // integration point and summing the results.
    auto expected_stiffness_matrix = UblasUtilities::CreateMatrix(
        {{1605, 1902, 2199, 2496}, {1854, 2196, 2538, 2880}, {2103, 2490, 2877, 3264}, {2352, 2784, 3216, 3648}});
    expected_stiffness_matrix*=3.0; // Sum of the integration weights
    KRATOS_CHECK_MATRIX_NEAR(actual_stiffness_matrix, expected_stiffness_matrix, 1e-4)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestsStiffnessForceContribution)
{
    constexpr std::size_t number_of_u_dof = 4;

    Geo::IntegrationCoefficientsGetter integration_coefficients_getter = []() {
        return std::vector{1.0, 2.0};
    };
    const auto b = UblasUtilities::CreateMatrix({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    Geo::BMatricesGetter     b_matrices_getter = [b]() { return std::vector{b, b}; };
    Geo::StrainVectorsGetter strains_getter    = []() {
        return std::vector{Vector{4, 0.5}, Vector{4, 0.6}};
    };

    const auto constitutive_matrix = UblasUtilities::CreateMatrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    std::vector<ConstitutiveLaw::Pointer> mock_laws{std::make_shared<MockLaw>(constitutive_matrix),
                                                    std::make_shared<MockLaw>(constitutive_matrix)};
    Geo::ConstitutiveLawsGetter           constitutive_laws_getter =
        [mock_laws]() -> const std::vector<ConstitutiveLaw::Pointer>& { return mock_laws; };
    const auto            properties        = Properties{};
    Geo::PropertiesGetter properties_getter = [&properties]() -> const Properties& {
        return properties;
    };

    StiffnessCalculator<number_of_u_dof>::InputProvider provider(
        b_matrices_getter, strains_getter, integration_coefficients_getter, properties_getter,
        GetProcessInfo, constitutive_laws_getter);
    StiffnessCalculator<number_of_u_dof> calculator(provider);

    const auto actual_stiffness_force = calculator.RHSContribution();

    // The expected matrix is obtained by calculating - B^T * stress vector * weight for each
    // integration point and summing the results.
    auto expected_stiffness_force = UblasUtilities::CreateVector({-38, -44, -50, -56});
    expected_stiffness_force*=3.0; // Sum of the integration weights
    KRATOS_CHECK_VECTOR_NEAR(actual_stiffness_force, expected_stiffness_force, 1e-4)
}

} // namespace Kratos::Testing