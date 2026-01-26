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
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

namespace
{

Kratos::ProcessInfo GetProcessInfo() { return Kratos::ProcessInfo(); }

class MockLaw : public Kratos::ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MockLaw);

    Kratos::Matrix& CalculateValue(Parameters& rParameterValues,
                                   const Kratos::Variable<boost::numeric::ublas::matrix<double>>& rThisVariable,
                                   Kratos::Matrix& rValue) override
    {
        constexpr std::size_t voigt_size = 4;
        rValue                           = Kratos::ScalarMatrix(voigt_size, voigt_size, 2.0);
        return rValue;
    }

    SizeType GetStrainSize() const override { return 4; }
};

} // namespace

namespace Kratos::Testing
{
TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestsStiffnessContribution)
{
    constexpr std::size_t voigt_size = 4;
    constexpr std::size_t n          = 3;

    Geo::IntegrationCoefficientsGetter integration_coefficients_getter = []() {
        return std::vector{1.0, 2.0};
    };
    const auto               b                 = Matrix(voigt_size, n, 0.5);
    Geo::BMatricesGetter     b_matrices_getter = [b]() { return std::vector{b, b}; };
    Geo::StrainVectorsGetter strains_getter    = []() {
        return std::vector{Vector{2, 0.5}, Vector{2, 0.6}};
    };
    std::vector<ConstitutiveLaw::Pointer> mock_laws{std::make_shared<MockLaw>(), std::make_shared<MockLaw>()};
    Geo::ConstitutiveLawsGetter constitutive_laws_getter =
        [mock_laws]() -> const std::vector<ConstitutiveLaw::Pointer>& { return mock_laws; };
    Geo::PropertiesGetter properties_getter = []() -> const Properties& { return Properties{}; };

    constexpr std::size_t                               number_of_u_dof = 4;
    StiffnessCalculator<number_of_u_dof>::InputProvider provider(
        b_matrices_getter, strains_getter, integration_coefficients_getter, properties_getter,
        GetProcessInfo, constitutive_laws_getter);
    StiffnessCalculator<number_of_u_dof> calculator(provider);

    const auto actual_stiffness_matrix = calculator.LHSContribution().value();

    const auto expected_stiffness_matrix = ScalarMatrix(n, n, 24.0);
    KRATOS_CHECK_MATRIX_NEAR(actual_stiffness_matrix, expected_stiffness_matrix, 1e-4)
}
} // namespace Kratos::Testing