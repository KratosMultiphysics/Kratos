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

#pragma once

#include "contribution_calculator.h"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_aliases.h"
#include "includes/constitutive_law.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

#include <utility>

namespace Kratos
{
template <unsigned int MatrixSize>
class StiffnessCalculator : public ContributionCalculator<MatrixSize>
{
public:
    using BaseType = ContributionCalculator<MatrixSize>;

    struct InputProvider {
        InputProvider(Geo::BMatricesGetter               GetBMatrices,
                      Geo::StrainVectorsGetter           GetStrains,
                      Geo::IntegrationCoefficientsGetter GetIntegrationCoefficients,
                      Geo::PropertiesGetter              GetElementProperties,
                      Geo::ProcessInfoGetter             GetProcessInfo,
                      Geo::ConstitutiveLawsGetter        GetConstitutiveLaws)
            : GetBMatrices(std::move(GetBMatrices)),
              GetStrains(std::move(GetStrains)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetElementProperties(std::move(GetElementProperties)),
              GetProcessInfo(std::move(GetProcessInfo)),
              GetConstitutiveLaws(std::move(GetConstitutiveLaws))
        {
        }

        Geo::BMatricesGetter               GetBMatrices;
        Geo::StrainVectorsGetter           GetStrains;
        Geo::IntegrationCoefficientsGetter GetIntegrationCoefficients;
        Geo::PropertiesGetter              GetElementProperties;
        Geo::ProcessInfoGetter             GetProcessInfo;
        Geo::ConstitutiveLawsGetter        GetConstitutiveLaws;
    };

    explicit StiffnessCalculator(InputProvider StiffnessInputProvider)
        : mInputProvider(std::move(StiffnessInputProvider))
    {
    }

    std::optional<typename BaseType::LHSMatrixType> LHSContribution() override
    {
        return std::make_optional(GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
            mInputProvider.GetBMatrices(), CalculateConstitutiveMatricesAtIntegrationPoints(),
            mInputProvider.GetIntegrationCoefficients()));
    }

    typename BaseType::RHSVectorType RHSContribution() override
    {
        const auto stresses = StressStrainUtilities::CalculateStressVectorsFromStrainVectors(
            mInputProvider.GetStrains(), mInputProvider.GetProcessInfo(),
            mInputProvider.GetElementProperties(), mInputProvider.GetConstitutiveLaws());
        auto result = typename BaseType::RHSVectorType{
            -1.0 * GeoEquationOfMotionUtilities::CalculateInternalForceVector(
                       mInputProvider.GetBMatrices(), stresses, mInputProvider.GetIntegrationCoefficients())};
        KRATOS_INFO("StiffnessCalculator::RHSContribution") << result << std::endl;
        return result;
    }

    std::pair<std::optional<typename BaseType::LHSMatrixType>, typename BaseType::RHSVectorType> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    }

private:
    std::vector<Matrix> CalculateConstitutiveMatricesAtIntegrationPoints()
    {
        const auto& r_properties                  = mInputProvider.GetElementProperties();
        const auto& r_process_info                = mInputProvider.GetProcessInfo();
        auto        calculate_constitutive_matrix = [&r_properties, &r_process_info](
                                                 const auto& rpConstitutiveLaw, auto rRelativeDisplacement) {
            auto result = Matrix{rpConstitutiveLaw->GetStrainSize(), rpConstitutiveLaw->GetStrainSize()};
            auto law_parameters = ConstitutiveLaw::Parameters{};
            law_parameters.SetMaterialProperties(r_properties);
            law_parameters.SetStrainVector(rRelativeDisplacement);
            law_parameters.SetProcessInfo(r_process_info);
            rpConstitutiveLaw->CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, result);
            return result;
        };
        auto result = std::vector<Matrix>{};
        result.reserve(mInputProvider.GetConstitutiveLaws().size());
        std::ranges::transform(mInputProvider.GetConstitutiveLaws(), mInputProvider.GetStrains(),
                               std::back_inserter(result), calculate_constitutive_matrix);

        return result;
    }

    InputProvider mInputProvider;
};

} // namespace Kratos
