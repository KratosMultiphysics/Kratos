// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Richard Faasse

#pragma once

#include "contribution_calculator.h"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
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
    struct InputProvider {
        InputProvider(std::function<std::vector<Matrix>()> GetBMatrices,
                      std::function<std::vector<Vector>()> GetRelativeDisplacements,
                      std::function<std::vector<double>()> GetIntegrationCoefficients,
                      std::function<const Properties&()>   GetElementProperties,
                      std::function<const ProcessInfo&()>  GetProcessInfo,
                      std::function<const std::vector<ConstitutiveLaw::Pointer>&()> GetConstitutiveLaws)
            : GetBMatrices(std::move(GetBMatrices)),
              GetRelativeDisplacements(std::move(GetRelativeDisplacements)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetElementProperties(std::move(GetElementProperties)),
              GetProcessInfo(std::move(GetProcessInfo)),
              GetConstitutiveLaws(std::move(GetConstitutiveLaws))
        {
        }

        std::function<std::vector<Matrix>()>                          GetBMatrices;
        std::function<std::vector<Vector>()>                          GetRelativeDisplacements;
        std::function<std::vector<double>()>                          GetIntegrationCoefficients;
        std::function<const Properties&()>                            GetElementProperties;
        std::function<const ProcessInfo&()>                           GetProcessInfo;
        std::function<const std::vector<ConstitutiveLaw::Pointer>&()> GetConstitutiveLaws;
    };

    explicit StiffnessCalculator(InputProvider AnInputProvider)
        : mInputProvider(std::move(AnInputProvider))
    {
    }

    std::optional<BoundedMatrix<double, MatrixSize, MatrixSize>> LHSContribution() override
    {
        return std::make_optional(GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
            mInputProvider.GetBMatrices(),
            CalculateConstitutiveMatricesAtIntegrationPoints(
                mInputProvider.GetConstitutiveLaws(), mInputProvider.GetElementProperties(),
                mInputProvider.GetRelativeDisplacements(), mInputProvider.GetProcessInfo()),
            mInputProvider.GetIntegrationCoefficients()));
    }

    BoundedVector<double, MatrixSize> RHSContribution() override
    {
        const auto local_b_matrices       = mInputProvider.GetBMatrices();
        const auto relative_displacements = mInputProvider.GetRelativeDisplacements();
        const auto tractions = StressStrainUtilities::CalculateStressFromStrain(
            relative_displacements, mInputProvider.GetProcessInfo(),
            mInputProvider.GetElementProperties(), mInputProvider.GetConstitutiveLaws());
        const auto integration_coefficients = mInputProvider.GetIntegrationCoefficients();
        return BoundedVector<double, MatrixSize>{-GeoEquationOfMotionUtilities::CalculateInternalForceVector(
            local_b_matrices, tractions, integration_coefficients)};
    }

    std::pair<std::optional<BoundedMatrix<double, MatrixSize, MatrixSize>>, BoundedVector<double, MatrixSize>> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    }

private:
    std::vector<Matrix> CalculateConstitutiveMatricesAtIntegrationPoints(
        const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws,
        const Properties&                            rProperties,
        const std::vector<Vector>&                   rRelativeDisplacements,
        const ProcessInfo&                           rProcessInfo)
    {
        auto get_constitutive_matrix = [&rProperties, &rProcessInfo](const auto& p_constitutive_law,
                                                                     auto rRelativeDisplacement) {
            auto result = Matrix{p_constitutive_law->GetStrainSize(), p_constitutive_law->GetStrainSize()};
            auto law_parameters = ConstitutiveLaw::Parameters{};
            law_parameters.SetMaterialProperties(rProperties);
            law_parameters.SetStrainVector(rRelativeDisplacement);
            law_parameters.SetProcessInfo(rProcessInfo);
            p_constitutive_law->CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, result);
            return result;
        };
        auto result = std::vector<Matrix>{};
        result.reserve(rConstitutiveLaws.size());
        std::transform(rConstitutiveLaws.begin(), rConstitutiveLaws.end(), rRelativeDisplacements.begin(),
                       std::back_inserter(result), get_constitutive_matrix);

        return result;
    }

    InputProvider mInputProvider;
};

} // namespace Kratos
