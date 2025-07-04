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
#include "custom_retention/retention_law.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

#include <utility>
#include <vector>

namespace Kratos
{

class CompressibilityCalculator : public ContributionCalculator
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<const Matrix&()>                             GetNContainer,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<double()>                        GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)
            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetNContainer(std::move(GetNContainer)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              GetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        std::function<const Properties&()>                         GetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws;
        std::function<const Matrix&()>                             GetNContainer;
        std::function<Vector()>                                    GetIntegrationCoefficients;
        std::function<double()>                                    GetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)>             GetNodalValues;
    };

    explicit CompressibilityCalculator(InputProvider rInputProvider);

    std::optional<Matrix>                    LHSContribution() override;
    Vector                                   RHSContribution() override;
    std::pair<std::optional<Matrix>, Vector> LocalSystemContribution() override;

private:
    [[nodiscard]] Matrix CalculateCompressibilityMatrix() const;
    [[nodiscard]] double CalculateBiotModulusInverse(const RetentionLaw::Pointer& rRetentionLaw,
                                                     double FluidPresssure) const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rCompressibilityMatrix) const;
    [[nodiscard]] Matrix LHSContribution(const Matrix& rCompressibilityMatrix) const;

    InputProvider mInputProvider;
};

} // namespace Kratos
