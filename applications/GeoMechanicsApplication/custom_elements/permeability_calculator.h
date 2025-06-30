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

#include <utility>

namespace Kratos
{

class PermeabilityCalculator : public ContributionCalculator
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients)
            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetNodalValues(std::move(GetNodalValuesOf)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients))
        {
        }

        std::function<const Properties&()>                           GetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Vector()>                                      GetIntegrationCoefficients;
        std::function<Vector(const Variable<double>&)>               GetNodalValues;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
    };

    explicit PermeabilityCalculator(InputProvider InputProvider);

    std::optional<Matrix>                    LHSContribution() override;
    Vector                                   RHSContribution() override;
    std::pair<std::optional<Matrix>, Vector> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;

    [[nodiscard]] Matrix CalculatePermeabilityMatrix() const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rPermeabilityMatrix) const;
};

} // namespace Kratos
