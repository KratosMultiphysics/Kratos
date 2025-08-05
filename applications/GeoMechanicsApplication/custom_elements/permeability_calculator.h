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

template <unsigned int TNumNodes>
class PermeabilityCalculator : public ContributionCalculator<TNumNodes> ///
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::vector<double>()> GetFluidPressures)
            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetNodalValues(std::move(GetNodalValuesOf)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetFluidPressures(std::move(GetFluidPressures))
        {
        }

        std::function<const Properties&()>                           GetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Vector()>                                      GetIntegrationCoefficients;
        std::function<Vector(const Variable<double>&)>               GetNodalValues;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
        std::function<std::vector<double>()>                         GetFluidPressures;
    };

    explicit PermeabilityCalculator(InputProvider InputProvider);

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override;
    BoundedVector<double, TNumNodes>                           RHSContribution() override;
    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;

    [[nodiscard]] BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix() const;
    [[nodiscard]] BoundedVector<double, TNumNodes> RHSContribution(const Matrix& rPermeabilityMatrix) const;
};

} // namespace Kratos
