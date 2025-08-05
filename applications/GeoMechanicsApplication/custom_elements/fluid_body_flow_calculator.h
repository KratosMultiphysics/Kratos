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

template <unsigned int TNumNodes>
class FluidBodyFlowCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<Vector()>              GetIntegrationCoefficients,
                      std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::size_t()>         GetLocalSpaceDimension,
                      std::function<std::vector<double>()> GetFluidPressures)

            : GetElementProperties(std::move(GetElementProperties)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityAtIntegrationPoints(std::move(GetProjectedGravityAtIntegrationPoints)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetLocalSpaceDimension(std::move(GetLocalSpaceDimension)),
              GetFluidPressures(std::move(GetFluidPressures))
        {
        }

        std::function<const Properties&()>   GetElementProperties;
        std::function<Vector()>              GetIntegrationCoefficients;
        std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
        std::function<std::size_t()>                                 GetLocalSpaceDimension;
        std::function<std::vector<double>()>                         GetFluidPressures;
    };

    explicit FluidBodyFlowCalculator(InputProvider AnInputProvider);

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override;
    BoundedVector<double, TNumNodes>                           RHSContribution() override;
    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;
    std::vector<double> CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const;
};

} // namespace Kratos
