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

class FluidBodyFlowCalculator : public ContributionCalculator
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<Vector()> GetIntegrationCoefficients,
                      std::function<Vector()> GetProjectedGravityForIntegrationPoints,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::size_t()> GetLocalSpaceDimension)

            : GetElementProperties(std::move(GetElementProperties)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityForIntegrationPoints(std::move(GetProjectedGravityForIntegrationPoints)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetLocalSpaceDimension(std::move(GetLocalSpaceDimension))
        {
        }

        std::function<const Properties&()> GetElementProperties;
        std::function<Vector()>            GetIntegrationCoefficients;
        std::function<Vector()>            GetProjectedGravityForIntegrationPoints;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
        std::function<std::size_t()>                                 GetLocalSpaceDimension;
    };

    explicit FluidBodyFlowCalculator(InputProvider AnInputProvider);

    std::optional<Matrix>                                   LHSContribution() override;
    std::optional<Vector>                                   RHSContribution() override;
    std::pair<std::optional<Matrix>, std::optional<Vector>> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;
};

} // namespace Kratos
