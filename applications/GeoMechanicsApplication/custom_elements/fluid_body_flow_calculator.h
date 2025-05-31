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
                      std::function<Vector()>              GetIntegrationCoefficients,
                      std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints,
                      std::function<const Matrix&()>       GetNContainer,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::size_t()>                   GetLocalSpaceDimension,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)

            : GetElementProperties(std::move(GetElementProperties)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityAtIntegrationPoints(std::move(GetProjectedGravityAtIntegrationPoints)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetNContainer(std::move(GetNContainer)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetLocalSpaceDimension(std::move(GetLocalSpaceDimension)),
              GetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        std::function<const Properties&()>   GetElementProperties;
        std::function<Vector()>              GetIntegrationCoefficients;
        std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<const Matrix&()>                               GetNContainer;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
        std::function<std::size_t()>                                 GetLocalSpaceDimension;
        std::function<Vector(const Variable<double>&)>               GetNodalValues;
    };

    explicit FluidBodyFlowCalculator(InputProvider AnInputProvider);

    std::optional<Matrix>                    LHSContribution() override;
    Vector                                   RHSContribution() override;
    std::pair<std::optional<Matrix>, Vector> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;
};

} // namespace Kratos
