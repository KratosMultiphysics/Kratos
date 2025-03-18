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
                      std::function<const Matrix&()>     GetNContainer,
                      std::function<Vector()>            GetIntegrationCoefficients,
                      std::function<Vector()>            GetProjectedGravityForIntegrationPoints,
                      std::function<double()>            GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)>             GetNodalValuesOf,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients)

            : GetElementProperties(std::move(GetElementProperties)),
              GetNContainer(std::move(GetNContainer)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityForIntegrationPoints(std::move(GetProjectedGravityForIntegrationPoints)),
              GetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              GetNodalValues(std::move(GetNodalValuesOf)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients))
        {
        }

        std::function<const Properties&()>             GetElementProperties;
        std::function<const Matrix&()>                 GetNContainer;
        std::function<Vector()>                        GetIntegrationCoefficients;
        std::function<Vector()>                        GetProjectedGravityForIntegrationPoints;
        std::function<double()>                        GetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)> GetNodalValues;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
    };

    explicit FluidBodyFlowCalculator(InputProvider AnInputProvider);

    Matrix                    LHSContribution() override;
    Vector                    RHSContribution() override;
    std::pair<Matrix, Vector> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;
};

} // namespace Kratos
