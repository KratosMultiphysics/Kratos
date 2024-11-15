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
    class InputProvider
    {
    public:
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients)
            : mGetElementProperties(std::move(GetElementProperties)),
              mGetRetentionLaws(std::move(GetRetentionLaws)),
              mGetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              mGetNodalValues(std::move(GetNodalValuesOf)),
              mGetShapeFunctionGradients(std::move(GetShapeFunctionGradients))
        {
        }

        [[nodiscard]] const Properties&                         GetElementProperties() const;
        [[nodiscard]] const std::vector<RetentionLaw::Pointer>& GetRetentionLaws() const;
        [[nodiscard]] Vector                                    GetIntegrationCoefficients() const;
        [[nodiscard]] Vector GetNodalValues(const Variable<double>& rVariable) const;
        [[nodiscard]] Geometry<Node>::ShapeFunctionsGradientsType GetShapeFunctionGradients() const;

    private:
        std::function<const Properties&()>                           mGetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   mGetRetentionLaws;
        std::function<Vector()>                                      mGetIntegrationCoefficients;
        std::function<Vector(const Variable<double>&)>               mGetNodalValues;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> mGetShapeFunctionGradients;
    };

    explicit PermeabilityCalculator(InputProvider InputProvider);

    Matrix                    LHSContribution() override;
    Vector                    RHSContribution() override;
    std::pair<Matrix, Vector> LocalSystemContribution() override;

private:
    InputProvider mInputProvider;

    [[nodiscard]] Matrix CalculatePermeabilityMatrix() const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rPermeabilityMatrix) const;
};

} // namespace Kratos
