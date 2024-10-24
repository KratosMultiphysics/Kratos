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

#include "calculator.h"

#include <custom_retention/retention_law.h>

#include <utility>

namespace Kratos
{

class PermeabilityCalculator : public Calculator
{
public:
    class InputProvider
    {
    public:
        InputProvider(std::function<const Properties&()>                         rElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> rRetentionLaws,
                      std::function<Vector()>                        rIntegrationCoefficients,
                      std::function<Vector(const Variable<double>&)> rNodalValuesOfDtWaterPressure,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> rShapeFunctionGradients)
            : mGetElementProperties(std::move(rElementProperties)),
              mGetRetentionLaws(std::move(rRetentionLaws)),
              mGetIntegrationCoefficients(std::move(rIntegrationCoefficients)),
              mGetNodalValues(std::move(rNodalValuesOfDtWaterPressure)),
              mGetShapeFunctionGradients(std::move(rShapeFunctionGradients))
        {
        }

        [[nodiscard]] const Properties& GetElementProperties() const
        {
            return mGetElementProperties();
        }

        [[nodiscard]] const std::vector<RetentionLaw::Pointer>& GetRetentionLaws() const
        {
            return mGetRetentionLaws();
        }

        [[nodiscard]] Vector GetIntegrationCoefficients() const
        {
            return mGetIntegrationCoefficients();
        }

        [[nodiscard]] Vector GetNodalValues(const Variable<double>& rVariable) const
        {
            return mGetNodalValues(rVariable);
        }

        [[nodiscard]] Geometry<Node>::ShapeFunctionsGradientsType GetShapeFunctionGradients() const
        {
            return mGetShapeFunctionGradients();
        }

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
    std::pair<Matrix, Vector> CalculateLeftAndRightHandSide() override;

private:
    InputProvider mInputProvider;

    [[nodiscard]] Matrix CalculatePermeabilityMatrix() const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rPermeabilityMatrix) const;
};

} // namespace Kratos
