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
            : mGetElementProperties(rElementProperties),
              mGetRetentionLaws(rRetentionLaws),
              mGetIntegrationCoefficients(rIntegrationCoefficients),
              mGetNodalValues(rNodalValuesOfDtWaterPressure),
              mGetShapeFunctionGradients(rShapeFunctionGradients)
        {
        }

        const Properties& GetElementProperties() const { return mGetElementProperties(); }

        const std::vector<RetentionLaw::Pointer>& GetRetentionLaws() const
        {
            return mGetRetentionLaws();
        }

        Vector GetIntegrationCoefficients() const { return mGetIntegrationCoefficients(); }

        Vector GetNodalValues(const Variable<double>& rVariable) const
        {
            return mGetNodalValues(rVariable);
        }

        Geometry<Node>::ShapeFunctionsGradientsType GetShapeFunctionGradients() const
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

    PermeabilityCalculator(InputProvider rInputProvider);

    Matrix LHSContribution() override;
    Vector                    RHSContribution() override;
    std::pair<Matrix, Vector> CalculateLeftAndRightHandSide() override;

private:
    InputProvider mInputProvider;

    Matrix CalculatePermeabilityMatrix() const;
};

} // namespace Kratos
