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
#include "custom_retention/retention_law.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

#include <vector>

namespace Kratos
{

class CompressibilityCalculator : public Calculator
{
public:
    class InputProvider
    {
    public:
        InputProvider(std::function<const Properties&()>                         rElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> rRetentionLaws,
                      std::function<const Matrix&()>                             rNContainer,
                      std::function<Vector()>                        rIntegrationCoefficients,
                      std::function<const double()>                  DtPressureCoefficient,
                      std::function<Vector(const Variable<double>&)> rNodalValuesOfDtWaterPressure)
            : mGetElementProperties(rElementProperties),
              mGetRetentionLaws(rRetentionLaws),
              mGetNContainer(rNContainer),
              mGetIntegrationCoefficients(rIntegrationCoefficients),
              mGetDtPressureCoefficient(DtPressureCoefficient),
              mGetNodalValues(rNodalValuesOfDtWaterPressure)
        {
        }

        const Properties& GetElementProperties() const { return mGetElementProperties(); }

        const std::vector<RetentionLaw::Pointer>& GetRetentionLaws() const
        {
            return mGetRetentionLaws();
        }

        const Matrix& GetNContainer() const { return mGetNContainer(); }

        Vector GetIntegrationCoefficients() const { return mGetIntegrationCoefficients(); }

        double GetDtPressureCoefficient() const { return mGetDtPressureCoefficient(); }

        Vector GetNodalValues(const Variable<double>& rVariable) const
        {
            return mGetNodalValues(rVariable);
        }

    private:
        std::function<const Properties&()>                         mGetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()> mGetRetentionLaws;
        std::function<const Matrix&()>                             mGetNContainer;
        std::function<Vector()>                                    mGetIntegrationCoefficients;
        std::function<const double()>                              mGetDtPressureCoefficient;
        std::function<Vector(const Variable<double>&)>             mGetNodalValues;
    };

    explicit CompressibilityCalculator(const InputProvider& rInputProvider);

    Matrix LHSContribution() override;
    Vector RHSContribution() override;
    void CalculateLeftAndRightHandSide(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector) override;

private:
    Matrix CompressibilityCalculator::CalculateCompressibilityMatrix(const Matrix& rNContainer,
                                                                     const Vector& rIntegrationCoefficients);
    double CompressibilityCalculator::CalculateBiotModulusInverse(const unsigned int integrationPointIndex) const;

    InputProvider mInputProvider;
};

} // namespace Kratos
