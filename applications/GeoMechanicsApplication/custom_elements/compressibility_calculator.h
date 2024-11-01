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

class CompressibilityCalculator : public ContributionCalculator
{
public:
    class InputProvider
    {
    public:
        InputProvider(std::function<const Properties&()>                         rElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> rRetentionLaws,
                      std::function<const Matrix&()>                             rNContainer,
                      std::function<Vector()>                        rIntegrationCoefficients,
                      std::function<double()>                        DtPressureCoefficient,
                      std::function<Vector(const Variable<double>&)> rNodalValuesOfDtWaterPressure)
            : mGetElementProperties(std::move(rElementProperties)),
              mGetRetentionLaws(std::move(rRetentionLaws)),
              mGetNContainer(std::move(rNContainer)),
              mGetIntegrationCoefficients(std::move(rIntegrationCoefficients)),
              mGetDtPressureCoefficient(std::move(DtPressureCoefficient)),
              mGetNodalValues(std::move(rNodalValuesOfDtWaterPressure))
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

        [[nodiscard]] const Matrix& GetNContainer() const { return mGetNContainer(); }

        [[nodiscard]] Vector GetIntegrationCoefficients() const
        {
            return mGetIntegrationCoefficients();
        }

        [[nodiscard]] double GetDtPressureCoefficient() const
        {
            return mGetDtPressureCoefficient();
        }

        [[nodiscard]] Vector GetNodalValues(const Variable<double>& rVariable) const
        {
            return mGetNodalValues(rVariable);
        }

    private:
        std::function<const Properties&()>                         mGetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()> mGetRetentionLaws;
        std::function<const Matrix&()>                             mGetNContainer;
        std::function<Vector()>                                    mGetIntegrationCoefficients;
        std::function<double()>                                    mGetDtPressureCoefficient;
        std::function<Vector(const Variable<double>&)>             mGetNodalValues;
    };

    explicit CompressibilityCalculator(InputProvider rInputProvider);

    Matrix                    LHSContribution() override;
    Vector                    RHSContribution() override;
    std::pair<Matrix, Vector> LocalSystemContribution() override;

private:
    [[nodiscard]] Matrix CalculateCompressibilityMatrix() const;
    [[nodiscard]] double CalculateBiotModulusInverse(const RetentionLaw::Pointer& rRetentionLaw) const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rCompressibilityMatrix) const;
    [[nodiscard]] Matrix LHSContribution(const Matrix& rCompressibilityMatrix) const;

    InputProvider mInputProvider;
};

} // namespace Kratos
