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
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<const Matrix&()>                             GetNContainer,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<double()>                        GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)
            : mGetElementProperties(std::move(GetElementProperties)),
              mGetRetentionLaws(std::move(GetRetentionLaws)),
              mGetNContainer(std::move(GetNContainer)),
              mGetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              mGetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              mGetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        [[nodiscard]] const Properties&                         GetElementProperties() const;
        [[nodiscard]] const std::vector<RetentionLaw::Pointer>& GetRetentionLaws() const;
        [[nodiscard]] const Matrix&                             GetNContainer() const;
        [[nodiscard]] Vector                                    GetIntegrationCoefficients() const;
        [[nodiscard]] double                                    GetMatrixScalarFactor() const;
        [[nodiscard]] Vector GetNodalValues(const Variable<double>& rVariable) const;

    private:
        std::function<const Properties&()>                         mGetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()> mGetRetentionLaws;
        std::function<const Matrix&()>                             mGetNContainer;
        std::function<Vector()>                                    mGetIntegrationCoefficients;
        std::function<double()>                                    mGetMatrixScalarFactor;
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
