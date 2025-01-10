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

class FilterCompressibilityCalculator : public ContributionCalculator
{
public:
    class InputProvider
    {
    public:
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const Matrix&()>                             GetNContainer,
                      std::function<Vector()> GetIntegrationCoefficients,
                      std::function<Vector()> GetProjectedGravityForIntegrationPoints,
                      std::function<double()> GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)
            : mGetElementProperties(std::move(GetElementProperties)),
              mGetNContainer(std::move(GetNContainer)),
              mGetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              mGetProjectedGravityForIntegrationPoints(std::move(GetProjectedGravityForIntegrationPoints)),
              mGetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              mGetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        [[nodiscard]] const Properties&                         GetElementProperties() const;
        [[nodiscard]] const Matrix&                             GetNContainer() const;
        [[nodiscard]] Vector                                    GetIntegrationCoefficients() const;
        [[nodiscard]] Vector GetProjectedGravityForIntegrationPoints() const;
        [[nodiscard]] double GetMatrixScalarFactor() const;
        [[nodiscard]] Vector GetNodalValues(const Variable<double>& rVariable) const;

    private:
        std::function<const Properties&()>                         mGetElementProperties;
        std::function<const Matrix&()>                             mGetNContainer;
        std::function<Vector()>                                    mGetIntegrationCoefficients;
        std::function<Vector()>                        mGetProjectedGravityForIntegrationPoints;
        std::function<double()>                        mGetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)> mGetNodalValues;
    };

    explicit FilterCompressibilityCalculator(InputProvider rInputProvider);

    Matrix                    LHSContribution() override;
    Vector                    RHSContribution() override;
    std::pair<Matrix, Vector> LocalSystemContribution() override;

private:
    [[nodiscard]] Matrix CalculateCompressibilityMatrix() const;
    [[nodiscard]] double CalculateBiotModulusInverse(double ProjectedGravity) const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rCompressibilityMatrix) const;
    [[nodiscard]] Matrix LHSContribution(const Matrix& rCompressibilityMatrix) const;

    InputProvider mInputProvider;
};

} // namespace Kratos
