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
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const Matrix&()>     GetNContainer,
                      std::function<Vector()>            GetIntegrationCoefficients,
                      std::function<std::vector<Vector>()>            GetProjectedGravityForIntegrationPoints,
                      std::function<double()>            GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)
            : GetElementProperties(std::move(GetElementProperties)),
              GetNContainer(std::move(GetNContainer)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityForIntegrationPoints(std::move(GetProjectedGravityForIntegrationPoints)),
              GetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              GetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        std::function<const Properties&()>             GetElementProperties;
        std::function<const Matrix&()>                 GetNContainer;
        std::function<Vector()>                        GetIntegrationCoefficients;
        std::function<std::vector<Vector>()>                        GetProjectedGravityForIntegrationPoints;
        std::function<double()>                        GetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)> GetNodalValues;
    };

    explicit FilterCompressibilityCalculator(InputProvider AnInputProvider);

    std::optional<Matrix>                                   LHSContribution() override;
    std::optional<Vector>                                   RHSContribution() override;
    std::pair<std::optional<Matrix>, std::optional<Vector>> LocalSystemContribution() override;

private:
    [[nodiscard]] Matrix CalculateCompressibilityMatrix() const;
    [[nodiscard]] double CalculateElasticCapacity(double ProjectedGravity) const;
    [[nodiscard]] Vector RHSContribution(const Matrix& rCompressibilityMatrix) const;
    [[nodiscard]] Matrix LHSContribution(const Matrix& rCompressibilityMatrix) const;

    InputProvider mInputProvider;
};

} // namespace Kratos
