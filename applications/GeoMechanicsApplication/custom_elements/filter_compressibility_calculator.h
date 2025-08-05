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
template <unsigned int TNumNodes>
class FilterCompressibilityCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()>   GetElementProperties,
                      std::function<const Matrix&()>       GetNContainer,
                      std::function<Vector()>              GetIntegrationCoefficients,
                      std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints,
                      std::function<double()>              GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)
            : GetElementProperties(std::move(GetElementProperties)),
              GetNContainer(std::move(GetNContainer)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityAtIntegrationPoints(std::move(GetProjectedGravityAtIntegrationPoints)),
              GetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              GetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        std::function<const Properties&()>             GetElementProperties;
        std::function<const Matrix&()>                 GetNContainer;
        std::function<Vector()>                        GetIntegrationCoefficients;
        std::function<std::vector<Vector>()>           GetProjectedGravityAtIntegrationPoints;
        std::function<double()>                        GetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)> GetNodalValues;
    };

    explicit FilterCompressibilityCalculator(InputProvider AnInputProvider);

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override;
    BoundedVector<double, TNumNodes>                           RHSContribution() override;
    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override;

private:
    [[nodiscard]] BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCompressibilityMatrix() const;
    [[nodiscard]] double CalculateElasticCapacity(double ProjectedGravity) const;
    [[nodiscard]] BoundedVector<double, TNumNodes> RHSContribution(const Matrix& rCompressibilityMatrix) const;
    [[nodiscard]] BoundedMatrix<double, TNumNodes, TNumNodes> LHSContribution(const Matrix& rCompressibilityMatrix) const;

    InputProvider mInputProvider;
};

} // namespace Kratos
