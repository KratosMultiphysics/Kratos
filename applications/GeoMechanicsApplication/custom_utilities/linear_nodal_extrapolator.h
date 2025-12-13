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

#include "geo_aliases.h"
#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"
#include "nodal_extrapolator.h"

#include <memory>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LinearNodalExtrapolator : public NodalExtrapolator
{
public:
    using NodalExtrapolator::GeometryType;

    [[nodiscard]] Matrix CalculateElementExtrapolationMatrix(
        const GeometryType& rGeometry, const Geo::IntegrationPointVectorType& rIntegrationPoints) const override;

private:
    void static CheckIfGeometryIsSupported(const GeometryType& rGeometry);
    [[nodiscard]] static std::unique_ptr<GeometryType> CreateLowerOrderGeometry(const GeometryType& rGeometry);
    static void AddRowsForMidsideNodes(const GeometryType& rGeometry, Matrix& rExtrapolationMatrix);
    [[nodiscard]] static Matrix CalculateExtrapolationMatrixForCornerNodes(const GeometryType& rGeometry,
                                                                           const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                                           const GeometryType& rCornerGeometry);
};

} // namespace Kratos
