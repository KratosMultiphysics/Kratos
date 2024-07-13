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

#include "geometries/geometry.h"
#include "includes/define.h"
#include "nodal_extrapolator.h"
#include <cstddef>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LinearNodalExtrapolator : public NodalExtrapolator
{
public:
    using NodalExtrapolator::GeometryType;
    using SizeType  = std::size_t;
    using IndexType = std::size_t;

    [[nodiscard]] Matrix CalculateElementExtrapolationMatrix(
        const GeometryType& rGeometry, const GeometryData::IntegrationMethod& rIntegrationMethod) const override;

private:
    void static CheckIfGeometryIsSupported(const GeometryType& rGeometry);
    [[nodiscard]] static std::unique_ptr<GeometryType> CreateLowerOrderGeometry(const GeometryType& rGeometry);
    static void AddRowsForMidsideNodes(const GeometryType& rGeometry, Matrix& rExtrapolationMatrix);
    [[nodiscard]] static Matrix CalculateExtrapolationMatrixForCornerNodes(const GeometryType& rGeometry,
                                                                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                                                                           const GeometryType& rCornerGeometry);
};

} // namespace Kratos
