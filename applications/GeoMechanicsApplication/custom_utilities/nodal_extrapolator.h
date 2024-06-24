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
#include "includes/node.h"
#include <cstddef>

namespace Kratos
{

class NodalExtrapolator
{
public:
    using GeometryType = Geometry<Node>;
    using SizeType     = std::size_t;
    using IndexType    = std::size_t;

    [[nodiscard]] Matrix CalculateElementExtrapolationMatrix(const GeometryType& rGeometry,
                                                             const GeometryData::IntegrationMethod& rIntegrationMethod) const;

private:
    void CheckIfGeometryIsSupported(const GeometryType& rGeometry) const;
    [[nodiscard]] std::unique_ptr<GeometryType> CreateLowerOrderGeometry(const GeometryType& rGeometry) const;
    void AddRowsForMidsideNodes(const GeometryType& rGeometry, Matrix& rExtrapolationMatrix) const;
    [[nodiscard]] Matrix CalculateExtrapolationMatrixForCornerNodes(const GeometryType& rGeometry,
                                                                    const GeometryData::IntegrationMethod& rIntegrationMethod,
                                                                    const GeometryType& rCornerGeometry) const;
};

} // namespace Kratos
