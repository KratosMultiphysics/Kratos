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
#include "geometries/geometry.h"
#include "includes/define.h"
#include "nodal_extrapolator.h"

#include <cstddef>

namespace Kratos
{

class Element;

class KRATOS_API(GEO_MECHANICS_APPLICATION) LinearNodalExtrapolator : public NodalExtrapolator
{
public:
    using NodalExtrapolator::GeometryType;
    using SizeType  = std::size_t;
    using IndexType = std::size_t;

    [[nodiscard]] Matrix CalculateElementExtrapolationMatrix(const Element& rElement) const override;

private:
    void static CheckIfGeometryIsSupported(const GeometryType& rGeometry);
    [[nodiscard]] static std::unique_ptr<GeometryType> CreateLowerOrderGeometry(const GeometryType& rGeometry);
    static void AddRowsForMidsideNodes(const GeometryType& rGeometry, Matrix& rExtrapolationMatrix);
    [[nodiscard]] static Matrix CalculateExtrapolationMatrixForCornerNodes(const GeometryType& rGeometry,
                                                                           const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                                           const GeometryType& rCornerGeometry);
};

} // namespace Kratos
