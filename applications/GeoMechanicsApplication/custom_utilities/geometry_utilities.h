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

#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "includes/kratos_export_api.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeometryUtilities
{
public:
    static Matrix Calculate2DRotationMatrixForLineGeometry(const Geometry<Node>& rGeometry,
                                                           const array_1d<double, 3>& rLocalCoordinate);
    static Matrix Calculate3DRotationMatrixForPlaneGeometry(const Geometry<Node>& rGeometry,
                                                            const array_1d<double, 3>& rLocalCoordinate);
    static void   ReverseNodes(PointerVector<Node>&                  rNodes,
                               GeometryData::KratosGeometryFamily    GeometryFamily,
                               GeometryData::KratosGeometryOrderType GeometryOrderType);
    static void   ReverseNodes(std::vector<std::size_t>&             rNodeIds,
                               GeometryData::KratosGeometryFamily    GeometryFamily,
                               GeometryData::KratosGeometryOrderType GeometryOrderType);
};

} // namespace Kratos
