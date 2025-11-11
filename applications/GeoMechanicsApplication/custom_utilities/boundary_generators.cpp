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
//                   Anne van de Graaf
//

#include "boundary_generators.h"

namespace Kratos
{

Geometry<Node>::GeometriesArrayType PointsGenerator::operator()(const Geometry<Node>& rGeometry) const
{
    return rGeometry.GeneratePoints();
}

Geometry<Node>::GeometriesArrayType EdgesGenerator::operator()(const Geometry<Node>& rGeometry) const
{
    return rGeometry.GenerateEdges();
}

Geometry<Node>::GeometriesArrayType FacesGenerator::operator()(const Geometry<Node>& rGeometry) const
{
    return rGeometry.GenerateFaces();
}

} // namespace Kratos
