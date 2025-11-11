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
    try {
        return rGeometry.GenerateEdges();
    }
    catch (const Exception& e) {
        // Since some of the surface geometries in core do not implement GenerateEdges, but we
        // know how to recover if the local dimension is 1, we do so here.
        if (rGeometry.LocalSpaceDimension() == 1) {
            Geometry<Node>::GeometriesArrayType result;
            result.push_back(std::make_shared<Geometry<Node>>(rGeometry));
            return result;
        }

        throw;
    }
}

Geometry<Node>::GeometriesArrayType FacesGenerator::operator()(const Geometry<Node>& rGeometry) const
{
    try {
        return rGeometry.GenerateFaces();
    }
    catch (const Exception& e) {
        // Since some of the surface geometries in core do not implement GenerateFaces, but we
        // know how to recover if the local dimension is 2, we do so here.
        if (rGeometry.LocalSpaceDimension() == 2) {
            Geometry<Node>::GeometriesArrayType result;
            result.push_back(std::make_shared<Geometry<Node>>(rGeometry));
            return result;
        }

        throw;
    }
}

} // namespace Kratos
