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

#pragma once

#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{

class BoundaryGenerator
{
public:
    virtual ~BoundaryGenerator() = default;
    virtual Geometry<Node>::GeometriesArrayType operator()(const Geometry<Node>& rGeometry) const = 0;
};

class EdgesGenerator : public BoundaryGenerator
{
public:
    Geometry<Node>::GeometriesArrayType operator()(const Geometry<Node>& rGeometry) const override;
};

class FacesGenerator : public BoundaryGenerator
{
public:
    Geometry<Node>::GeometriesArrayType operator()(const Geometry<Node>& rGeometry) const override;
};

} // namespace Kratos
