//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//

// Project includes
#include "qef_utility.h"

namespace Kratos {
    array_1d<double,3> QEF::CalculateCenter(const GeometryType& rVoxel) {
        PointsArrayType nodes = rVoxel.Points();
        double x = (nodes[0].X() + nodes[1].X())/2.0;
        double y = (nodes[1].Y() + nodes[2].Y())/2.0;
        double z = (nodes[0].Z() + nodes[4].Z())/2.0;
        array_1d<double,3> center({x,y,z});
        return center;
    }
}