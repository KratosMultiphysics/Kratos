//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alireza 
//

#ifndef KRATOS_INTERSECTION_POINTS_CONTAINER_H_INCLUDED
#define KRATOS_INTERSECTION_POINTS_CONTAINER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "containers/array_1d.h"
#include <vector>

namespace Kratos
{
    struct IntersectionPointData {
        unsigned int elementId;
        unsigned int pointId;
        array_1d<double, 3> coordinates;
    };

    namespace KratosDropletDynamics 
    {
        extern std::vector<IntersectionPointData> g_IntersectionPointsContainer;
        
    }
}

#endif // KRATOS_INTERSECTION_POINTS_CONTAINER_H_INCLUDED
