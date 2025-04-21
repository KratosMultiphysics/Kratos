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

    struct InterfaceAverageData
{
    int elementId;
    array_1d<double, 3> averageCoordinates;
    array_1d<double, 3> averageNormal;
    double interfaceArea;
    int numberOfPoints;
    
    InterfaceAverageData() : 
        elementId(0), 
        interfaceArea(0.0), 
        numberOfPoints(0) 
    {
        averageCoordinates = ZeroVector(3);
        averageNormal = ZeroVector(3);
    }
};

    namespace KratosDropletDynamics 
    {
        extern std::vector<IntersectionPointData> g_IntersectionPointsContainer;
        extern std::vector<InterfaceAverageData> mInterfaceAverageContainer;
        
    }
}

#endif // KRATOS_INTERSECTION_POINTS_CONTAINER_H_INCLUDED
