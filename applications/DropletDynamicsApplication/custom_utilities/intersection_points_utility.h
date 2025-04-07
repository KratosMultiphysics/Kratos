//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alireza
//
//

// intersection_points_utility.h
#ifndef KRATOS_INTERSECTION_POINTS_UTILITY_H_INCLUDED
#define KRATOS_INTERSECTION_POINTS_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/divide_geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "intersection_points_container.h"

namespace Kratos
{
namespace KratosDropletDynamics
{
    class IntersectionPointsUtility
    {
    public:
        // Collect intersection points from a specific element
        static void CollectElementIntersectionPoints(Element::Pointer pElement);
        
        // Clear all intersection points
        static void ClearIntersectionPoints();
        
        // Get reference to all intersection points
        static const std::vector<IntersectionPointData>& GetIntersectionPoints();
        
        // Save all intersection points to a file
        static void SaveIntersectionPointsToFile(const std::string& filename);

        //////////////////////////////////////
        // static void AddIntersectionPoint(int elementId, int pointId, const array_1d<double, 3>& coordinates);
        static void DiagnosticOutput(const ModelPart& rModelPart);
        static void ExtractIntersectionPointsFromSplitter(DivideGeometry<Node>* p_splitter, int elementId);
        //////////////////////////////////////
    };
}
}

#endif // KRATOS_INTERSECTION_POINTS_UTILITY_H_INCLUDED
