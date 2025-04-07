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

// Forward declaration of ElementBasedCurveFitter
class ElementBasedCurveFitter;

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
    
    // Add this function declaration inside the IntersectionPointsUtility class
    static void ProcessIntersectionPointsAndFitCurves(const std::string& output_file = "element_curves.txt");
    //////////////////////////////////////
};

class ElementBasedCurveFitter
{
public:
struct ElementIntersection
{
    int elementId;
    std::vector<IntersectionPointData> points;
    
    // Add a default constructor
    ElementIntersection() : elementId(-1) {}
    
    // Keep the existing constructor
    ElementIntersection(int id) : elementId(id) {}
};
    
    struct QuadraticCoefficients
    {
        double a, b, c;  // y = ax² + bx + c
        QuadraticCoefficients() : a(0.0), b(0.0), c(0.0) {}
    };
    
    // Process all intersection points and fit curves using element connectivity
    static void ProcessIntersectionPointsAndFitCurves(
        const std::vector<IntersectionPointData>& points,
        const ModelPart& modelPart,
        const std::string& output_file = "element_curves.txt");
        
private:
    // Group points by element
    static void GroupPointsByElement(
        const std::vector<IntersectionPointData>& all_points,
        std::map<int, ElementIntersection>& element_intersections);
        
    // Find neighboring elements that share intersection points
    static void FindConnectedElements(
        const std::map<int, ElementIntersection>& element_intersections,
        const ModelPart& modelPart,
        std::map<int, std::vector<int>>& element_neighbors);
        
    // Build continuous curves using element connectivity
    static void BuildContinuousCurves(
        const std::map<int, ElementIntersection>& element_intersections,
        const std::map<int, std::vector<int>>& element_neighbors,
        std::map<int, std::vector<IntersectionPointData>>& curves);
        
    // Fit quadratic curves to points
    static void FitQuadraticCurves(
        const std::map<int, std::vector<IntersectionPointData>>& curves,
        std::map<int, QuadraticCoefficients>& fits);
        
    // Save results to file
    static void SaveToFile(
        const std::map<int, std::vector<IntersectionPointData>>& curves,
        const std::map<int, QuadraticCoefficients>& fits,
        const std::string& filename);
};

}
}
#endif // KRATOS_INTERSECTION_POINTS_UTILITY_H_INCLUDED