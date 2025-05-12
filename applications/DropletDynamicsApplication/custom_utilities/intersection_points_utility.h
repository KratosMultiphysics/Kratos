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
    
    static void ExtractIntersectionPointsFromSplitter(DivideGeometry<Node>* p_splitter, int elementId);
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
    
    struct ConicCoefficients
    {
        double a, b, c, d, e;  // y² + ax² + bxy + cy + dx + e = 0
        ConicCoefficients() : a(0.0), b(0.0), c(0.0), d(0.0), e(0.0) {}
    };
    
    // Process all intersection points and fit curves using element connectivity
    static void ProcessIntersectionPointsAndFitCurves(
        const std::vector<IntersectionPointData>& points,
        const ModelPart& modelPart,
        const std::string& output_file = "element_curves.txt");

    static void ProcessIntersectionPointsAndFitCurvesparabola(
        const std::vector<IntersectionPointData>& points,
        const ModelPart& modelPart,
        const std::string& output_file = "element_curves.txt");
        
    static void ProcessIntersectionPointsAndFitGeneralConic(
        const std::vector<IntersectionPointData>& points,
        const ModelPart& modelPart,
        const std::string& output_file = "element_curves_general_conic.txt");
        
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
        
    // Fit general conic section curves to points
    static void FitGeneralConicCurves(
        const std::map<int, std::vector<IntersectionPointData>>& curves,
        std::map<int, ConicCoefficients>& fits);
        
    // Save results to file
    static void SaveToFile(
        const std::map<int, std::vector<IntersectionPointData>>& curves,
        const std::map<int, QuadraticCoefficients>& fits,
        const std::string& filename);
        
    // Save results for general conic fits to file
    static void SaveConicFitsToFile(
        const std::map<int, std::vector<IntersectionPointData>>& curves,
        const std::map<int, ConicCoefficients>& fits,
        const std::string& filename);
};

/// Utility for computing interface averages
class InterfaceAveragesUtility
{
public:
    typedef Matrix MatrixType;
    
    /// Compute interface averages for a given element
    static void CollectElementInterfaceAverages(Element::Pointer pElement);
    
    /// Compute interface averages for all elements in a model part
    static void ComputeModelPartInterfaceAverages(const ModelPart& rModelPart);
    
    /// Clear all stored interface average data
    static void ClearInterfaceAverages();
    
    /// Get the global container of interface average data
    static const std::vector<InterfaceAverageData>& GetInterfaceAverages();
    
    /// Apply interface averages to model part elements as element variables
    static void ApplyInterfaceAveragesToModelPart(ModelPart& rModelPart, const std::string& variable_name);

private:
    /// Global container for interface average data
    static std::vector<InterfaceAverageData> mInterfaceAverageContainer;
};

// Functions for combined intersection length and normal data
void ClearIntersectionDataWithNormal();
const std::vector<IntersectionDataWithNormal>& GetIntersectionDataWithNormal();
int CollectIntersectionDataWithNormal(ModelPart& rModelPart);
void SaveIntersectionDataWithNormalToFile(const std::string& filename);
int SaveIntersectionLengthsAndNormalsToFile(const ModelPart& rModelPart, const std::string& filename);

// Helper functions for element variables
double GetElementIntersectionLength(const Element& rElement);
double GetElementCutNormalX(const Element& rElement);
double GetElementCutNormalY(const Element& rElement);
double GetElementCutNormalZ(const Element& rElement);

// Calculate and store intersection lengths
int CalculateAndStoreElementIntersectionLengths(ModelPart& rModelPart);
void SaveIntersectionLengthsToFile(const std::map<int, double>& rIntersectionLengths, const std::string& filename);

// Element cut normal functions
int SetElementCutNormals(ModelPart& rModelPart);
void SaveElementAverageNormalsToFile(const ModelPart& rModelPart, const std::string& filename);

// Helper functions finding elements and neighbors
const InterfaceAverageData* FindElementInInterfaceAverages(
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    int ElementId);

std::vector<int> GetElementNeighbors(const ModelPart& rModelPart, int ElementId);

// Normal fitting functions
array_1d<double, 3> FitLinearNormal(
    const ModelPart& rModelPart,
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    int ElementId,
    double& a0, double& a1, double& a2,
    double& b0, double& b1, double& b2);

// Save functions
void SaveFittedNormalsToFile(
    const ModelPart& rModelPart,
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    const std::string& Filename);

// Averaged normal utility functions
void ClearAveragedNormals(
    ModelPart& rModelPart,
    const std::string& VariableName = "ELEMENT_CUT_NORMAL_AVERAGED");

int ComputeAndStoreAveragedNormals(
    ModelPart& rModelPart,
    int NeighborLevels = 1,
    const std::string& VariableName = "ELEMENT_CUT_NORMAL_AVERAGED");

array_1d<double, 3> ComputeAveragedElementNormal(
    const ModelPart& rModelPart,
    int ElementId,
    int NeighborLevels = 1);

// Save averaged normals to file
void SaveAveragedNormalsToFile(
    const ModelPart& rModelPart,
    const std::string& Filename,
    const std::string& VariableName = "ELEMENT_CUT_NORMAL_AVERAGED");

}
}
#endif // KRATOS_INTERSECTION_POINTS_UTILITY_H_INCLUDED