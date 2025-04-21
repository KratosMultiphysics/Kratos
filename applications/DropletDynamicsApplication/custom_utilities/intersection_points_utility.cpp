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

// intersection_points_utility.cpp
#include "intersection_points_utility.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "utilities/divide_geometry.h"
#include "droplet_dynamics_application_variables.h"
#include <fstream>
#include "custom_elements/droplet_dynamics_element.h"
#include "droplet_dynamics_application_variables.h"
#include "../../FluidDynamicsApplication/custom_utilities/two_fluid_navier_stokes_data.h"

namespace Kratos
{
namespace KratosDropletDynamics
{
    // Define the global container for intersection points
    std::vector<IntersectionPointData> g_IntersectionPointsContainer;

void IntersectionPointsUtility::CollectElementIntersectionPoints(Element::Pointer pElement)
{
    // Get the geometry and distance values from the element
    auto p_geom = pElement->pGetGeometry();
    
    // Only proceed if the element is properly initialized
    if (!p_geom) return;
    
    // Get the distance values from the element's nodes
    Vector nodal_distances;
    nodal_distances.resize(p_geom->size());
    
    for (unsigned int i = 0; i < p_geom->size(); ++i) {
        nodal_distances[i] = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE);
    }

    // Check if the element is actually split by the interface
    bool is_split = false;
    const double sign_threshold = 1e-14;
    int pos_count = 0, neg_count = 0;
    for (unsigned int i = 0; i < p_geom->size(); ++i) {
        if (nodal_distances[i] > sign_threshold) {
            pos_count++;
        } else if (nodal_distances[i] < -sign_threshold) {
            neg_count++;
        }
    }
    
    // Element is split only if it has both positive and negative distance values
    is_split = (pos_count > 0 && neg_count > 0);
    
    // Only proceed with intersection calculations if the element is actually split
    if (!is_split) {
        return; // Skip this element as it's not split by the interface
    }
    
    // Structure nodes info
    Vector structure_node_id = ZeroVector(p_geom->size());
    // for (unsigned int i_node = 0; i_node < p_geom->size(); i_node++) {
    //     if ((*p_geom)[i_node].Is(BOUNDARY)) {
    //         structure_node_id[i_node] = 1.0;
    //     }
    // }
    
    // Create the modified shape functions utility
    ModifiedShapeFunctions::Pointer p_modified_sh_func;
    
    // Create the appropriate modified shape functions based on geometry type
    if (p_geom->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
        p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, nodal_distances, structure_node_id);
    } 
    else if (p_geom->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
        p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, nodal_distances, structure_node_id);
    }
    
    if (p_modified_sh_func) {
        // Get the splitting utility
        auto p_splitting_util = p_modified_sh_func->pGetSplittingUtil();
        
        if (p_splitting_util) {
            try {
                // Force generation of the intersection skin
                if (p_geom->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
                    auto p_triangle_splitter = dynamic_cast<DivideTriangle2D3<Node>*>(p_splitting_util.get());
                    if (p_triangle_splitter) {
                        p_triangle_splitter->GenerateIntersectionsSkin();
                        
                        // Extract real intersection points
                        ExtractIntersectionPointsFromSplitter(p_triangle_splitter, pElement->Id());
                    }
                }
                else if (p_geom->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
                    auto p_tetra_splitter = dynamic_cast<DivideTetrahedra3D4<Node>*>(p_splitting_util.get());
                    if (p_tetra_splitter) {
                        p_tetra_splitter->GenerateIntersectionsSkin();
                        
                        // Extract real intersection points
                        ExtractIntersectionPointsFromSplitter(p_tetra_splitter, pElement->Id());
                    }
                }
                
                std::cout << "Processed interface points for element " << pElement->Id() << std::endl;
            }
            catch (std::exception& e) {
                std::cerr << "Error processing element " << pElement->Id() 
                          << ": " << e.what() << std::endl;
            }
        }
    }
}
    
    void IntersectionPointsUtility::ClearIntersectionPoints()
    {
        g_IntersectionPointsContainer.clear();
    }
    
    const std::vector<IntersectionPointData>& IntersectionPointsUtility::GetIntersectionPoints()
    {
        return g_IntersectionPointsContainer;
    }
    
void IntersectionPointsUtility::SaveIntersectionPointsToFile(const std::string& filename)
{
    std::cout << "Saving " << g_IntersectionPointsContainer.size() << " intersection points to file: " << filename << std::endl;
    
    // Rest of your original function...
    std::ofstream outFile(filename);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    outFile << std::fixed << std::setprecision(15);  
    outFile << "Element_ID\tPoint_ID\tX\tY\tZ" << std::endl;
    
    for (const auto& point : g_IntersectionPointsContainer) {
        outFile << point.elementId << "\t" 
               << point.pointId << "\t"
               << point.coordinates[0] << "\t" 
               << point.coordinates[1] << "\t" 
               << point.coordinates[2] << std::endl;
    }
    
    outFile.close();
    std::cout << "Successfully wrote " << g_IntersectionPointsContainer.size() << " points to " << filename << std::endl;
}


void IntersectionPointsUtility::DiagnosticOutput(const ModelPart& rModelPart)
{
    int total_elements = rModelPart.NumberOfElements();
    int split_elements = 0;
    
    for (auto& element : rModelPart.Elements()) {
        auto p_geom = element.pGetGeometry();
        if (!p_geom) continue;
        
        // Get the distance values
        Vector nodal_distances;
        nodal_distances.resize(p_geom->size());
        for (unsigned int i = 0; i < p_geom->size(); ++i) {
            nodal_distances[i] = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE);
        }
        
        // Check if element is split
        const double sign_threshold = 1e-14;
        int pos_count = 0, neg_count = 0;
        for (unsigned int i = 0; i < p_geom->size(); ++i) {
            if (nodal_distances[i] > sign_threshold) {
                pos_count++;
            } else if (nodal_distances[i] < -sign_threshold) {
                neg_count++;
            }
        }
        
        bool is_split = (pos_count > 0 && neg_count > 0);
        
        if (is_split) {
            split_elements++;
        }
    }
    
    std::cout << "Diagnostic Output:" << std::endl;
    std::cout << "Total Elements: " << total_elements << std::endl;
    std::cout << "Split Elements: " << split_elements << std::endl;
    std::cout << "Points Collected: " << g_IntersectionPointsContainer.size() << std::endl;
}
void IntersectionPointsUtility::ExtractIntersectionPointsFromSplitter(DivideGeometry<Node>* p_splitter, int elementId)
{
    if (!p_splitter) return;
    
    try {
        // Try to get interface points using the GetInterfacePoints method
        auto interface_points = p_splitter->GetInterfacePoints();
        
        if(interface_points.size() > 0) {
            std::cout << "Found " << interface_points.size() << " interface points for element " << elementId << std::endl;
            
            // Add each interface point to our container
            int point_count = 0;
            for (size_t i = 0; i < interface_points.size(); ++i) {
                if (interface_points[i]) {  // Make sure the pointer is valid
                    IntersectionPointData point;
                    point.elementId = elementId;
                    point.pointId = point_count++;
                    
                    // Copy coordinates from the IndexedPoint
                    const auto& coords = interface_points[i]->Coordinates();
                    point.coordinates[0] = coords[0];
                    point.coordinates[1] = coords[1];
                    point.coordinates[2] = coords[2];
                    
                    g_IntersectionPointsContainer.push_back(point);
                }
            }
            
            std::cout << "Added " << point_count << " real intersection points from element " << elementId << std::endl;
        } else {
            std::cout << "No interface points found for element " << elementId << std::endl;
            

               }
        } catch (std::exception& e) {
        std::cerr << "Error extracting interface points: " << e.what() << std::endl;
        
       
        }
}

#include <iomanip>  // For std::setprecision

void IntersectionPointsUtility::ProcessIntersectionPointsAndFitCurves(const std::string& output_file)
{
    // Get all intersection points
    const auto& points = g_IntersectionPointsContainer;
    
    if (points.empty()) {
        std::cout << "No intersection points available for circle fitting." << std::endl;
        return;
    }
    
    // Configuration parameters
    const int MIN_POINTS_FOR_CIRCLE_FIT = 3;  // Absolute minimum needed for circle
    const int NEIGHBOR_EXPANSION_LEVEL = 4;   // Expand to n-hop neighbors
    
    std::cout << "Starting circle fitting with " << points.size() << " intersection points." << std::endl;
    std::cout << "Using all available points from 2-hop neighborhoods." << std::endl;
    
    // Group points by element
    std::map<int, std::vector<IntersectionPointData>> element_points;
    // Create a map of points by their coordinates
    std::map<std::pair<double, double>, std::vector<int>> point_to_elements;
    
    for (const auto& point : points) {
        int elemId = point.elementId;
        
        // Round coordinates to handle floating point precision
        double x = std::round(point.coordinates[0] * 1.0E15) / 1.0E15;
        double y = std::round(point.coordinates[1] * 1.0E15) / 1.0E15;
        std::pair<double, double> coord_key(x, y);
        
        // Add this element to the list for this point
        point_to_elements[coord_key].push_back(elemId);
        
        // Add this point to the element's list
        element_points[elemId].push_back(point);
    }
    
    // Find element neighbors (elements that share intersection points)
    std::map<int, std::set<int>> element_neighbors;
    
    for (const auto& [coord, elements] : point_to_elements) {
        // If this point belongs to multiple elements, they are neighbors
        for (size_t i = 0; i < elements.size(); ++i) {
            for (size_t j = i+1; j < elements.size(); ++j) {
                element_neighbors[elements[i]].insert(elements[j]);
                element_neighbors[elements[j]].insert(elements[i]);
            }
        }
    }
    
    // Expand the neighborhood to n-hop neighbors
    std::cout << "Expanding neighborhood with " << NEIGHBOR_EXPANSION_LEVEL << " hops..." << std::endl;
    std::map<int, std::set<int>> expanded_neighbors = element_neighbors;
    
    for (int hop = 2; hop <= NEIGHBOR_EXPANSION_LEVEL; hop++) {
        std::map<int, std::set<int>> next_level_neighbors = expanded_neighbors;
        
        for (const auto& [elemId, current_neighbors] : expanded_neighbors) {
            for (int neighbor : current_neighbors) {
                // for (int next_hop : expanded_neighbors[neighbor]) {
                for (int next_hop : element_neighbors[neighbor]) {
                    if (next_hop != elemId && !expanded_neighbors[elemId].count(next_hop)) {
                        next_level_neighbors[elemId].insert(next_hop);
                    }
                }
            }
        }
        
        expanded_neighbors = next_level_neighbors;
        std::cout << "Completed " << hop << "-hop neighborhood expansion." << std::endl;
    }
    
    // Structure to hold circle fit coefficients
    struct CircleCoefficients {
        double a;  // x-center
        double b;  // y-center
        double c;  // radius squared
    };
    
    // Maps to store results
    std::map<int, CircleCoefficients> elementFits;
    std::map<int, int> elementTotalPoints;
    
    // Process each element
    for (const auto& [elemId, neighbors] : expanded_neighbors) {
        // Get original points for this element
        std::vector<IntersectionPointData> original_points = element_points[elemId];
        int original_point_count = original_points.size();
        
        // Create a pool of neighbor points
        std::vector<IntersectionPointData> neighbor_points;
        for (int neighborId : neighbors) {
            neighbor_points.insert(neighbor_points.end(), 
                                 element_points[neighborId].begin(), 
                                 element_points[neighborId].end());
        }
        
        // Remove duplicates and points shared with original set
        std::map<std::pair<double, double>, IntersectionPointData> unique_neighbor_points;
        for (const auto& point : neighbor_points) {
            double x = std::round(point.coordinates[0] * 1.0E15) / 1.0E15;
            double y = std::round(point.coordinates[1] * 1.0E15) / 1.0E15;
            std::pair<double, double> key(x, y);
            
            // Skip points that are in the original set
            bool is_in_original = false;
            for (const auto& orig_point : original_points) {
                double ox = std::round(orig_point.coordinates[0] * 1.0E15) / 1.0E15;
                double oy = std::round(orig_point.coordinates[1] * 1.0E15) / 1.0E15;
                if (ox == x && oy == y) {
                    is_in_original = true;
                    break;
                }
            }
            
            if (!is_in_original) {
                unique_neighbor_points[key] = point;
            }
        }
        
        // Create a vector of unique neighbor points
        neighbor_points.clear();
        for (const auto& [_, point] : unique_neighbor_points) {
            neighbor_points.push_back(point);
        }
        
        // Build the set of points for circle fitting
        std::vector<IntersectionPointData> combined_points = original_points;
        
        // Add all neighbor points
        combined_points.insert(combined_points.end(), neighbor_points.begin(), neighbor_points.end());
        
        int points_from_neighbors = neighbor_points.size();
        
        // Only fit if we have enough points
        if (combined_points.size() >= MIN_POINTS_FOR_CIRCLE_FIT) {
            std::cout << "Element " << elemId 
                      << " has " << combined_points.size() 
                      << " points for circle fitting (" 
                      << original_point_count << " original + " 
                      << points_from_neighbors << " from neighbors)." << std::endl;
            
            // Prepare matrices for least squares fitting
            double sum_x = 0.0, sum_y = 0.0;
            double sum_x2 = 0.0, sum_y2 = 0.0;
            double sum_xy = 0.0;
            double sum_x2y2 = 0.0;  // sum of (x^2 + y^2)
            double sum_x3 = 0.0, sum_xy2 = 0.0;
            double sum_x2y = 0.0, sum_y3 = 0.0;
            
            for (const auto& point : combined_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                
                double x2 = x * x;
                double y2 = y * y;
                
                sum_x += x;
                sum_y += y;
                sum_x2 += x2;
                sum_y2 += y2;
                sum_xy += x * y;
                sum_x2y2 += (x2 + y2);
                sum_x3 += x * x2;
                sum_xy2 += x * y2;
                sum_x2y += x2 * y;
                sum_y3 += y * y2;
            }
            
            int n = combined_points.size();
            
            // Set up the system of equations: x^2 + y^2 = 2ax + 2by + d
            Matrix A(3, 3);
            Vector b(3);
            
            A(0, 0) = sum_x2;    A(0, 1) = sum_xy;    A(0, 2) = sum_x;
            A(1, 0) = sum_xy;    A(1, 1) = sum_y2;    A(1, 2) = sum_y;
            A(2, 0) = sum_x;     A(2, 1) = sum_y;     A(2, 2) = n;
            
            b[0] = sum_x3 + sum_xy2;  // sum of x * (x^2 + y^2)
            b[1] = sum_x2y + sum_y3;  // sum of y * (x^2 + y^2)
            b[2] = sum_x2y2;          // sum of (x^2 + y^2)
            
            // Solve using Cramer's rule
            double det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
                         A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                         A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));

            
            Matrix A1 = A, A2 = A, A3 = A;
            
            for (int i = 0; i < 3; i++) {
                A1(i, 0) = b[i];
                A2(i, 1) = b[i];
                A3(i, 2) = b[i];
            }
            
            double det1 = A1(0, 0) * (A1(1, 1) * A1(2, 2) - A1(2, 1) * A1(1, 2)) -
                          A1(0, 1) * (A1(1, 0) * A1(2, 2) - A1(1, 2) * A1(2, 0)) +
                          A1(0, 2) * (A1(1, 0) * A1(2, 1) - A1(1, 1) * A1(2, 0));
                          
            double det2 = A2(0, 0) * (A2(1, 1) * A2(2, 2) - A2(2, 1) * A2(1, 2)) -
                          A2(0, 1) * (A2(1, 0) * A2(2, 2) - A2(1, 2) * A2(2, 0)) +
                          A2(0, 2) * (A2(1, 0) * A2(2, 1) - A2(1, 1) * A2(2, 0));
                          
            double det3 = A3(0, 0) * (A3(1, 1) * A3(2, 2) - A3(2, 1) * A3(1, 2)) -
                          A3(0, 1) * (A3(1, 0) * A3(2, 2) - A3(1, 2) * A3(2, 0)) +
                          A3(0, 2) * (A3(1, 0) * A3(2, 1) - A3(1, 1) * A3(2, 0));
            
            // Solve for parameters
            double twoA = det1 / det;
            double twoB = det2 / det;
            double d = det3 / det;
            
            CircleCoefficients fit;
            fit.a = twoA / 2.0;  // x-center
            fit.b = twoB / 2.0;  // y-center
            fit.c = d + fit.a * fit.a + fit.b * fit.b;  // radius squared
            
            // Store results
            elementFits[elemId] = fit;
            elementTotalPoints[elemId] = combined_points.size();
            
            // Calculate error on original points
            double radius = std::sqrt(fit.c);
            double total_error = 0.0;
            
            for (const auto& point : original_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                double dist_squared = (x - fit.a) * (x - fit.a) + (y - fit.b) * (y - fit.b);
                total_error += std::abs(dist_squared - fit.c);
            }
            
            double avg_error = total_error / (original_points.empty() ? 1.0 : original_points.size());
            
            std::cout << "Element " << elemId 
                      << " fitted with " << combined_points.size() 
                      << " points: (x-" << fit.a 
                      << ")² + (y-" << fit.b << ")² = " << fit.c 
                      << " (radius = " << radius << ")" 
                      << std::endl;
            std::cout << "    Average fit error on original points: " << avg_error << std::endl;
        } else {
            std::cout << "Element " << elemId 
                      << " has only " << combined_points.size() 
                      << " unique points (less than " << MIN_POINTS_FOR_CIRCLE_FIT 
                      << " required) - cannot perform circle fitting." << std::endl;
        }
    }
    
    // Write results to file
    std::ofstream outFile(output_file);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << " for writing." << std::endl;
        return;
    }
    
    outFile << "Element_ID\tNum_Original_Points\tTotal_Points\ta(center_x)\tb(center_y)\tc(radius_squared)\tradius\tavg_error\n";
    
    for (const auto& [elemId, fit] : elementFits) {
        int numPoints = element_points[elemId].size();
        int totalPoints = elementTotalPoints[elemId];
        double radius = std::sqrt(fit.c);
        
        // Calculate error
        double total_error = 0.0;
        for (const auto& point : element_points[elemId]) {
            double x = point.coordinates[0];
            double y = point.coordinates[1];
            double dist_squared = (x - fit.a) * (x - fit.a) + (y - fit.b) * (y - fit.b);
            total_error += std::abs(dist_squared - fit.c);
        }
        
        double avg_error = total_error / (numPoints > 0 ? numPoints : 1.0);
        
        outFile << elemId << "\t" 
                << numPoints << "\t"
                << totalPoints << "\t"
                << fit.a << "\t" 
                << fit.b << "\t" 
                << fit.c << "\t"
                << radius << "\t"
                << avg_error << "\n";
    }
    
    outFile.close();
    
    std::cout << "Saved " << elementFits.size() << " element circle fits to " << output_file << std::endl;
    std::cout << "Used all available points from 2-hop neighborhoods." << std::endl;
}

void IntersectionPointsUtility::ProcessIntersectionPointsAndFitCurvesparabola(const std::string& output_file)
{
    // Get all intersection points
    const auto& points = g_IntersectionPointsContainer;
    
    if (points.empty()) {
        std::cout << "No intersection points available for curve fitting." << std::endl;
        return;
    }
    
    // Configuration parameters
    const int MIN_POINTS_FOR_CURVE_FIT = 3;  // Absolute minimum needed for quadratic
    const int NEIGHBOR_EXPANSION_LEVEL = 4;   // Expand to n-hop neighbors
    
    std::cout << "Starting quadratic curve fitting with " << points.size() << " intersection points." << std::endl;
    std::cout << "Using all available points from 2-hop neighborhoods." << std::endl;
    
    // Group points by element
    std::map<int, std::vector<IntersectionPointData>> element_points;
    // Create a map of points by their coordinates
    std::map<std::pair<double, double>, std::vector<int>> point_to_elements;
    
    for (const auto& point : points) {
        int elemId = point.elementId;
        
        // Round coordinates to handle floating point precision
        double x = std::round(point.coordinates[0] * 1.0E15) / 1.0E15;
        double y = std::round(point.coordinates[1] * 1.0E15) / 1.0E15;
        std::pair<double, double> coord_key(x, y);
        
        // Add this element to the list for this point
        point_to_elements[coord_key].push_back(elemId);
        
        // Add this point to the element's list
        element_points[elemId].push_back(point);
    }
    
    // Find element neighbors (elements that share intersection points)
    std::map<int, std::set<int>> element_neighbors;
    
    for (const auto& [coord, elements] : point_to_elements) {
        // If this point belongs to multiple elements, they are neighbors
        for (size_t i = 0; i < elements.size(); ++i) {
            for (size_t j = i+1; j < elements.size(); ++j) {
                element_neighbors[elements[i]].insert(elements[j]);
                element_neighbors[elements[j]].insert(elements[i]);
            }
        }
    }
    
    // Expand the neighborhood to n-hop neighbors
    std::cout << "Expanding neighborhood with " << NEIGHBOR_EXPANSION_LEVEL << " hops..." << std::endl;
    std::map<int, std::set<int>> expanded_neighbors = element_neighbors;
    
    for (int hop = 2; hop <= NEIGHBOR_EXPANSION_LEVEL; hop++) {
        std::map<int, std::set<int>> next_level_neighbors = expanded_neighbors;
        
        for (const auto& [elemId, current_neighbors] : expanded_neighbors) {
            for (int neighbor : current_neighbors) {
                // for (int next_hop : expanded_neighbors[neighbor]) {
                for (int next_hop : element_neighbors[neighbor]) {
                    if (next_hop != elemId && !expanded_neighbors[elemId].count(next_hop)) {
                        next_level_neighbors[elemId].insert(next_hop);
                    }
                }
            }
        }
        
        expanded_neighbors = next_level_neighbors;
        std::cout << "Completed " << hop << "-hop neighborhood expansion." << std::endl;
    }
    
    // Structure to hold quadratic curve fit coefficients (y = ax² + bx + c)
    struct QuadraticCoefficients {
        double a;  // coefficient of x²
        double b;  // coefficient of x
        double c;  // constant term
    };
    
    // Maps to store results
    std::map<int, QuadraticCoefficients> elementFits;
    std::map<int, int> elementTotalPoints;
    
    // Process each element
    for (const auto& [elemId, neighbors] : expanded_neighbors) {
        // Get original points for this element
        std::vector<IntersectionPointData> original_points = element_points[elemId];
        int original_point_count = original_points.size();
        
        // Create a pool of neighbor points
        std::vector<IntersectionPointData> neighbor_points;
        for (int neighborId : neighbors) {
            neighbor_points.insert(neighbor_points.end(), 
                                 element_points[neighborId].begin(), 
                                 element_points[neighborId].end());
        }
        
        // Remove duplicates and points shared with original set
        std::map<std::pair<double, double>, IntersectionPointData> unique_neighbor_points;
        for (const auto& point : neighbor_points) {
            double x = std::round(point.coordinates[0] * 1.0E15) / 1.0E15;
            double y = std::round(point.coordinates[1] * 1.0E15) / 1.0E15;
            std::pair<double, double> key(x, y);
            
            // Skip points that are in the original set
            bool is_in_original = false;
            for (const auto& orig_point : original_points) {
                double ox = std::round(orig_point.coordinates[0] * 1.0E15) / 1.0E15;
                double oy = std::round(orig_point.coordinates[1] * 1.0E15) / 1.0E15;
                if (ox == x && oy == y) {
                    is_in_original = true;
                    break;
                }
            }
            
            if (!is_in_original) {
                unique_neighbor_points[key] = point;
            }
        }
        
        // Create a vector of unique neighbor points
        neighbor_points.clear();
        for (const auto& [_, point] : unique_neighbor_points) {
            neighbor_points.push_back(point);
        }
        
        // Build the set of points for curve fitting
        std::vector<IntersectionPointData> combined_points = original_points;
        
        // Add all neighbor points
        combined_points.insert(combined_points.end(), neighbor_points.begin(), neighbor_points.end());
        
        int points_from_neighbors = neighbor_points.size();
        
        // Only fit if we have enough points
        if (combined_points.size() >= MIN_POINTS_FOR_CURVE_FIT) {
            std::cout << "Element " << elemId 
                      << " has " << combined_points.size() 
                      << " points for quadratic fitting (" 
                      << original_point_count << " original + " 
                      << points_from_neighbors << " from neighbors)." << std::endl;
            
            // Prepare matrices for least squares fitting of y = ax² + bx + c
            double sum_x = 0.0, sum_y = 0.0;
            double sum_x2 = 0.0, sum_x3 = 0.0, sum_x4 = 0.0;
            double sum_xy = 0.0, sum_x2y = 0.0;
            
            for (const auto& point : combined_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                
                double x2 = x * x;
                
                sum_x += x;
                sum_y += y;
                sum_x2 += x2;
                sum_x3 += x2 * x;
                sum_x4 += x2 * x2;
                sum_xy += x * y;
                sum_x2y += x2 * y;
            }
            
            int n = combined_points.size();
            
            // Set up the system of equations for quadratic fit: y = ax² + bx + c
            Matrix A(3, 3);
            Vector b(3);
            
            A(0, 0) = sum_x4;    A(0, 1) = sum_x3;    A(0, 2) = sum_x2;
            A(1, 0) = sum_x3;    A(1, 1) = sum_x2;    A(1, 2) = sum_x;
            A(2, 0) = sum_x2;    A(2, 1) = sum_x;     A(2, 2) = n;
            
            b[0] = sum_x2y;
            b[1] = sum_xy;
            b[2] = sum_y;
            
            // Solve using Cramer's rule
            double det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
                         A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                         A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));         
            Matrix A1 = A, A2 = A, A3 = A;
            
            for (int i = 0; i < 3; i++) {
                A1(i, 0) = b[i];
                A2(i, 1) = b[i];
                A3(i, 2) = b[i];
            }
            
            double det1 = A1(0, 0) * (A1(1, 1) * A1(2, 2) - A1(2, 1) * A1(1, 2)) -
                          A1(0, 1) * (A1(1, 0) * A1(2, 2) - A1(1, 2) * A1(2, 0)) +
                          A1(0, 2) * (A1(1, 0) * A1(2, 1) - A1(1, 1) * A1(2, 0));
                          
            double det2 = A2(0, 0) * (A2(1, 1) * A2(2, 2) - A2(2, 1) * A2(1, 2)) -
                          A2(0, 1) * (A2(1, 0) * A2(2, 2) - A2(1, 2) * A2(2, 0)) +
                          A2(0, 2) * (A2(1, 0) * A2(2, 1) - A2(1, 1) * A2(2, 0));
                          
            double det3 = A3(0, 0) * (A3(1, 1) * A3(2, 2) - A3(2, 1) * A3(1, 2)) -
                          A3(0, 1) * (A3(1, 0) * A3(2, 2) - A3(1, 2) * A3(2, 0)) +
                          A3(0, 2) * (A3(1, 0) * A3(2, 1) - A3(1, 1) * A3(2, 0));
            
            // Solve for quadratic parameters
            QuadraticCoefficients fit;
            fit.a = det1 / det;  // coefficient of x²
            fit.b = det2 / det;  // coefficient of x
            fit.c = det3 / det;  // constant term
            
            // Store results
            elementFits[elemId] = fit;
            elementTotalPoints[elemId] = combined_points.size();
            
            // Calculate error on original points
            double total_error = 0.0;
            
            for (const auto& point : original_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                
                // Calculate y value from fitted curve
                double fitted_y = fit.a * x * x + fit.b * x + fit.c;
                
                // Error is the vertical distance between point and curve
                double error = std::abs(y - fitted_y);
                total_error += error;
            }
            
            double avg_error = total_error / (original_points.empty() ? 1.0 : original_points.size());
            
            std::cout << "Element " << elemId 
                      << " fitted with " << combined_points.size() 
                      << " points: y = " << fit.a
                      << "x² + " << fit.b << "x + " << fit.c
                      << std::endl;
            std::cout << "    Average fit error on original points: " << avg_error << std::endl;
        } else {
            std::cout << "Element " << elemId 
                      << " has only " << combined_points.size() 
                      << " unique points (less than " << MIN_POINTS_FOR_CURVE_FIT 
                      << " required) - cannot perform quadratic curve fitting." << std::endl;
        }
    }
    
    // Write results to file
    std::ofstream outFile(output_file);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << " for writing." << std::endl;
        return;
    }
    
    outFile << "Element_ID\tNum_Original_Points\tTotal_Points\ta(x²)\tb(x)\tc(const)\tavg_error\n";
    
    for (const auto& [elemId, fit] : elementFits) {
        int numPoints = element_points[elemId].size();
        int totalPoints = elementTotalPoints[elemId];
        
        // Calculate error
        double total_error = 0.0;
        for (const auto& point : element_points[elemId]) {
            double x = point.coordinates[0];
            double y = point.coordinates[1];
            
            // Calculate y value from fitted curve
            double fitted_y = fit.a * x * x + fit.b * x + fit.c;
            
            // Error is the vertical distance
            double error = std::abs(y - fitted_y);
            total_error += error;
        }
        
        double avg_error = total_error / (numPoints > 0 ? numPoints : 1.0);
        
        outFile << elemId << "\t" 
                << numPoints << "\t"
                << totalPoints << "\t"
                << fit.a << "\t" 
                << fit.b << "\t" 
                << fit.c << "\t"
                << avg_error << "\n";
    }
    
    outFile.close();
    
    std::cout << "Saved " << elementFits.size() << " element quadratic curve fits to " << output_file << std::endl;
    std::cout << "Used all available points from 2-hop neighborhoods." << std::endl;
}

void IntersectionPointsUtility::ProcessIntersectionPointsAndFitGeneralConic(const std::string& output_file)
{
    // Get all intersection points
    const auto& points = g_IntersectionPointsContainer;
    
    if (points.empty()) {
        std::cout << "No intersection points available for curve fitting." << std::endl;
        return;
    }
    
    // Configuration parameters
    const int MIN_POINTS_FOR_CURVE_FIT = 5;  // Minimum needed for general conic section
    const int NEIGHBOR_EXPANSION_LEVEL = 6;   // Expand to n-hop neighbors
    
    std::cout << "Starting general conic section fitting with " << points.size() << " intersection points." << std::endl;
    std::cout << "Using all available points from 2-hop neighborhoods." << std::endl;
    
    // Group points by element
    std::map<int, std::vector<IntersectionPointData>> element_points;
    // Create a map of points by their coordinates
    std::map<std::pair<double, double>, std::vector<int>> point_to_elements;
    
    for (const auto& point : points) {
        int elemId = point.elementId;
        
        // Round coordinates to handle floating point precision
        double x = std::round(point.coordinates[0] * 1.0E15) / 1.0E15;
        double y = std::round(point.coordinates[1] * 1.0E15) / 1.0E15;
        std::pair<double, double> coord_key(x, y);
        
        // Add this element to the list for this point
        point_to_elements[coord_key].push_back(elemId);
        
        // Add this point to the element's list
        element_points[elemId].push_back(point);
    }
    
    // Find element neighbors (elements that share intersection points)
    std::map<int, std::set<int>> element_neighbors;
    
    for (const auto& [coord, elements] : point_to_elements) {
        // If this point belongs to multiple elements, they are neighbors
        for (size_t i = 0; i < elements.size(); ++i) {
            for (size_t j = i+1; j < elements.size(); ++j) {
                element_neighbors[elements[i]].insert(elements[j]);
                element_neighbors[elements[j]].insert(elements[i]);
            }
        }
    }
    
    // Expand the neighborhood to n-hop neighbors
    std::cout << "Expanding neighborhood with " << NEIGHBOR_EXPANSION_LEVEL << " hops..." << std::endl;
    std::map<int, std::set<int>> expanded_neighbors = element_neighbors;
    
    for (int hop = 2; hop <= NEIGHBOR_EXPANSION_LEVEL; hop++) {
        std::map<int, std::set<int>> next_level_neighbors = expanded_neighbors;
        
        for (const auto& [elemId, current_neighbors] : expanded_neighbors) {
            for (int neighbor : current_neighbors) {
                // for (int next_hop : expanded_neighbors[neighbor]) {
                for (int next_hop : element_neighbors[neighbor]) {
                    if (next_hop != elemId && !expanded_neighbors[elemId].count(next_hop)) {
                        next_level_neighbors[elemId].insert(next_hop);
                    }
                }
            }
        }
        
        expanded_neighbors = next_level_neighbors;
        std::cout << "Completed " << hop << "-hop neighborhood expansion." << std::endl;
    }
    
    // Structure to hold general conic section coefficients (y² + ax² + bxy + cy + dx + e = 0)
    struct ConicCoefficients {
        double a;  // coefficient of x²
        double b;  // coefficient of xy
        double c;  // coefficient of y
        double d;  // coefficient of x
        double e;  // constant term
    };
    
    // Maps to store results
    std::map<int, ConicCoefficients> elementFits;
    std::map<int, int> elementTotalPoints;
    
    // Process each element
    for (const auto& [elemId, neighbors] : expanded_neighbors) {
        // Get original points for this element
        std::vector<IntersectionPointData> original_points = element_points[elemId];
        int original_point_count = original_points.size();
        
        // Create a pool of neighbor points
        std::vector<IntersectionPointData> neighbor_points;
        for (int neighborId : neighbors) {
            neighbor_points.insert(neighbor_points.end(), 
                                 element_points[neighborId].begin(), 
                                 element_points[neighborId].end());
        }
        
        // Remove duplicates and points shared with original set
        std::map<std::pair<double, double>, IntersectionPointData> unique_neighbor_points;
        for (const auto& point : neighbor_points) {
            double x = std::round(point.coordinates[0] * 1.0E15) / 1.0E15;
            double y = std::round(point.coordinates[1] * 1.0E15) / 1.0E15;
            std::pair<double, double> key(x, y);
            
            // Skip points that are in the original set
            bool is_in_original = false;
            for (const auto& orig_point : original_points) {
                double ox = std::round(orig_point.coordinates[0] * 1.0E15) / 1.0E15;
                double oy = std::round(orig_point.coordinates[1] * 1.0E15) / 1.0E15;
                if (ox == x && oy == y) {
                    is_in_original = true;
                    break;
                }
            }
            
            if (!is_in_original) {
                unique_neighbor_points[key] = point;
            }
        }
        
        // Create a vector of unique neighbor points
        neighbor_points.clear();
        for (const auto& [_, point] : unique_neighbor_points) {
            neighbor_points.push_back(point);
        }
        
        // Build the set of points for curve fitting
        std::vector<IntersectionPointData> combined_points = original_points;
        
        // Add all neighbor points
        combined_points.insert(combined_points.end(), neighbor_points.begin(), neighbor_points.end());
        
        int points_from_neighbors = neighbor_points.size();
        
        // Only fit if we have enough points
        if (combined_points.size() >= MIN_POINTS_FOR_CURVE_FIT) {
            std::cout << "Element " << elemId 
                      << " has " << combined_points.size() 
                      << " points for general conic fitting (" 
                      << original_point_count << " original + " 
                      << points_from_neighbors << " from neighbors)." << std::endl;
            
            // Prepare matrices for least squares fitting of y² + ax² + bxy + cy + dx + e = 0
            // We'll use Matrix class (as in the original function)
            
            // For a general conic, we need to solve for 5 parameters (a, b, c, d, e)
            // We'll rearrange as: y² = -ax² - bxy - cy - dx - e
            Matrix A(5, 5);
            Vector b(5);
            
            // Initialize matrices with zeros
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    A(i, j) = 0.0;
                }
                b[i] = 0.0;
            }
            
            // Fill matrices by summing contributions from each point
            for (const auto& point : combined_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                
                double x2 = x * x;
                double y2 = y * y;
                double xy = x * y;
                
                // Row 0: x^4 x^3y x^2y x^3 x^2 | x^2y^2
                A(0, 0) += x2 * x2;       // x^4
                A(0, 1) += x2 * x * y;     // x^3y
                A(0, 2) += x2 * y;         // x^2y
                A(0, 3) += x2 * x;         // x^3
                A(0, 4) += x2;             // x^2
                b[0] += x2 * y2;           // x^2y^2
                
                // Row 1: x^3y x^2y^2 xy^2 x^2y xy | xy^3
                A(1, 0) += x2 * x * y;     // x^3y
                A(1, 1) += x2 * y2;        // x^2y^2
                A(1, 2) += x * y2;         // xy^2
                A(1, 3) += x2 * y;         // x^2y
                A(1, 4) += x * y;          // xy
                b[1] += x * y2 * y;        // xy^3
                
                // Row 2: x^2y xy^2 y^2 xy y | y^3
                A(2, 0) += x2 * y;         // x^2y
                A(2, 1) += x * y2;         // xy^2
                A(2, 2) += y2;             // y^2
                A(2, 3) += x * y;          // xy
                A(2, 4) += y;              // y
                b[2] += y2 * y;            // y^3
                
                // Row 3: x^3 x^2y xy x^2 x | xy^2
                A(3, 0) += x2 * x;         // x^3
                A(3, 1) += x2 * y;         // x^2y
                A(3, 2) += x * y;          // xy
                A(3, 3) += x2;             // x^2
                A(3, 4) += x;              // x
                b[3] += x * y2;            // xy^2
                
                // Row 4: x^2 xy y x 1 | y^2
                A(4, 0) += x2;             // x^2
                A(4, 1) += x * y;          // xy
                A(4, 2) += y;              // y
                A(4, 3) += x;              // x
                A(4, 4) += 1.0;            // 1
                b[4] += y2;                // y^2
            }
            
            // Solve the system for the conic coefficients
            // For simplicity, we'll use Cramer's rule as in the original code
            
            double det = determinant5x5(A);
            
            // Create copies of A for solving using Cramer's rule
            Matrix A1 = A, A2 = A, A3 = A, A4 = A, A5 = A;
            
            // Replace columns with b vector
            for (int i = 0; i < 5; i++) {
                A1(i, 0) = b[i];
                A2(i, 1) = b[i];
                A3(i, 2) = b[i];
                A4(i, 3) = b[i];
                A5(i, 4) = b[i];
            }
            
            // Calculate determinants
            double det1 = determinant5x5(A1);
            double det2 = determinant5x5(A2);
            double det3 = determinant5x5(A3);
            double det4 = determinant5x5(A4);
            double det5 = determinant5x5(A5);
            
            // Solve for conic parameters
            ConicCoefficients fit;
            fit.a = -det1 / det;  // coefficient of x²
            fit.b = -det2 / det;  // coefficient of xy
            fit.c = -det3 / det;  // coefficient of y
            fit.d = -det4 / det;  // coefficient of x
            fit.e = -det5 / det;  // constant term
            
            // Store results
            elementFits[elemId] = fit;
            elementTotalPoints[elemId] = combined_points.size();
            
            // Calculate error on original points
            double total_error = 0.0;
            
            for (const auto& point : original_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                
                // Calculate error as deviation from the conic equation: y² + ax² + bxy + cy + dx + e = 0
                double equation_value = y*y + fit.a*x*x + fit.b*x*y + fit.c*y + fit.d*x + fit.e;
                
                // Error is the absolute value of the equation
                double error = std::abs(equation_value);
                total_error += error;
            }
            
            double avg_error = total_error / (original_points.empty() ? 1.0 : original_points.size());
            
            std::cout << "Element " << elemId 
                      << " fitted with " << combined_points.size() 
                      << " points: y² + " << fit.a << "x² + " 
                      << fit.b << "xy + " << fit.c << "y + " 
                      << fit.d << "x + " << fit.e << " = 0"
                      << std::endl;
            std::cout << "    Average fit error on original points: " << avg_error << std::endl;
        } else {
            std::cout << "Element " << elemId 
                      << " has only " << combined_points.size() 
                      << " unique points (less than " << MIN_POINTS_FOR_CURVE_FIT 
                      << " required) - cannot perform general conic curve fitting." << std::endl;
        }
    }
    
    // Write results to file
    std::ofstream outFile(output_file);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << " for writing." << std::endl;
        return;
    }
    
    outFile << "Element_ID\tNum_Original_Points\tTotal_Points\ta(x²)\tb(xy)\tc(y)\td(x)\te(const)\tavg_error\n";
    
    for (const auto& [elemId, fit] : elementFits) {
        int numPoints = element_points[elemId].size();
        int totalPoints = elementTotalPoints[elemId];
        
        // Calculate error
        double total_error = 0.0;
        for (const auto& point : element_points[elemId]) {
            double x = point.coordinates[0];
            double y = point.coordinates[1];
            
            // Calculate error as deviation from the conic equation
            double equation_value = y*y + fit.a*x*x + fit.b*x*y + fit.c*y + fit.d*x + fit.e;
            
            // Error is the absolute value of the equation
            double error = std::abs(equation_value);
            total_error += error;
        }
        
        double avg_error = total_error / (numPoints > 0 ? numPoints : 1.0);
        
        outFile << elemId << "\t" 
                << numPoints << "\t"
                << totalPoints << "\t"
                << fit.a << "\t" 
                << fit.b << "\t" 
                << fit.c << "\t"
                << fit.d << "\t" 
                << fit.e << "\t"
                << avg_error << "\n";
    }
    
    outFile.close();
    
    std::cout << "Saved " << elementFits.size() << " element general conic curve fits to " << output_file << std::endl;
    std::cout << "Used all available points from 2-hop neighborhoods." << std::endl;
}

// Helper function to calculate determinant of a 5x5 matrix
double IntersectionPointsUtility::determinant5x5(const Matrix& A) {
    // For a general 5x5 determinant, we'll expand along the first row
    // This is not the most efficient method but it's straightforward to implement
    
    double det = 0.0;
    
    for (int j = 0; j < 5; j++) {
        // Create a 4x4 submatrix by excluding row 0 and column j
        Matrix subMatrix(4, 4);
        
        // Fill the submatrix
        for (int r = 0; r < 4; r++) {
            int row = r + 1;  // Skip row 0
            int subCol = 0;
            
            for (int c = 0; c < 5; c++) {
                if (c != j) {
                    subMatrix(r, subCol) = A(row, c);
                    subCol++;
                }
            }
        }
        
        // Calculate sign: (-1)^(i+j)
        double sign = ((j % 2) == 0) ? 1.0 : -1.0;
        
        // Calculate determinant recursively
        det += sign * A(0, j) * determinant4x4(subMatrix);
    }
    
    return det;
}

// Helper function to calculate determinant of a 4x4 matrix
double IntersectionPointsUtility::determinant4x4(const Matrix& A) {
    // Calculate determinant using cofactor expansion
    
    double det = 0.0;
    
    for (int j = 0; j < 4; j++) {
        // Create a 3x3 submatrix by excluding row 0 and column j
        Matrix subMatrix(3, 3);
        
        // Fill the submatrix
        for (int r = 0; r < 3; r++) {
            int row = r + 1;  // Skip row 0
            int subCol = 0;
            
            for (int c = 0; c < 4; c++) {
                if (c != j) {
                    subMatrix(r, subCol) = A(row, c);
                    subCol++;
                }
            }
        }
        
        // Calculate sign: (-1)^(i+j)
        double sign = ((j % 2) == 0) ? 1.0 : -1.0;
        
        // Calculate determinant recursively
        det += sign * A(0, j) * determinant3x3(subMatrix);
    }
    
    return det;
}

// Helper function to calculate determinant of a 3x3 matrix
double IntersectionPointsUtility::determinant3x3(const Matrix& A) {
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
           A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
           A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
}
/////////////////////////////////////////////////////////
    // Initialize the static container
    std::vector<InterfaceAverageData> InterfaceAveragesUtility::mInterfaceAverageContainer;

    void InterfaceAveragesUtility::ClearInterfaceAverages()
    {
        mInterfaceAverageContainer.clear();
    }
    
    const std::vector<InterfaceAverageData>& InterfaceAveragesUtility::GetInterfaceAverages()
    {
        return mInterfaceAverageContainer;
    }

    void InterfaceAveragesUtility::CollectElementInterfaceAverages(Element::Pointer pElement)
    {
        // Get the geometry and distance values from the element
        auto p_geom = pElement->pGetGeometry();
        
        // Only proceed if the element is properly initialized
        if (!p_geom) return;
        
        // Get the distance values from the element's nodes
        Vector nodal_distances;
        nodal_distances.resize(p_geom->size());
        
        for (unsigned int i = 0; i < p_geom->size(); ++i) {
            nodal_distances[i] = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Check if the element is actually split by the interface
        bool is_split = false;
        const double sign_threshold = 1e-14;
        int pos_count = 0, neg_count = 0;
        for (unsigned int i = 0; i < p_geom->size(); ++i) {
            if (nodal_distances[i] > sign_threshold) {
                pos_count++;
            } else if (nodal_distances[i] < -sign_threshold) {
                neg_count++;
            }
        }
        
        // Element is split only if it has both positive and negative distance values
        is_split = (pos_count > 0 && neg_count > 0);
        
        // Only proceed with calculations if the element is actually split
        if (!is_split) {
            return; // Skip this element as it's not split by the interface
        }
        
        // Structure nodes info (if needed)
        Vector structure_node_id = ZeroVector(p_geom->size());
        
        // Create the modified shape functions utility
        ModifiedShapeFunctions::Pointer p_modified_sh_func;
        
        // Create the appropriate modified shape functions based on geometry type
        if (p_geom->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, nodal_distances, structure_node_id);
        } 
        else if (p_geom->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
            p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, nodal_distances, structure_node_id);
        }
        
        if (p_modified_sh_func) {
            // Prepare variables to store interface data
            Matrix interface_shape_function_neg;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType interface_shape_derivatives_neg;
            Vector interface_weights_neg;
            std::vector<array_1d<double,3>> interface_normals_neg;
            
            try {
                // Calculate interface shape functions, weights, and normals
                p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                    interface_shape_function_neg,
                    interface_shape_derivatives_neg,
                    interface_weights_neg,
                    GeometryData::IntegrationMethod::GI_GAUSS_2);
                    
                p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
                    interface_normals_neg,
                    GeometryData::IntegrationMethod::GI_GAUSS_2);
            }
            catch (std::exception& e) {
                std::cerr << "Error in element " << pElement->Id() << " computing interface: " << e.what() << std::endl;
                return;
            }
                
            // Normalize interface normals
            for (unsigned int gp = 0; gp < interface_normals_neg.size(); ++gp) {
                const double normal_norm = norm_2(interface_normals_neg[gp]);
                if (normal_norm > 1e-10) {
                    interface_normals_neg[gp] /= normal_norm;
                }
            }
            
            // Only proceed if we have interface points
            if (interface_weights_neg.size() > 0) {
                // Create and compute the average interface data
                InterfaceAverageData avg_data;
                avg_data.elementId = pElement->Id();
                avg_data.numberOfPoints = interface_weights_neg.size();
                avg_data.interfaceArea = 0.0;
                
                // Initialize average coordinates and normal to zero
                avg_data.averageCoordinates = ZeroVector(3);
                avg_data.averageNormal = ZeroVector(3);
                
                // Calculate the global coordinates of each interface Gauss point
                std::vector<array_1d<double,3>> interface_global_coords;
                interface_global_coords.resize(interface_weights_neg.size());
                
                // Initialize to zero
                for (unsigned int gp = 0; gp < interface_weights_neg.size(); ++gp) {
                    interface_global_coords[gp] = ZeroVector(3);
                }
                
                // Calculate global coordinates using shape functions
                for (unsigned int gp = 0; gp < interface_weights_neg.size(); ++gp) {
                    for (unsigned int i_node = 0; i_node < p_geom->size(); ++i_node) {
                        const array_1d<double, 3>& r_node_coords = (*p_geom)[i_node].Coordinates();
                        for (unsigned int d = 0; d < 3; ++d) {
                            interface_global_coords[gp][d] += interface_shape_function_neg(gp, i_node) * r_node_coords[d];
                        }
                    }
                }
                
                // Compute total interface area
                for (unsigned int gp = 0; gp < interface_weights_neg.size(); ++gp) {
                    avg_data.interfaceArea += interface_weights_neg[gp];
                }
                
                // Compute weighted average coordinates and normal in a simpler way
                for (unsigned int gp = 0; gp < interface_weights_neg.size(); ++gp) {
                    avg_data.averageCoordinates += interface_weights_neg[gp] * interface_global_coords[gp];
                    avg_data.averageNormal += interface_weights_neg[gp] * interface_normals_neg[gp];
                }
                
                // Finalize the averages by dividing by the total interface area
                if (avg_data.interfaceArea > 1e-10) {
                    avg_data.averageCoordinates /= avg_data.interfaceArea;
                    
                    // Normalize the average normal
                    const double avg_normal_norm = norm_2(avg_data.averageNormal);
                    if (avg_normal_norm > 1e-10) {
                        avg_data.averageNormal /= avg_normal_norm;
                    }
                }
                
                // Add the data to the container
                mInterfaceAverageContainer.push_back(avg_data);
                
                std::cout << "Processed interface averages for element " << pElement->Id() 
                          << " with " << interface_weights_neg.size() << " interface points" << std::endl;
            }
        }
    }
    
    void InterfaceAveragesUtility::ComputeModelPartInterfaceAverages(const ModelPart& rModelPart)
    {
        // Clear previous data
        ClearInterfaceAverages();
        
        // Process each element in the model part
        for (ModelPart::ElementsContainerType::const_iterator it = rModelPart.ElementsBegin(); 
             it != rModelPart.ElementsEnd(); ++it) {
            Element::Pointer p_element = *(it.base());
            CollectElementInterfaceAverages(p_element);
        }
        
        std::cout << "Computed interface averages for " << mInterfaceAverageContainer.size() 
                  << " elements in model part " << rModelPart.Name() << std::endl;
    }
    
    void InterfaceAveragesUtility::SaveInterfaceAveragesToFile(const std::string& filename)
    {
        std::cout << "Saving " << mInterfaceAverageContainer.size() << " interface averages to file: " << filename << std::endl;
        
        std::ofstream outFile(filename);
        
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        outFile << std::fixed << std::setprecision(15);  
        outFile << "Element_ID\tX_avg\tY_avg\tZ_avg\tNormal_X\tNormal_Y\tNormal_Z\tInterface_Area\tNum_Points" << std::endl;
        
        for (const auto& avg : mInterfaceAverageContainer) {
            outFile << avg.elementId << "\t" 
                   << avg.averageCoordinates[0] << "\t" 
                   << avg.averageCoordinates[1] << "\t" 
                   << avg.averageCoordinates[2] << "\t"
                   << avg.averageNormal[0] << "\t"
                   << avg.averageNormal[1] << "\t"
                   << avg.averageNormal[2] << "\t"
                   << avg.interfaceArea << "\t"
                   << avg.numberOfPoints << std::endl;
        }
        
        outFile.close();
        std::cout << "Successfully wrote " << mInterfaceAverageContainer.size() << " interface averages to " << filename << std::endl;
    }
    
    void InterfaceAveragesUtility::DiagnosticOutput(const ModelPart& rModelPart)
    {
        int total_elements = rModelPart.NumberOfElements();
        int split_elements = 0;
        int elements_with_interface = mInterfaceAverageContainer.size();
        double total_interface_area = 0.0;
        
        for (const auto& avg : mInterfaceAverageContainer) {
            total_interface_area += avg.interfaceArea;
        }
        
        for (ModelPart::ElementsContainerType::const_iterator it = rModelPart.ElementsBegin(); 
             it != rModelPart.ElementsEnd(); ++it) {
            auto p_geom = it->pGetGeometry();
            if (!p_geom) continue;
            
            // Get the distance values
            Vector nodal_distances;
            nodal_distances.resize(p_geom->size());
            for (unsigned int i = 0; i < p_geom->size(); ++i) {
                nodal_distances[i] = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE);
            }
            
            // Check if element is split
            const double sign_threshold = 1e-14;
            int pos_count = 0, neg_count = 0;
            for (unsigned int i = 0; i < p_geom->size(); ++i) {
                if (nodal_distances[i] > sign_threshold) {
                    pos_count++;
                } else if (nodal_distances[i] < -sign_threshold) {
                    neg_count++;
                }
            }
            
            bool is_split = (pos_count > 0 && neg_count > 0);
            
            if (is_split) {
                split_elements++;
            }
        }
        
        std::cout << "Interface Averages Diagnostic Output:" << std::endl;
        std::cout << "Total Elements: " << total_elements << std::endl;
        std::cout << "Split Elements: " << split_elements << std::endl;
        std::cout << "Elements with Interface Data: " << elements_with_interface << std::endl;
        std::cout << "Total Interface Area: " << total_interface_area << std::endl;
        
        if (elements_with_interface > 0) {
            std::cout << "Average Interface Area per Element: " << total_interface_area / elements_with_interface << std::endl;
        }
    }
    
    // void InterfaceAveragesUtility::ApplyInterfaceAveragesToModelPart(ModelPart& rModelPart, const std::string& variable_name)
    // {
    //     // This is an example of how you might use the interface averages
    //     // to set a variable on the elements
        
    //     // For example, set the average interface coordinates or normals as element variables
    //     for (const auto& avg_data : mInterfaceAverageContainer) {
    //         // Find the element
    //         auto elem_iterator = rModelPart.Elements().find(avg_data.elementId);
    //         if (elem_iterator != rModelPart.Elements().end()) {
    //             Element& r_element = *elem_iterator;
                
    //             // Example: Store the average normal as an element variable
    //             if (variable_name == "INTERFACE_NORMAL") {
    //                 r_element.SetValue(INTERFACE_NORMAL, avg_data.averageNormal);
    //             }
    //             // Example: Store the average coordinate as an element variable
    //             else if (variable_name == "INTERFACE_POINT") {
    //                 r_element.SetValue(INTERFACE_POINT, avg_data.averageCoordinates);
    //             }
    //             // Example: Store the interface area as an element variable
    //             else if (variable_name == "INTERFACE_AREA") {
    //                 r_element.SetValue(INTERFACE_AREA, avg_data.interfaceArea);
    //             }
    //         }
    //     }
    // }

/**
 * @brief Find an element in the interface averages container by ID
 * @param rInterfaceAverages Container of interface average data
 * @param ElementId ID of the element to find
 * @return Pointer to the interface average data, or nullptr if not found
 */
inline const InterfaceAverageData* FindElementInInterfaceAverages(
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    int ElementId)
{
    for (const auto& avg : rInterfaceAverages) {
        if (avg.elementId == ElementId) {
            return &avg;
        }
    }
    return nullptr;
}

/**
 * @brief Get neighbors of an element through node connectivity
 * @param rModelPart The model part containing elements
 * @param ElementId ID of the element to get neighbors for
 * @return Vector of neighbor element IDs
 */
inline std::vector<int> GetElementNeighbors(const ModelPart& rModelPart, int ElementId)
{
    std::vector<int> neighbors;
    
    // Try to find the element
    auto elem_it = rModelPart.Elements().find(ElementId);
    if (elem_it == rModelPart.Elements().end()) {
        return neighbors;
    }
    
    Element& r_element = *elem_it;
    
    // Get element neighbors through nodes
    for (auto& node : r_element.GetGeometry()) {
        GlobalPointersVector<Element>& r_neighbor_elements = node.GetValue(NEIGHBOUR_ELEMENTS);
        for (auto& neighbor_elem : r_neighbor_elements) {
            if (neighbor_elem.Id() != ElementId && 
                std::find(neighbors.begin(), neighbors.end(), neighbor_elem.Id()) == neighbors.end()) {
                neighbors.push_back(neighbor_elem.Id());
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Get three points (target element and two neighbors) for normal fitting
 * @param rModelPart The model part containing elements
 * @param rInterfaceAverages Container of interface average data
 * @param ElementId ID of the target element
 * @param rPoints Output vector to hold the three points (target + 2 neighbors)
 * @return True if three points were found, false otherwise
 */
inline bool GetThreePointsForFitting(
    const ModelPart& rModelPart,
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    int ElementId,
    std::vector<const InterfaceAverageData*>& rPoints)
{
    // Clear the vector
    rPoints.clear();
    
    // Find the target element data
    const InterfaceAverageData* target_data = FindElementInInterfaceAverages(rInterfaceAverages, ElementId);
    if (!target_data) {
        std::cout << "Element " << ElementId << " not found in interface averages" << std::endl;
        return false;
    }
    
    // Add the target element to the points
    rPoints.push_back(target_data);
    
    // Get immediate neighbors with interface data
    std::vector<int> neighbors = GetElementNeighbors(rModelPart, ElementId);
    std::vector<const InterfaceAverageData*> interface_neighbors;
    
    for (int neighbor_id : neighbors) {
        const InterfaceAverageData* neighbor_data = FindElementInInterfaceAverages(rInterfaceAverages, neighbor_id);
        if (neighbor_data) {
            interface_neighbors.push_back(neighbor_data);
            if (interface_neighbors.size() >= 2) {
                // We found two neighbors with interface data, add them to points
                rPoints.push_back(interface_neighbors[0]);
                rPoints.push_back(interface_neighbors[1]);
                return true;
            }
        }
    }
    
    // If we don't have enough immediate neighbors, try neighbors of neighbors
    if (interface_neighbors.size() < 2) {
        std::set<int> processed_neighbors;
        for (int id : neighbors) {
            processed_neighbors.insert(id);
        }
        processed_neighbors.insert(ElementId);
        
        for (int neighbor_id : neighbors) {
            std::vector<int> second_level_neighbors = GetElementNeighbors(rModelPart, neighbor_id);
            
            for (int nn_id : second_level_neighbors) {
                // Skip already processed elements
                if (processed_neighbors.find(nn_id) != processed_neighbors.end()) {
                    continue;
                }
                
                const InterfaceAverageData* nn_data = FindElementInInterfaceAverages(rInterfaceAverages, nn_id);
                if (nn_data) {
                    interface_neighbors.push_back(nn_data);
                    processed_neighbors.insert(nn_id);
                    
                    if (interface_neighbors.size() >= 2) {
                        // We found two neighbors with interface data (including second level), add them to points
                        rPoints.push_back(interface_neighbors[0]);
                        rPoints.push_back(interface_neighbors[1]);
                        return true;
                    }
                }
            }
        }
    }
    
    // If we still don't have enough points, add what we have
   if (interface_neighbors.size() == 1) {
        rPoints.push_back(interface_neighbors[0]);
        std::cout << "Warning: Only two points available for element " << ElementId << std::endl;
        return false;
    }
    
    std::cout << "Warning: Only one point available for element " << ElementId << std::endl;
    return false;
}

// /**
//  * @brief Fit a normal vector using exactly three points (target element and two neighbors)
//  * @param rModelPart The model part containing elements
//  * @param rInterfaceAverages Container of interface average data
//  * @param ElementId ID of the target element
//  * @return The fitted normal vector, or the original normal if fitting fails
//  */
// inline array_1d<double, 3> FitLinearNormal(
//     const ModelPart& rModelPart,
//     const std::vector<InterfaceAverageData>& rInterfaceAverages,
//     int ElementId,
//     double& a0, double& a1, double& a2,   // Coefficients for nx = a0 + a1*x + a2*y
//     double& b0, double& b1, double& b2)   // Coefficients for ny = b0 + b1*x + b2*y
// {
//     // Get three points for fitting
//     std::vector<const InterfaceAverageData*> points;
//     bool success = GetThreePointsForFitting(rModelPart, rInterfaceAverages, ElementId, points);
    
//     // If we don't have enough points or fitting failed, return the original normal
//     if (!success || points.size() < 3) {
//         const InterfaceAverageData* target_data = points[0]; // This should be the target element
//         return target_data->averageNormal;
//     }
    
//     // Extract coordinates and normal components from the three points
//     const double x1 = points[0]->averageCoordinates[0];
//     const double y1 = points[0]->averageCoordinates[1];
//     const double nx1 = points[0]->averageNormal[0];
//     const double ny1 = points[0]->averageNormal[1];
    
//     const double x2 = points[1]->averageCoordinates[0];
//     const double y2 = points[1]->averageCoordinates[1];
//     const double nx2 = points[1]->averageNormal[0];
//     const double ny2 = points[1]->averageNormal[1];
    
//     const double x3 = points[2]->averageCoordinates[0];
//     const double y3 = points[2]->averageCoordinates[1];
//     const double nx3 = points[2]->averageNormal[0];
//     const double ny3 = points[2]->averageNormal[1];
    
//     // Determine the linear equations for nx and ny components:
//     // nx = a0 + a1*x + a2*y
//     // ny = b0 + b1*x + b2*y
    
//     // Set up matrices for the nx component
//     Matrix A_nx(3, 3);
//     Vector b_nx(3);
    
//     A_nx(0, 0) = 1.0; A_nx(0, 1) = x1; A_nx(0, 2) = y1;
//     A_nx(1, 0) = 1.0; A_nx(1, 1) = x2; A_nx(1, 2) = y2;
//     A_nx(2, 0) = 1.0; A_nx(2, 1) = x3; A_nx(2, 2) = y3;
    
//     b_nx[0] = nx1;
//     b_nx[1] = nx2;
//     b_nx[2] = nx3;
    
//     // Set up matrices for the ny component
//     Matrix A_ny(3, 3);
//     Vector b_ny(3);
    
//     A_ny(0, 0) = 1.0; A_ny(0, 1) = x1; A_ny(0, 2) = y1;
//     A_ny(1, 0) = 1.0; A_ny(1, 1) = x2; A_ny(1, 2) = y2;
//     A_ny(2, 0) = 1.0; A_ny(2, 1) = x3; A_ny(2, 2) = y3;
    
//     b_ny[0] = ny1;
//     b_ny[1] = ny2;
//     b_ny[2] = ny3;
    
//     // Solve the linear systems using Cramer's rule
//     double det = A_nx(0, 0) * (A_nx(1, 1) * A_nx(2, 2) - A_nx(1, 2) * A_nx(2, 1)) -
//                  A_nx(0, 1) * (A_nx(1, 0) * A_nx(2, 2) - A_nx(1, 2) * A_nx(2, 0)) +
//                  A_nx(0, 2) * (A_nx(1, 0) * A_nx(2, 1) - A_nx(1, 1) * A_nx(2, 0));
    
//     // Check if determinant is near zero (singular matrix)
//     if (std::abs(det) < 1.0e-14) {
//         std::cout << "Error: Singular matrix when fitting normal for element " << ElementId << std::endl;
//         return points[0]->averageNormal;  // Return original normal
//     }
    
//     // Create matrices with replaced columns for Cramer's rule
//     Matrix A_nx1 = A_nx;
//     Matrix A_nx2 = A_nx;
//     Matrix A_nx3 = A_nx;
    
//     for (int i = 0; i < 3; i++) {
//         A_nx1(i, 0) = b_nx[i];
//         A_nx2(i, 1) = b_nx[i];
//         A_nx3(i, 2) = b_nx[i];
//     }
    
//     double det_nx1 = A_nx1(0, 0) * (A_nx1(1, 1) * A_nx1(2, 2) - A_nx1(1, 2) * A_nx1(2, 1)) -
//                       A_nx1(0, 1) * (A_nx1(1, 0) * A_nx1(2, 2) - A_nx1(1, 2) * A_nx1(2, 0)) +
//                       A_nx1(0, 2) * (A_nx1(1, 0) * A_nx1(2, 1) - A_nx1(1, 1) * A_nx1(2, 0));
                      
//     double det_nx2 = A_nx2(0, 0) * (A_nx2(1, 1) * A_nx2(2, 2) - A_nx2(1, 2) * A_nx2(2, 1)) -
//                       A_nx2(0, 1) * (A_nx2(1, 0) * A_nx2(2, 2) - A_nx2(1, 2) * A_nx2(2, 0)) +
//                       A_nx2(0, 2) * (A_nx2(1, 0) * A_nx2(2, 1) - A_nx2(1, 1) * A_nx2(2, 0));
                      
//     double det_nx3 = A_nx3(0, 0) * (A_nx3(1, 1) * A_nx3(2, 2) - A_nx3(1, 2) * A_nx3(2, 1)) -
//                       A_nx3(0, 1) * (A_nx3(1, 0) * A_nx3(2, 2) - A_nx3(1, 2) * A_nx3(2, 0)) +
//                       A_nx3(0, 2) * (A_nx3(1, 0) * A_nx3(2, 1) - A_nx3(1, 1) * A_nx3(2, 0));
    
//     // Compute coefficients for nx: nx = a0 + a1*x + a2*y
// a0 = det_nx1 / det;
// a1 = det_nx2 / det;
// a2 = det_nx3 / det;
    
//     // Repeat for ny coefficients
//     Matrix A_ny1 = A_ny;
//     Matrix A_ny2 = A_ny;
//     Matrix A_ny3 = A_ny;
    
//     for (int i = 0; i < 3; i++) {
//         A_ny1(i, 0) = b_ny[i];
//         A_ny2(i, 1) = b_ny[i];
//         A_ny3(i, 2) = b_ny[i];
//     }
    
//     double det_ny1 = A_ny1(0, 0) * (A_ny1(1, 1) * A_ny1(2, 2) - A_ny1(1, 2) * A_ny1(2, 1)) -
//                       A_ny1(0, 1) * (A_ny1(1, 0) * A_ny1(2, 2) - A_ny1(1, 2) * A_ny1(2, 0)) +
//                       A_ny1(0, 2) * (A_ny1(1, 0) * A_ny1(2, 1) - A_ny1(1, 1) * A_ny1(2, 0));
                      
//     double det_ny2 = A_ny2(0, 0) * (A_ny2(1, 1) * A_ny2(2, 2) - A_ny2(1, 2) * A_ny2(2, 1)) -
//                       A_ny2(0, 1) * (A_ny2(1, 0) * A_ny2(2, 2) - A_ny2(1, 2) * A_ny2(2, 0)) +
//                       A_ny2(0, 2) * (A_ny2(1, 0) * A_ny2(2, 1) - A_ny2(1, 1) * A_ny2(2, 0));
                      
//     double det_ny3 = A_ny3(0, 0) * (A_ny3(1, 1) * A_ny3(2, 2) - A_ny3(1, 2) * A_ny3(2, 1)) -
//                       A_ny3(0, 1) * (A_ny3(1, 0) * A_ny3(2, 2) - A_ny3(1, 2) * A_ny3(2, 0)) +
//                       A_ny3(0, 2) * (A_ny3(1, 0) * A_ny3(2, 1) - A_ny3(1, 1) * A_ny3(2, 0));
    
//     // Compute coefficients for ny: ny = b0 + b1*x + b2*y
// b0 = det_ny1 / det;
// b1 = det_ny2 / det;
// b2 = det_ny3 / det;
    
//     // Compute the fitted normal at the target element position
//     double x0 = points[0]->averageCoordinates[0];
//     double y0 = points[0]->averageCoordinates[1];
    
//     double nx = a0 + a1 * x0 + a2 * y0;
//     double ny = b0 + b1 * x0 + b2 * y0;

//     // Print the coefficients
//     std::cout << "Normal fitting coefficients for element " << ElementId << ":" << std::endl;
//     std::cout << "  nx = " << a0 << " + " << a1 << "*x + " << a2 << "*y" << std::endl;
//     std::cout << "  ny = " << b0 << " + " << b1 << "*x + " << b2 << "*y" << std::endl;
    
//     // Compute nz to ensure unit normal
//     // Choose sign to match original normal
//     double nz_sign = (points[0]->averageNormal[2] >= 0) ? 1.0 : -1.0;
//     double nz = nz_sign * std::sqrt(std::max(0.0, 1.0 - nx*nx - ny*ny));
    
//     // Create the final normal vector
//     array_1d<double, 3> normal = ZeroVector(3);
//     normal[0] = nx;
//     normal[1] = ny;
//     normal[2] = nz;
    
//     // Normalize to ensure unit length
//     double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
//     if (norm > 1.0e-10) {
//         normal[0] /= norm;
//         normal[1] /= norm;
//         normal[2] /= norm;
//     }
    
//     return normal;
// }

/**
 * @brief Fit a normal vector using exactly three points with focus on stability
 * @param rModelPart The model part containing elements
 * @param rInterfaceAverages Container of interface average data
 * @param ElementId ID of the target element
 * @return The fitted normal vector, or the original normal if fitting fails
 */
inline array_1d<double, 3> FitLinearNormal(
    const ModelPart& rModelPart,
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    int ElementId,
    double& a0, double& a1, double& a2,   // Coefficients for nx = a0 + a1*x + a2*y
    double& b0, double& b1, double& b2)   // Coefficients for ny = b0 + b1*x + b2*y
{
    // Find target element data
    const InterfaceAverageData* target_data = nullptr;
    for (const auto& avg : rInterfaceAverages) {
        if (avg.elementId == ElementId) {
            target_data = &avg;
            break;
        }
    }
    
    // If target element not found, return a default normal
    if (!target_data) {
        std::cout << "Error: Element " << ElementId << " not found in interface averages" << std::endl;
        array_1d<double, 3> default_normal = ZeroVector(3);
        default_normal[0] = 1.0;  // Default to x-direction
        
        // Set coefficients to zero
        a0 = 1.0; a1 = 0.0; a2 = 0.0;
        b0 = 0.0; b1 = 0.0; b2 = 0.0;
        
        return default_normal;
    }
    
    // Store the target element coordinates and normal
    const double target_x = target_data->averageCoordinates[0];
    const double target_y = target_data->averageCoordinates[1];
    const double target_nx = target_data->averageNormal[0];
    const double target_ny = target_data->averageNormal[1];
    const double target_nz = target_data->averageNormal[2];
    
    // Get neighbors without using the GetThreePointsForFitting function
    std::vector<const InterfaceAverageData*> neighbors;
    
    // Get direct neighbors from ModelPart
    std::vector<int> neighbor_ids = GetElementNeighbors(rModelPart, ElementId);
    
    // Find matching neighbors in interface averages
    for (int id : neighbor_ids) {
        for (const auto& avg : rInterfaceAverages) {
            if (avg.elementId == id) {
                neighbors.push_back(&avg);
                break;
            }
        }
        
        // Break if we found enough neighbors
        if (neighbors.size() >= 2) {
            break;
        }
    }
    
    // If we don't have enough neighbors, search for neighbors of neighbors
    if (neighbors.size() < 2) {
        for (int id : neighbor_ids) {
            std::vector<int> second_level = GetElementNeighbors(rModelPart, id);
            
            for (int id2 : second_level) {
                // Skip if it's the target or already processed
                if (id2 == ElementId || std::find(neighbor_ids.begin(), neighbor_ids.end(), id2) != neighbor_ids.end()) {
                    continue;
                }
                
                for (const auto& avg : rInterfaceAverages) {
                    if (avg.elementId == id2) {
                        neighbors.push_back(&avg);
                        break;
                    }
                }
                
                // Break if we found enough neighbors
                if (neighbors.size() >= 2) {
                    break;
                }
            }
            
            if (neighbors.size() >= 2) {
                break;
            }
        }
    }
    
    // If we still don't have enough neighbors, just return the original normal
    if (neighbors.size() < 2) {
        std::cout << "Warning: Not enough neighbors for element " << ElementId << ". Using original normal." << std::endl;
        
        array_1d<double, 3> original_normal = ZeroVector(3);
        original_normal[0] = target_nx;
        original_normal[1] = target_ny;
        original_normal[2] = target_nz;
        
        // Set coefficients (constants only, no derivatives)
        a0 = target_nx; a1 = 0.0; a2 = 0.0;
        b0 = target_ny; b1 = 0.0; b2 = 0.0;
        
        return original_normal;
    }
    
    // We now have the target and two neighbors
    const double x1 = target_x;
    const double y1 = target_y;
    const double nx1 = target_nx;
    const double ny1 = target_ny;
    
    const double x2 = neighbors[0]->averageCoordinates[0];
    const double y2 = neighbors[0]->averageCoordinates[1];
    const double nx2 = neighbors[0]->averageNormal[0];
    const double ny2 = neighbors[0]->averageNormal[1];
    
    const double x3 = neighbors[1]->averageCoordinates[0];
    const double y3 = neighbors[1]->averageCoordinates[1];
    const double nx3 = neighbors[1]->averageNormal[0];
    const double ny3 = neighbors[1]->averageNormal[1];
    
    // Calculate coefficients directly using determinants
    // This avoids potential issues with Gaussian elimination
    
    // Compute the determinant of the coefficient matrix
    double det = (x2*y3 - x3*y2) - x1*(y3 - y2) + y1*(x3 - x2);
    
    if (std::abs(det) < 1.0e-10) {
        std::cout << "Warning: Singular matrix for element " << ElementId << ". Using original normal." << std::endl;
        
        array_1d<double, 3> original_normal = ZeroVector(3);
        original_normal[0] = target_nx;
        original_normal[1] = target_ny;
        original_normal[2] = target_nz;
        
        // Set coefficients (constants only, no derivatives)
        a0 = target_nx; a1 = 0.0; a2 = 0.0;
        b0 = target_ny; b1 = 0.0; b2 = 0.0;
        
        return original_normal;
    }
    
    // Compute coefficients for nx
    a0 = ((nx1*(x2*y3 - x3*y2)) + (nx2*(x3*y1 - x1*y3)) + (nx3*(x1*y2 - x2*y1))) / det;
    a1 = ((nx1*(y2 - y3)) + (nx2*(y3 - y1)) + (nx3*(y1 - y2))) / det;
    a2 = ((nx1*(x3 - x2)) + (nx2*(x1 - x3)) + (nx3*(x2 - x1))) / det;
    
    // Compute coefficients for ny
    b0 = ((ny1*(x2*y3 - x3*y2)) + (ny2*(x3*y1 - x1*y3)) + (ny3*(x1*y2 - x2*y1))) / det;
    b1 = ((ny1*(y2 - y3)) + (ny2*(y3 - y1)) + (ny3*(y1 - y2))) / det;
    b2 = ((ny1*(x3 - x2)) + (ny2*(x1 - x3)) + (ny3*(x2 - x1))) / det;
    
    // Evaluate the fitted normal at the target point (which should match the original normal)
    double fitted_nx = a0 + a1*x1 + a2*y1;
    double fitted_ny = b0 + b1*x1 + b2*y1;
    
    // Create the fitted normal vector
    array_1d<double, 3> fitted_normal = ZeroVector(3);
    fitted_normal[0] = fitted_nx;
    fitted_normal[1] = fitted_ny;
    fitted_normal[2] = target_nz; // Preserve the original z-component
    
    // Normalize the vector
    double norm = std::sqrt(fitted_normal[0]*fitted_normal[0] + 
                           fitted_normal[1]*fitted_normal[1] + 
                           fitted_normal[2]*fitted_normal[2]);
    
    if (norm > 1.0e-10) {
        fitted_normal[0] /= norm;
        fitted_normal[1] /= norm;
        fitted_normal[2] /= norm;
    }
    
    // Calculate curvature
    double curvature = a1 + b2;
    
    // Output information
    std::cout << "Normal fitting for element " << ElementId << ":" << std::endl;
    std::cout << "  nx = " << a0 << " + " << a1 << "*x + " << a2 << "*y" << std::endl;
    std::cout << "  ny = " << b0 << " + " << b1 << "*x + " << b2 << "*y" << std::endl;
    std::cout << "  nz = " << target_nz << " (preserved from original)" << std::endl;
    std::cout << "  Curvature (a1 + b2): " << curvature << std::endl;
    
    // Verify fit at target element
    double nx_error = std::abs(fitted_nx - target_nx);
    double ny_error = std::abs(fitted_ny - target_ny);
    
    if (nx_error > 1.0e-10 || ny_error > 1.0e-10) {
        std::cout << "Warning: Target element " << ElementId << " normal not exactly fitted:" << std::endl;
        std::cout << "  Original: [" << target_nx << ", " << target_ny << "]" << std::endl;
        std::cout << "  Fitted:   [" << fitted_nx << ", " << fitted_ny << "]" << std::endl;
        std::cout << "  Error:    [" << nx_error << ", " << ny_error << "]" << std::endl;
    }
    
    return fitted_normal;
}

/**
 * @brief Save fitted normals to a file with robust error handling
 * @param rModelPart The model part containing elements
 * @param rInterfaceAverages Container of interface average data
 * @param Filename Filename to save the results
 */
inline void SaveFittedNormalsToFile(
    const ModelPart& rModelPart,
    const std::vector<InterfaceAverageData>& rInterfaceAverages,
    const std::string& Filename)
{
    try {
        // Open output file
        std::ofstream outFile(Filename);
        
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << Filename << " for writing." << std::endl;
            return;
        }
        
        // Write header
        outFile << std::fixed << std::setprecision(15);
        outFile << "Element_ID\tX\tY\tZ\tOrig_NX\tOrig_NY\tOrig_NZ\tFit_NX\tFit_NY\tFit_NZ\tDiff\t"
                << "a0\ta1\ta2\tb0\tb1\tb2\tCurvature" << std::endl;
        
        // Process each element one by one
        int successful_elements = 0;
        
        for (const auto& avg : rInterfaceAverages) {
            int element_id = avg.elementId;
            double a0 = 0.0, a1 = 0.0, a2 = 0.0;
            double b0 = 0.0, b1 = 0.0, b2 = 0.0;
            array_1d<double, 3> fitted_normal;
            double curvature = 0.0;
            double diff_magnitude = 0.0;
            
            try {
                // Calculate a linearly fitted normal
                fitted_normal = FitLinearNormal(rModelPart, rInterfaceAverages, element_id, 
                                             a0, a1, a2, b0, b1, b2);
                
                // Calculate difference between fitted and original normal
                array_1d<double, 3> normal_diff = fitted_normal - avg.averageNormal;
                diff_magnitude = std::sqrt(normal_diff[0]*normal_diff[0] + 
                                        normal_diff[1]*normal_diff[1] + 
                                        normal_diff[2]*normal_diff[2]);
                
                // Calculate curvature
                curvature = a1 + b2;
                successful_elements++;
            }
            catch (const std::exception& e) {
                std::cerr << "Error processing element " << element_id << ": " << e.what() << std::endl;
                
                // Use original normal as fallback
                fitted_normal = avg.averageNormal;
                a0 = avg.averageNormal[0];
                b0 = avg.averageNormal[1];
                curvature = 0.0;
                diff_magnitude = 0.0;
            }
            
            // Write to file
            outFile << element_id << "\t"
                    << avg.averageCoordinates[0] << "\t"
                    << avg.averageCoordinates[1] << "\t"
                    << avg.averageCoordinates[2] << "\t"
                    << avg.averageNormal[0] << "\t"
                    << avg.averageNormal[1] << "\t"
                    << avg.averageNormal[2] << "\t"
                    << fitted_normal[0] << "\t"
                    << fitted_normal[1] << "\t"
                    << fitted_normal[2] << "\t"
                    << diff_magnitude << "\t"
                    << a0 << "\t" << a1 << "\t" << a2 << "\t"
                    << b0 << "\t" << b1 << "\t" << b2 << "\t"
                    << curvature << std::endl;
        }
        
        outFile.close();
        std::cout << "Saved " << rInterfaceAverages.size() << " fitted normals with curvature to " << Filename 
                  << " (" << successful_elements << " successfully fitted)" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error in SaveFittedNormalsToFile: " << e.what() << std::endl;
    }
}

/**
 * @brief Apply fitted normals to all interface elements in the model part
 * @param rModelPart The model part containing elements
 * @param rInterfaceAverages Container of interface average data
 * @param StoreOriginalNormal Whether to store the original normal for comparison
 */
// inline void ApplyFittedNormalsToModelPart(
//     ModelPart& rModelPart,
//     const std::vector<InterfaceAverageData>& rInterfaceAverages,
//     bool StoreOriginalNormal = true)
// {
//     for (const auto& avg : rInterfaceAverages) {
//         int element_id = avg.elementId;
        
//         // Calculate a linearly fitted normal
//         array_1d<double, 3> fitted_normal = FitLinearNormal(rModelPart, rInterfaceAverages, element_id);
        
        // // Find the element
        // auto elem_it = rModelPart.Elements().find(element_id);
        // if (elem_it != rModelPart.Elements().end()) {
        //     Element& r_element = *elem_it;
            
        //     // Store the fitted normal on the element
        //     r_element.SetValue(NORMAL, fitted_normal);
            
        //     // Store the original normal for comparison if needed
        //     if (StoreOriginalNormal) {
        //         r_element.SetValue(ORIGINAL_NORMAL, avg.averageNormal);
        //     }
            
        //     // Calculate and store the difference for analysis
        //     array_1d<double, 3> normal_diff = fitted_normal - avg.averageNormal;
        //     double diff_magnitude = std::sqrt(normal_diff[0]*normal_diff[0] + 
        //                                     normal_diff[1]*normal_diff[1] + 
        //                                     normal_diff[2]*normal_diff[2]);
            
        //     r_element.SetValue(NORMAL_SMOOTHING_DIFF, diff_magnitude);
        // }
    // }
    
//     std::cout << "Applied fitted normals to " << rInterfaceAverages.size() 
//                 << " interface elements" << std::endl;
// }

}
}
