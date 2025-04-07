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

//     void IntersectionPointsUtility::AddIntersectionPoint(int elementId, int pointId, const array_1d<double, 3>& coordinates)
// {
//     IntersectionPointData point;
//     point.elementId = elementId;
//     point.pointId = pointId;
    
//     point.coordinates[0] = coordinates[0];
//     point.coordinates[1] = coordinates[1];
//     point.coordinates[2] = coordinates[2];
    
//     g_IntersectionPointsContainer.push_back(point);
// }

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
            
    //         // Fallback to test points if no interface points were found
    //         for (int i = 0; i < 2; i++) {
    //             IntersectionPointData point;
    //             point.elementId = elementId;
    //             point.pointId = i;
                
    //             // Sample coordinates
    //             point.coordinates[0] = 0.5 + 0.1*i;  // X
    //             point.coordinates[1] = 0.5 - 0.1*i;  // Y
    //             point.coordinates[2] = 0.0;          // Z (0 for 2D)
                
    //             g_IntersectionPointsContainer.push_back(point);
    //         }
            
    //         std::cout << "Used test points for element " << elementId << " (no interface points found)" << std::endl;
               }
        } catch (std::exception& e) {
        std::cerr << "Error extracting interface points: " << e.what() << std::endl;
        
        // // Fallback to test points if there's an error
        // for (int i = 0; i < 2; i++) {
        //     IntersectionPointData point;
        //     point.elementId = elementId;
        //     point.pointId = i;
            
        //     // Sample coordinates
        //     point.coordinates[0] = 0.5 + 0.1*i;  // X
        //     point.coordinates[1] = 0.5 - 0.1*i;  // Y
        //     point.coordinates[2] = 0.0;          // Z (0 for 2D)
            
        //     g_IntersectionPointsContainer.push_back(point);
        // }
        
        // std::cout << "Fell back to test points for element " << elementId << " due to error: " << e.what() << std::endl;
        }
}
////////////////////////////////////////////////
void IntersectionPointsUtility::ProcessIntersectionPointsAndFitCurves(const std::string& output_file)
{
    // Get all intersection points
    const auto& points = g_IntersectionPointsContainer;
    
    if (points.empty()) {
        std::cout << "No intersection points available for curve fitting." << std::endl;
        return;
    }
    
    std::cout << "Starting curve fitting with " << points.size() << " intersection points." << std::endl;
    
    // First, create a map of points by their coordinates
    // This will help us identify which points are shared between elements
    std::map<std::pair<double, double>, std::vector<int>> point_to_elements;
    std::map<int, std::vector<IntersectionPointData>> element_points;
    
    // Group points by element and build point->elements mapping
    for (const auto& point : points) {
        int elemId = point.elementId;
        
        // Round coordinates to handle floating point precision
        double x = std::round(point.coordinates[0] * 10000000.0) / 10000000.0;
        double y = std::round(point.coordinates[1] * 10000000.0) / 10000000.0;
        std::pair<double, double> coord_key(x, y);
        
        // Add this element to the list for this point
        point_to_elements[coord_key].push_back(elemId);
        
        // Add this point to the element's list
        element_points[elemId].push_back(point);
    }
    
    // Find connected element pairs (elements that share intersection points)
    std::map<int, std::set<int>> element_neighbors;
    
    for (const auto& [coord, elements] : point_to_elements) {
        // If this point belongs to multiple elements, they are neighbors
        for (size_t i = 0; i < elements.size(); ++i) {
            for (size_t j = i+1; j < elements.size(); ++j) {
                int elem1 = elements[i];
                int elem2 = elements[j];
                
                // Mark as neighbors
                element_neighbors[elem1].insert(elem2);
                element_neighbors[elem2].insert(elem1);
            }
        }
    }
    
    // Now, for each element, fit a curve using its points and its neighbors' points
    std::map<int, ElementBasedCurveFitter::QuadraticCoefficients> elementFits;
    
    for (const auto& [elemId, neighbors] : element_neighbors) {
        // Collect all points from this element and its neighbors
        std::vector<IntersectionPointData> combined_points = element_points[elemId];
        
        for (int neighborId : neighbors) {
            // Add neighbor's points
            combined_points.insert(combined_points.end(), 
                                 element_points[neighborId].begin(), 
                                 element_points[neighborId].end());
        }
        
        // Remove duplicate points
        std::map<std::pair<double, double>, IntersectionPointData> unique_points;
        for (const auto& point : combined_points) {
            double x = std::round(point.coordinates[0] * 10000000.0) / 10000000.0;
            double y = std::round(point.coordinates[1] * 10000000.0) / 10000000.0;
            std::pair<double, double> key(x, y);
            unique_points[key] = point;
        }
        
        // Convert back to vector
        combined_points.clear();
        for (const auto& [_, point] : unique_points) {
            combined_points.push_back(point);
        }
        
        // Only fit if we have enough points
        if (combined_points.size() >= 3) {
            // Fit curve for this element group
            ElementBasedCurveFitter::QuadraticCoefficients fit;
            
            // Setup matrices for least squares fitting
            double sum_x = 0.0, sum_y = 0.0;
            double sum_x2 = 0.0, sum_x3 = 0.0, sum_x4 = 0.0;
            double sum_xy = 0.0, sum_x2y = 0.0;
            
            for (const auto& point : combined_points) {
                double x = point.coordinates[0];
                double y = point.coordinates[1];
                
                sum_x += x;
                sum_y += y;
                sum_x2 += x * x;
                sum_x3 += x * x * x;
                sum_x4 += x * x * x * x;
                sum_xy += x * y;
                sum_x2y += x * x * y;
            }
            
            int n = combined_points.size();
            
            // Create the system of equations
            Matrix A(3, 3);
            Vector b(3);
            
            A(0, 0) = sum_x4;    A(0, 1) = sum_x3;    A(0, 2) = sum_x2;
            A(1, 0) = sum_x3;    A(1, 1) = sum_x2;    A(1, 2) = sum_x;
            A(2, 0) = sum_x2;    A(2, 1) = sum_x;     A(2, 2) = n;
            
            b[0] = sum_x2y;
            b[1] = sum_xy;
            b[2] = sum_y;
            
            // Solve using Cramer's rule for simplicity
            double det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
                         A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                         A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
                         

            // Solve using Cramer's rule
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
                
            fit.a = det1 / det;
            fit.b = det2 / det;
            fit.c = det3 / det;
            
            elementFits[elemId] = fit;
            
            std::cout << "Element " << elemId 
                      << " with " << combined_points.size() 
                      << " points (including neighbors): y = " 
                      << fit.a << "x² + " << fit.b << "x + " << fit.c << std::endl;
        } else {
            std::cout << "Element " << elemId 
                      << " still has only " << combined_points.size() 
                      << " unique points (less than 3) - skipping curve fitting." << std::endl;
        }
    }
    
    // Save results to file
    std::ofstream outFile(output_file);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << " for writing." << std::endl;
        return;
    }
    
    outFile << "Element_ID\tNum_Points\ta\tb\tc\n";
    
    for (const auto& fit_pair : elementFits) {
        int elemId = fit_pair.first;
        const auto& fit = fit_pair.second;
        int numPoints = element_points[elemId].size();
        
        outFile << elemId << "\t" 
                << numPoints << "\t"
                << fit.a << "\t" 
                << fit.b << "\t" 
                << fit.c << "\n";
    }
    
    outFile.close();
    
    std::cout << "Saved " << elementFits.size() << " element curve fits to " << output_file << std::endl;
}
/////////////////////////////////////////////////////////
}
}
