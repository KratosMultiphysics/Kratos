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

#include <iomanip>  // For std::setprecision
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
    // Define the global containers
    std::vector<IntersectionPointData> g_IntersectionPointsContainer;
    std::vector<InterfaceAverageData> InterfaceAveragesUtility::mInterfaceAverageContainer;
    std::vector<IntersectionDataWithNormal> g_IntersectionDataWithNormalContainer;

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
    const double sign_threshold = 1e-15;
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

/**
 * @brief Calculate and store intersection lengths for 2D elements
 * @param rModelPart The model part containing elements
 * @return Number of elements with valid intersection lengths
 */
int CalculateAndStoreElementIntersectionLengths(ModelPart& rModelPart)
{
    // Count of elements with valid intersection lengths
    int count = 0;
    
    // Clear any previous intersection points
    KratosDropletDynamics::IntersectionPointsUtility::ClearIntersectionPoints();
    
    // Collect intersection points for all elements
    for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it) {
        Element::Pointer pElement = *(it.base());
        KratosDropletDynamics::IntersectionPointsUtility::CollectElementIntersectionPoints(pElement);
    }
    
    // Get all intersection points
    const auto& points = KratosDropletDynamics::IntersectionPointsUtility::GetIntersectionPoints();
    
    // Group points by element ID
    std::map<int, std::vector<array_1d<double, 3>>> elementPoints;
    for (const auto& point : points) {
        elementPoints[point.elementId].push_back(point.coordinates);
    }
    
    // Initialize INTERSECTION_LENGTH_2D to zero for all elements
    for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it) {
        it->SetValue(INTERSECTION_LENGTH_2D, 0.0);
    }
    
    // Calculate length for each element and store it
    for (const auto& [elementId, elementPointList] : elementPoints) {
        if (elementPointList.size() == 2) { // For 2D elements
            // Calculate distance between the two points
            double dx = elementPointList[0][0] - elementPointList[1][0];
            double dy = elementPointList[0][1] - elementPointList[1][1];
            double dz = elementPointList[0][2] - elementPointList[1][2];
            
            double length = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            // Find the element and store the value
            auto elem_iterator = rModelPart.Elements().find(elementId);
            if (elem_iterator != rModelPart.Elements().end()) {
                elem_iterator->SetValue(INTERSECTION_LENGTH_2D, length);
                count++;
            }
        }
    }
    
    return count;
}

/**
 * @brief Get intersection length from an element
 * @param rElement Element to get the intersection length from
 * @return Intersection length value
 */
double GetElementIntersectionLength(const Element& rElement)
{
    return rElement.GetValue(INTERSECTION_LENGTH_2D);
}

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


/**
 * @brief Find an element in the interface averages container by ID
 * @param rInterfaceAverages Container of interface average data
 * @param ElementId ID of the element to find
 * @return Pointer to the interface average data, or nullptr if not found
 */
const InterfaceAverageData* FindElementInInterfaceAverages(
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
// std::vector<int> GetElementNeighbors(const ModelPart& rModelPart, int ElementId)
// {
//     std::vector<int> neighbors;
    
//     // Try to find the element
//     auto elem_it = rModelPart.Elements().find(ElementId);
//     if (elem_it == rModelPart.Elements().end()) {
//         return neighbors;
//     }
    
//     Element& r_element = *elem_it;
    
//     // Get element neighbors through nodes
//     for (auto& node : r_element.GetGeometry()) {
//         GlobalPointersVector<Element>& r_neighbor_elements = node.GetValue(NEIGHBOUR_ELEMENTS);
//         for (auto& neighbor_elem : r_neighbor_elements) {
//             if (neighbor_elem.Id() != ElementId && 
//                 std::find(neighbors.begin(), neighbors.end(), neighbor_elem.Id()) == neighbors.end()) {
//                 neighbors.push_back(neighbor_elem.Id());
//             }
//         }
//     }
    
//     return neighbors;
// }

/**
 * @brief Get neighbors of an element through elemental connectivity
 * @param rModelPart The model part containing elements
 * @param ElementId ID of the element to get neighbors for
 * @return Vector of neighbor element IDs
 */
std::vector<int> GetElementNeighbors(const ModelPart& rModelPart, int ElementId)
{
    std::vector<int> neighbors;
    
    // Find the element
    auto elem_it = rModelPart.Elements().find(ElementId);
    if (elem_it == rModelPart.Elements().end()) {
        return neighbors;
    }
    
    Element& r_element = *elem_it;
    
    // Check if the element has NEIGHBOUR_ELEMENTS variable
    if (r_element.Has(NEIGHBOUR_ELEMENTS)) {
        // Get direct element neighbors
        GlobalPointersVector<Element>& r_neighbor_elements = r_element.GetValue(NEIGHBOUR_ELEMENTS);
        
        // Loop through neighbors with null check
        for (unsigned int i = 0; i < r_neighbor_elements.size(); i++) {
            // Check for null pointers
            if (r_neighbor_elements(i).get() != nullptr) {
                neighbors.push_back(r_neighbor_elements(i)->Id());
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
bool GetThreePointsForFitting(
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

/**
 * @brief Fit a normal vector using exactly three points with focus on stability
 * @param rModelPart The model part containing elements
 * @param rInterfaceAverages Container of interface average data
 * @param ElementId ID of the target element
 * @return The fitted normal vector, or the original normal if fitting fails
 */
array_1d<double, 3> FitLinearNormal(
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
void SaveFittedNormalsToFile(
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
 * @brief Set the ELEMENT_CUT_NORMAL variable for all elements that are cut by the interface
 * @param rModelPart The model part containing elements
 * @return Number of elements where normal was set
 */
int SetElementCutNormals(ModelPart& rModelPart)
{
    int count = 0;
    
    // First compute interface averages if not already done
    if (InterfaceAveragesUtility::GetInterfaceAverages().empty()) {
        InterfaceAveragesUtility::ComputeModelPartInterfaceAverages(rModelPart);
    }
    
    // Get interface averages
    const auto& interface_averages = InterfaceAveragesUtility::GetInterfaceAverages();
    
    // Set normal for each element in the interface averages
    for (const auto& avg : interface_averages) {
        auto elem_it = rModelPart.Elements().find(avg.elementId);
        if (elem_it != rModelPart.Elements().end()) {
            Element& r_element = *elem_it;
            
            // Set the normal vector as element variable
            r_element.SetValue(ELEMENT_CUT_NORMAL, avg.averageNormal);
            count++;
        }
    }
    
    return count;
}

/**
 * @brief Get X component of element cut normal
 * @param rElement Element to get the normal from
 * @return X component of the normal
 */
double GetElementCutNormalX(const Element& rElement)
{
    return rElement.GetValue(ELEMENT_CUT_NORMAL)[0];
}

/**
 * @brief Get Y component of element cut normal
 * @param rElement Element to get the normal from
 * @return Y component of the normal
 */
double GetElementCutNormalY(const Element& rElement)
{
    return rElement.GetValue(ELEMENT_CUT_NORMAL)[1];
}

/**
 * @brief Get Z component of element cut normal
 * @param rElement Element to get the normal from
 * @return Z component of the normal
 */
double GetElementCutNormalZ(const Element& rElement)
{
    return rElement.GetValue(ELEMENT_CUT_NORMAL)[2];
}

/**
 * @brief Clear the intersection data with normal container
 */
void ClearIntersectionDataWithNormal()
{
    g_IntersectionDataWithNormalContainer.clear();
}

/**
 * @brief Get the container of intersection data with normals
 * @return Reference to the container
 */
const std::vector<IntersectionDataWithNormal>& GetIntersectionDataWithNormal()
{
    return g_IntersectionDataWithNormalContainer;
}

/**
 * @brief Populate the intersection data with normal container
 * @param rModelPart The model part containing elements
 * @return Number of elements processed
 */
int CollectIntersectionDataWithNormal(ModelPart& rModelPart)
{
    // Clear previous data
    ClearIntersectionDataWithNormal();
    
    // Make sure we have length data
    int num_lengths = CalculateAndStoreElementIntersectionLengths(rModelPart);
    
    // Make sure we have interface averages
    if (InterfaceAveragesUtility::GetInterfaceAverages().empty()) {
        InterfaceAveragesUtility::ComputeModelPartInterfaceAverages(rModelPart);
    }
    
    // Get interface averages
    const auto& interface_averages = InterfaceAveragesUtility::GetInterfaceAverages();
    
    // Create a map for easy lookup of interface data by element ID
    std::map<int, const InterfaceAverageData*> avg_map;
    for (const auto& avg : interface_averages) {
        avg_map[avg.elementId] = &avg;
    }
    
    // Process each element with intersection length data
    int count = 0;
    for (ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); 
         it != rModelPart.ElementsEnd(); ++it) {
        Element& r_element = *it;
        double length = r_element.GetValue(INTERSECTION_LENGTH_2D);
        
        // Only include elements with non-zero intersection length
        if (length > 0.0) {
            int element_id = r_element.Id();
            
            // Create the data object
            IntersectionDataWithNormal data;
            data.elementId = element_id;
            data.intersectionLength = length;
            
            // Try to find normal data
            auto avg_it = avg_map.find(element_id);
            if (avg_it != avg_map.end()) {
                const InterfaceAverageData* avg = avg_it->second;
                
                // Copy normal and coordinates
                data.normal = avg->averageNormal;
                data.coordinates = avg->averageCoordinates;
            }
            
            // Add to the container
            g_IntersectionDataWithNormalContainer.push_back(data);
            count++;
        }
    }
    
    std::cout << "Collected " << count << " elements with intersection data and normals" << std::endl;
    return count;
}

/**
 * @brief Save the intersection data with normals to file
 * @param filename Filename to save the results
 */
void SaveIntersectionDataWithNormalToFile(const std::string& filename)
{
    std::ofstream outFile(filename);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    
    outFile << std::fixed << std::setprecision(15);
    outFile << "Element_ID\tIntersection_Length\tNormal_X\tNormal_Y\tNormal_Z\tCoord_X\tCoord_Y\tCoord_Z" << std::endl;
    
    for (const auto& data : g_IntersectionDataWithNormalContainer) {
        outFile << data.elementId << "\t"
               << data.intersectionLength << "\t"
               << data.normal[0] << "\t"
               << data.normal[1] << "\t" 
               << data.normal[2] << "\t"
               << data.coordinates[0] << "\t"
               << data.coordinates[1] << "\t"
               << data.coordinates[2] << std::endl;
    }
    
    outFile.close();
    std::cout << "Successfully wrote " << g_IntersectionDataWithNormalContainer.size() 
              << " elements with intersection data and normals to " << filename << std::endl;
}
//////////////////////////////////////////////////////////////////////
/**
 * @brief Compute an averaged normal for an element considering neighboring cut elements
 * @param rModelPart The model part containing elements
 * @param ElementId ID of the target element
 * @param NeighborLevels Number of neighbor levels to consider (0 = no neighbors, 1 = direct neighbors only)
 * @return Weighted average normal vector
 */
array_1d<double, 3> ComputeAveragedElementNormal(
    const ModelPart& rModelPart,
    int ElementId,
    int NeighborLevels)
{
    // Initialize output normal
    array_1d<double, 3> averaged_normal = ZeroVector(3);
    double total_weight = 0.0;
    
    // Get the element
    auto elem_it = rModelPart.Elements().find(ElementId);
    if (elem_it == rModelPart.Elements().end()) {
        return averaged_normal; // Return zero vector if element not found
    }
    
    Element& r_element = *elem_it;
    
    // Add the current element to the calculation
    double elem_length = r_element.GetValue(INTERSECTION_LENGTH_2D);
    array_1d<double, 3> elem_normal = r_element.GetValue(ELEMENT_CUT_NORMAL);
    
    // Skip if element is not cut
    if (elem_length <= 0.0) {
        return elem_normal; // Return element's normal directly if not cut
    }
    
    // Special case: if neighbor hops is 0, return only the element's own normal
    if (NeighborLevels == 0) {
        return elem_normal;
    }
    
    // Add the element's contribution
    averaged_normal += elem_length * elem_normal;
    total_weight += elem_length;
    
    // Find and collect all neighbors up to specified levels
    std::set<int> processed_elements;
    std::set<int> current_level_ids;
    std::set<int> next_level_ids;
    
    // Start with current element
    processed_elements.insert(ElementId);
    current_level_ids.insert(ElementId);
    
    // Traverse neighbor levels
    for (int level = 0; level < NeighborLevels; level++) {
        next_level_ids.clear();

        ///////////////////
        // Calculate weight scaling factor for this level (2/(1+e^(level))
        double level_weight = 2/(1+exp(level+1));  // 0.537883, 0.238406, etc.
        ///////////////////
        // Process each element in current level
        for (int current_id : current_level_ids) {
            // Find direct neighbors
            std::vector<int> neighbors = GetElementNeighbors(rModelPart, current_id);
            //KRATOS_INFO("element Neighbors") <<"level= "<< level+1 <<", "<<"level_weight= "<< level_weight << ", "<< ElementId << ": " <<neighbors << std::endl;
            
            for (int neighbor_id : neighbors) {
                // Skip if already processed
                if (processed_elements.find(neighbor_id) != processed_elements.end()) {
                    continue;
                }
                
                // Get neighbor element
                auto neighbor_it = rModelPart.Elements().find(neighbor_id);
                if (neighbor_it == rModelPart.Elements().end()) {
                    continue;
                }
                
                Element& r_neighbor = *neighbor_it;
                
                // Get neighbor intersection length and normal
                double neighbor_length = r_neighbor.GetValue(INTERSECTION_LENGTH_2D);
                
                // Only include if it's a cut element
                if (neighbor_length > 0.0) {
                    array_1d<double, 3> neighbor_normal = r_neighbor.GetValue(ELEMENT_CUT_NORMAL);
                    
                    // Make sure the normal orientation is consistent
                    // (dot product > 0 means normals point in roughly same direction)
                    if (inner_prod(elem_normal, neighbor_normal) < 0) {
                        neighbor_normal = -neighbor_normal; // Flip the normal if pointing in opposite direction
                    }
                    
                    // // Add contribution
                    // averaged_normal += neighbor_length * neighbor_normal;
                    // total_weight += neighbor_length;
                    //////////////////////////////
                    // Add contribution with level-based weight scaling
                    averaged_normal += level_weight * neighbor_length * neighbor_normal;
                    total_weight += level_weight * neighbor_length;
                    /////////////////////////////
                }
                
                // Add to next level
                next_level_ids.insert(neighbor_id);
                processed_elements.insert(neighbor_id);
            }
        }
        
        // Update for next level
        current_level_ids = next_level_ids;
    }
    
    // Normalize the result if we have non-zero weight
    if (total_weight > 1e-12) {
        averaged_normal /= total_weight;
        
        // Normalize the vector
        double norm = norm_2(averaged_normal);
        if (norm > 1e-12) {
            averaged_normal /= norm;
        }
    } else {
        // If no weights, return original normal
        averaged_normal = elem_normal;
    }
    
    return averaged_normal;
}

/**
 * @brief Compute and store averaged normals for all cut elements in the model part
 * @param rModelPart The model part containing elements
 * @param NeighborLevels Number of neighbor levels to consider
 * @param VariableName The variable name to store the averaged normal (default: ELEMENT_CUT_NORMAL_AVERAGED)
 * @return Number of elements processed
 */
int ComputeAndStoreAveragedNormals(
    ModelPart& rModelPart,
    int NeighborLevels,
    const std::string& VariableName)
{
    ClearAveragedNormals(rModelPart, VariableName);
    
    // Make sure we have cut normals set
    if (InterfaceAveragesUtility::GetInterfaceAverages().empty()) {
        InterfaceAveragesUtility::ComputeModelPartInterfaceAverages(rModelPart);
    }
    
    int num_elements = SetElementCutNormals(rModelPart);
    std::cout << "Set cut normals for " << num_elements << " elements" << std::endl;
    
    // Count of processed elements
    int count = 0;
    
    // Process each element
    for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it) {
        Element& r_element = *it;
        double length = r_element.GetValue(INTERSECTION_LENGTH_2D);
        
        // Only process cut elements
        if (length > 0.0) {
            // Compute averaged normal
            array_1d<double, 3> averaged_normal = ComputeAveragedElementNormal(
                rModelPart, r_element.Id(), NeighborLevels);
            
            // Store it
            r_element.SetValue(ELEMENT_CUT_NORMAL_AVERAGED, averaged_normal);
            count++;
        }
    }
    
    std::cout << "Computed averaged normals for " << count << " elements using " 
              << NeighborLevels << " neighbor levels" << std::endl;
    
    return count;
}
/**
 * @brief Clear averaged normal values from all elements
 * @param rModelPart The model part containing elements
 * @param VariableName The variable name storing the averaged normal
 */
void ClearAveragedNormals(
    ModelPart& rModelPart,
    const std::string& VariableName)
{
    // Check if the variable exists
    if (!KratosComponents<Variable<array_1d<double, 3>>>::Has(VariableName)) {
        std::cout << "Warning: Variable " << VariableName << " not found. Nothing to clear." << std::endl;
        return;
    }
    
    // Zero vector to reset values
    array_1d<double, 3> zero_vector = ZeroVector(3);
    
    // Clear all element values
    for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it) {
        it->SetValue(ELEMENT_CUT_NORMAL_AVERAGED, zero_vector);
    }
    
    std::cout << "Cleared " << VariableName << " values for all elements" << std::endl;
}

/**
 * @brief Save the averaged normals to file
 * @param rModelPart The model part containing elements
 * @param Filename Filename to save the results
 * @param VariableName Name of the variable storing the averaged normals
 */
void SaveAveragedNormalsToFile(
    const ModelPart& rModelPart,
    const std::string& Filename,
    const std::string& VariableName)
{
    // Open output file
    std::ofstream outFile(Filename);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << Filename << " for writing." << std::endl;
        return;
    }
    
    // Write header
    outFile << std::fixed << std::setprecision(15);
    outFile << "Element_ID\tNormal_X\tNormal_Y\tNormal_Z\tLength" << std::endl;
    
    // Count of elements with valid averaged normals
    int count = 0;
    
    // Process each element
    for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it) {
        Element& r_element = *it;
        
        // Check if element has an intersection length
        double intersection_length = r_element.GetValue(INTERSECTION_LENGTH_2D);
        
        // Skip elements that are not cut by the interface
        if (intersection_length <= 0.0) {
            continue;
        }
        
        // Get the averaged normal vector - directly use ELEMENT_CUT_NORMAL_AVERAGED
        array_1d<double, 3> averaged_normal = r_element.GetValue(ELEMENT_CUT_NORMAL_AVERAGED);
        
        // Skip elements with zero normal (not calculated or invalid)
        if (norm_2(averaged_normal) < 1e-10) {
            continue;
        }
        
        // Write to file
        outFile << std::fixed << std::setprecision(15); 
        outFile << r_element.Id() << "\t"
                << averaged_normal[0] << "\t"
                << averaged_normal[1] << "\t"
                << averaged_normal[2] << "\t"
                << intersection_length << std::endl;
        
        count++;
    }
    
    outFile.close();
    std::cout << "Successfully wrote " << count << " elements with averaged normals to " << Filename << std::endl;
}
}
}