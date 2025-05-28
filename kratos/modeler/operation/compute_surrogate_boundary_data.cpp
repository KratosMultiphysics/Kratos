//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes Danes
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "compute_surrogate_boundary_data.h"

namespace Kratos {

Parameters ComputeSurrogateBoundaryData::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "compute_surrogate_boundary_data",
        "model_part_name": "Undefined",
        "outside_color": 1,
        "inside_color": -1,
        "output_path" : "distances.txt"
    })");
}

void ComputeSurrogateBoundaryData::ValidateParameters()
{
    Parameters parameters = GetParameters();
    parameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

void ComputeSurrogateBoundaryData::Apply() const
{   
    /// The data of every node in the voxel mesh, with its corresponging SB information
    std::vector<SurrogateBoundaryNode> surrogate_boundary_data;

    Parameters parameters = GetParameters();
    int outside_color = parameters["outside_color"].GetInt();
    int inside_color = parameters["inside_color"].GetInt();
    if (outside_color == inside_color) inside_color = outside_color == 0? -1 : -outside_color;
    ModelPart& r_skin_part = GetModelPart(parameters["model_part_name"].GetString());

    std::string output_path = parameters["output_path"].GetString();

    GeometricalObjectsBins elements_bin(r_skin_part.ElementsBegin(),r_skin_part.ElementsEnd());
    GeometricalObjectsBins conditions_bin(r_skin_part.ConditionsBegin(),r_skin_part.ConditionsEnd());
    array_1d<std::size_t, 3> number_of_divisions = GetNumberOfDivisions();
    surrogate_boundary_data.resize(number_of_divisions[0]*number_of_divisions[1]*number_of_divisions[2]);

    // Compute colors for every node. 
    array_1d< std::size_t, 3 > min_ray_position{0,0,0};
    auto colors = GetMeshColors();
    colors.InitializeRays(min_ray_position, number_of_divisions, "nodes");
    for(auto& r_element : r_skin_part.Elements()) {
        auto& r_geometry = r_element.GetGeometry();
        colors.AddGeometry(r_geometry, false);
    }
    for(auto& r_condition : r_skin_part.Conditions()) {
        auto& r_geometry = r_condition.GetGeometry();
        colors.AddGeometry(r_geometry, false);
    }
    colors.CalculateNodalRayColors(min_ray_position, number_of_divisions, inside_color, outside_color);

    // Compute distance vector and signed distance for every node.
    for (std::size_t i = 0; i < number_of_divisions[0]; i++) 
    {
        for (std::size_t j = 0; j < number_of_divisions[1]; j++) 
        {
            for (std::size_t k = 0; k < number_of_divisions[2]; k++) 
            {
                const std::size_t node_index = GetNodeIndex(i,j,k);
                if (colors.GetNodalColor(i,j,k) == inside_color) 
                {
                    surrogate_boundary_data[node_index].IsInside() = true;
                }

                auto& nodal_data = GetNodalData(i,j,k);
                auto node_pointer = nodal_data.pGetNode();
                if (node_pointer) 
                {
                    surrogate_boundary_data[node_index].IsActive() = true;
                    surrogate_boundary_data[node_index].SetNodePointer(node_pointer);

                    Point point = *node_pointer;
                    GeometricalObjectsBins::ResultType search_result = elements_bin.SearchNearest(point);
                    if (!search_result.IsObjectFound()) 
                    {
                        search_result = conditions_bin.SearchNearest(point);
                    }
                    if(search_result.IsObjectFound()) // It should always be found
                    {
                        auto p_result = search_result.Get().get();
                        GeometryPtrType geometry = p_result->pGetGeometry();
                        const PointsArrayType& points = geometry->Points();
                        Point closest_point = NearestPointUtilities::TriangleNearestPoint(point, points[0],points[1],points[2]);
                        
                        for (std::size_t ii = 0; ii < 3; ii++) 
                        {
                            surrogate_boundary_data[node_index].GetVectorDistance()[ii] = point[ii] - closest_point[ii];
                        }
                    } else {
                        KRATOS_WARNING("SurrogateBoundaryModeler") << "Input geometry has no elements or conditions. Unable to compute distance to skin." << std::endl;
                    } 
                    
                    double d = norm_2(surrogate_boundary_data[node_index].GetVectorDistance());
                    if (colors.GetNodalColor(i,j,k) == inside_color) 
                    {
                        surrogate_boundary_data[node_index].GetSignedDistance() = -d;
                    } else {
                        surrogate_boundary_data[node_index].GetSignedDistance() = d;
                    }
                }
            } 
        }  
    }

    std::ofstream output_file(output_path);
    
    // Check if the file was opened successfully
    if (output_file.is_open()) {
        output_file <<  PrintSurrogateBoundaryData(surrogate_boundary_data);
        output_file.close();
    } else  
        KRATOS_ERROR << "ComputeSurrogateBoundaryData: Data could not be saved in file" << std::endl;
} 

/// Print object's data.
std::string ComputeSurrogateBoundaryData::PrintSurrogateBoundaryData(std::vector<SurrogateBoundaryNode>& rSurrogateBoundaryData) const 
{
    std::stringstream rOStream;
    const auto n = GetNumberOfDivisions();

    rOStream << "=== Meshing Data ===" << std::endl;
    rOStream << "NumberOfDivisions: " << n[0] << " " << n[1] << " " << n[2] << std::endl;

    rOStream << "VoxelSize: " 
            << GetKeyPlanes(0)[1] - GetKeyPlanes(0)[0] << " " 
            << GetKeyPlanes(1)[1] - GetKeyPlanes(1)[0] << " " 
            << GetKeyPlanes(2)[1] - GetKeyPlanes(2)[0] << " " 
            << std::endl;

    rOStream << "=== Node Data ===" << std::endl;

    for (std::size_t i = 0; i < n[0]; ++i) {
        for (std::size_t j = 0; j < n[1]; ++j) {
            for (std::size_t k = 0; k < n[2]; ++k) {
                SurrogateBoundaryNode& node = rSurrogateBoundaryData[GetNodeIndex(i,j,k)];;
                auto& v = node.GetVectorDistance();
                rOStream << node.IsActive() << " " << node.IsInside() << " "
                        << v[0] << " " << v[1] << " " << v[2] << std::endl;
            }
        }
    }
    return rOStream.str();
}

}