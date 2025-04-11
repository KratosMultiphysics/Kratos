//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes Danes
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/geometry.h"
#include "voxel_mesh_generator_modeler.h"

#include "spatial_containers/geometrical_objects_bins.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"

namespace Kratos
{ 
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class VoxelMeshGeneratorModelerWithMinDist
 * @ingroup KratosCore
 * @brief Generates the voxel mesh plus stores the min distance from every node to the skin of the model.
 */
class VoxelMeshGeneratorModelerWithMinDist : public VoxelMeshGeneratorModeler
{
public:

    ///@name Type Definitions
    ///@{

    typedef Node NodeType;
    typedef Node::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    /// Pointer definition of VoxelMeshGeneratorModelerWithMinDist
    KRATOS_CLASS_POINTER_DEFINITION( VoxelMeshGeneratorModelerWithMinDist );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
     VoxelMeshGeneratorModelerWithMinDist(){}

    //Constructor
    VoxelMeshGeneratorModelerWithMinDist(
        Model& rModel, Parameters rParameters = Parameters()) 
        : VoxelMeshGeneratorModeler(rModel, rParameters ) {

            // bool cells_in_touch = false;
            // for(auto parameters : rParameters["coloring_settings_list"]){
            //     cells_in_touch = cells_in_touch || (parameters["type"].GetString() == "cells_in_touch");
            //     mInsideColor = parameters["color"].GetDouble();
            // }

            // KRATOS_ERROR_IF(!cells_in_touch) 
            //     << "Min Distance Function only works including cells_in_touch in the colouring settings" << std::endl;
            
            // //FindVectorDistanceToSurface();

            mInsideColor = 14;
        }

    /// Destructor
    virtual ~VoxelMeshGeneratorModelerWithMinDist(){}

    ///@}
    ///@name Operations
    ///@{


    void FindDistanceToSkin() {
        GeometricalObjectsBins triangle_bin(mpInputModelPart->ElementsBegin(),mpInputModelPart->ElementsEnd());

        array_1d<std::size_t, 3> number_of_divisions = mMeshingData.GetNumberOfDivisions();
        mDistances.resize(number_of_divisions[0]*number_of_divisions[1]*number_of_divisions[2]);

        for (std::size_t i = 0; i < number_of_divisions[0]; i++) 
        {
            for (std::size_t j = 0; j < number_of_divisions[1]; j++) 
            {
                for (std::size_t k = 0; k < number_of_divisions[2]; k++) 
                {
                    CartesianNodalData& nodal_data = mMeshingData.GetNodalData(i,j,k);
                    Node::Pointer node_pointer = nodal_data.pGetNode();
                    int node_index = mMeshingData.GetNodeIndex(i,j,k);
                    if (node_pointer)  //falla cada cop que k = 0 xd;
                    {
                        Point point = *node_pointer;
                        const auto& search_result = triangle_bin.SearchNearest(point);

                        if(search_result.IsObjectFound()) 
                        {
                            auto p_result = search_result.Get().get();
                            GeometryPtrType geometry = p_result->pGetGeometry();
                            const PointsArrayType& points = geometry->Points();
                            Point closest_point = NearestPointUtilities::TriangleNearestPoint(point, points[0],points[1],points[2]);
                            
                            for (std::size_t ii = 0; ii < 3; ii++) 
                            {
                                mDistances[node_index][ii] = point[ii] - closest_point[ii];
                            }

                            KRATOS_INFO("ModelerDistances") << "Point: " << point << " Distance to closest point: " << mDistances[node_index] << std::endl; 
                        }   
                    } else std::cout << "Problemas..." << std::endl;
                } 
            }  
        }
        std::cout << number_of_divisions << std::endl; 
    }

    void ApplyColoringToNodes()
    {
        /* auto parameters = GetParameters();
        auto& r_model_part = GetModelPart(parameters["model_part_name"].GetString());
        auto& r_colors = GetMeshColors();
        const int inside_color = parameters["color"].GetInt();
        const int outside_color = parameters["outside_color"].GetInt();
        const std::string input_entities = parameters["input_entities"].GetString(); */
        array_1d<std::size_t, 3> number_of_divisions = mMeshingData.GetNumberOfDivisions();
        mIsInside.resize(number_of_divisions[0]*number_of_divisions[1]*number_of_divisions[2]);
        const int inside_color = mInsideColor;
        const int outside_color = 1;
        const std::string input_entities = "elements";

        array_1d< std::size_t, 3 > min_ray_position{0,0,0};
        mColors.InitializeRays(min_ray_position, number_of_divisions, "nodes");

        if(input_entities == "elements") {
            for(auto& r_element : mpInputModelPart->Elements()) {
                Element::GeometryType& r_geometry = r_element.GetGeometry();
                //CheckGeometryType(r_geometry.GetGeometryType());
                mColors.AddGeometry(r_geometry, false);
            }
        } /* else if(input_entities == "conditions"){
            for(auto& r_condition : r_model_part.Conditions()) {
                Condition::GeometryType& r_geometry = r_condition.GetGeometry();
                CheckGeometryType(r_geometry.GetGeometryType());
                r_colors.AddGeometry(r_geometry, false);
            }
        } else if(input_entities == "geometries"){
            for(auto& r_geometry : r_model_part.Geometries()) {
                CheckGeometryType(r_geometry.GetGeometryType());
                r_colors.AddGeometry(r_geometry, false);
            }
        } else{
            KRATOS_ERROR << "The input_entities " << parameters["input_entities"] << " is not supported. The supported input_entities are geometries, elements and conditions" << std::endl;
        } */

        mColors.CalculateNodalRayColors(min_ray_position, number_of_divisions, inside_color, outside_color);
        int nodes_inside = 0;
        for (std::size_t i = 0; i < number_of_divisions[0]; i++) 
        {
            for (std::size_t j = 0; j < number_of_divisions[1]; j++) 
            {
                for (std::size_t k = 0; k < number_of_divisions[2]; k++) 
                {
                    CartesianNodalData& nodal_data = mMeshingData.GetNodalData(i,j,k);
                    Node::Pointer node_pointer = nodal_data.pGetNode();
                    int node_index = mMeshingData.GetNodeIndex(i,j,k);
                    mIsInside[node_index] = mColors.GetNodalColor(i,j,k) == mInsideColor? 1:0; 
                    if (node_pointer)  //falla cada cop que k = 0 xd;
                    {
                        Point point = *node_pointer;
                        KRATOS_INFO("ColorNodes") << "Point: " << point << " \n Color: " << mColors.GetNodalColor(i,j,k) << std::endl;  
                        if(mColors.GetNodalColor(i,j,k) == mInsideColor) nodes_inside++;
                    } else std::cout << "DRAAAMAA" << " Color: " << mColors.GetNodalColor(i,j,k) << std::endl;
                }
            }
        }
        std::cout << "Nodes Inside: " << nodes_inside << std::endl; 
    }

    std::vector<int> &GetNodalColors() 
    { 
        return mIsInside; 
    }

    std::vector<array_1d<double,3>> &GetDistances()
    {
        return mDistances;
    }

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

        /// The color assigned to the nodes inside the volume
        double mInsideColor;

        /// @brief For each node, 1 if it is inside the volume, 0 if outside.
        std::vector<int> mIsInside; //bools

        /// @brief Distance from each node to the closest skin point.
        std::vector<array_1d<double,3>> mDistances;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /* bool CheckAllNeighbourQEFAreInside(int i, int j, int k) {
        bool res = (mColors.GetElementalColor(i-1,j-1,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j-1,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i-1,j,k-1) == mInsideColor) &&
                    (mColors.GetElementalColor(i-1,j-1,k) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j-1,k) == mInsideColor) &&
                    (mColors.GetElementalColor(i,j,k) == mInsideColor) &&
                    (mColors.GetElementalColor(i-1,j,k) == mInsideColor);
        return res;
    } */

}; /* Class DualCountouringMesher */

}  /* namespace Kratos.*/