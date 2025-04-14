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
 * @class SurrogateBoundaryModeler
 * @ingroup KratosCore
 * @brief Generates the voxel mesh plus stores the min distance from every node to the skin of the model.
 */
class SurrogateBoundaryModeler : public VoxelMeshGeneratorModeler
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

    /// Pointer definition of SurrogateBoundaryModeler
    KRATOS_CLASS_POINTER_DEFINITION( SurrogateBoundaryModeler );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    SurrogateBoundaryModeler(){}

    //Constructor
    SurrogateBoundaryModeler(
        Model& rModel, Parameters rParameters = Parameters()) 
        : VoxelMeshGeneratorModeler(rModel, rParameters ) {

            // bool cells_in_touch = false;
            // for(auto parameters : rParameters["coloring_settings_list"]){
            //     cells_in_touch = cells_in_touch || (parameters["type"].GetString() == "cells_in_touch");
            //     mInsideColor = parameters["color"].GetDouble();
            // }

            // KRATOS_ERROR_IF(!cells_in_touch) 
            //     << "Surrogate Boundary only works properly including cells_in_touch in the colouring settings" << std::endl;
            
            // //FindVectorDistanceToSurface();

            mInsideColor = 14;
        }

    /// Destructor
    virtual ~SurrogateBoundaryModeler(){}

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * 
     */
    class SurrogateBoundaryNode 
    {
        bool mIsActive = false; 
        double mSignedDistance = 0;
        array_1d<double,3> mDistanceToSkin{0,0,0};
        Node::Pointer node_pointer = nullptr;

        public: 
        bool IsActive() { return mIsActive; }
        void SetActive(bool active) { mIsActive = active; }
        double GetSignedDistance() { return mSignedDistance; }
        void SetSignedDistance(double new_distance) { mSignedDistance = new_distance; }
        bool IsInside() { return mSignedDistance >= 0; }
        array_1d<double,3>& GetVectorDistance() { return mDistanceToSkin; }
        Node::Pointer GetNodePtr(){ return node_pointer; }
        void SetNodePointer(Node::Pointer pNode) { node_pointer = pNode; }
    };

    void ComputeSurrogateBoundary() 
    {
        GeometricalObjectsBins triangle_bin(mpInputModelPart->ElementsBegin(),mpInputModelPart->ElementsEnd());

        array_1d<std::size_t, 3> number_of_divisions = mMeshingData.GetNumberOfDivisions();
        mSurrogateBoundaryData.resize(number_of_divisions[0]*number_of_divisions[1]*number_of_divisions[2]);

        for(auto& r_element : mpInputModelPart->Elements()) {
            Element::GeometryType& r_geometry = r_element.GetGeometry();
            mColors.AddGeometry(r_geometry, false);
        }
        array_1d< std::size_t, 3 > min_ray_position{0,0,0};
        mColors.InitializeRays(min_ray_position, number_of_divisions, "nodes");
        mColors.CalculateNodalRayColors(min_ray_position, number_of_divisions, mInsideColor, 0);

        for (std::size_t i = 0; i < number_of_divisions[0]; i++) 
        {
            for (std::size_t j = 0; j < number_of_divisions[1]; j++) 
            {
                for (std::size_t k = 0; k < number_of_divisions[2]; k++) 
                {
                    CartesianNodalData& nodal_data = mMeshingData.GetNodalData(i,j,k);
                    Node::Pointer node_pointer = nodal_data.pGetNode();
                    if (node_pointer) 
                    {
                        int node_index = mMeshingData.GetNodeIndex(i,j,k);
                        mSurrogateBoundaryData[node_index].SetActive(true);
                        mSurrogateBoundaryData[node_index].SetNodePointer(node_pointer);

                        Point point = *node_pointer;
                        const auto& search_result = triangle_bin.SearchNearest(point);

                        if(search_result.IsObjectFound()) // It should always be found
                        {
                            auto p_result = search_result.Get().get();
                            GeometryPtrType geometry = p_result->pGetGeometry();
                            const PointsArrayType& points = geometry->Points();
                            Point closest_point = NearestPointUtilities::TriangleNearestPoint(point, points[0],points[1],points[2]);
                            
                            for (std::size_t ii = 0; ii < 3; ii++) 
                            {
                                mSurrogateBoundaryData[node_index].GetVectorDistance()[ii] = point[ii] - closest_point[ii];
                            }
                        }   
                        
                        double d = norm_2(mSurrogateBoundaryData[node_index].GetVectorDistance());
                        if (mColors.GetNodalColor(i,j,k) == mInsideColor) 
                        {
                            mSurrogateBoundaryData[node_index].SetSignedDistance( d);
                        } else {
                            mSurrogateBoundaryData[node_index].SetSignedDistance(-d);
                        }

                        KRATOS_INFO("ModelerDistances") << "Point: " << point << " Signed distance to closest point: " << mSurrogateBoundaryData[node_index].GetSignedDistance()  << std::endl; 
                    }
                } 
            }  
        }
    }

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
                    if (node_pointer) 
                    {
                        Point point = *node_pointer;
                        const auto& search_result = triangle_bin.SearchNearest(point);

                        if(search_result.IsObjectFound()) // It should always be found
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
                    }
                } 
            }  
        }
    }

    void ApplyColoringToNodes()
    {
        //auto parameters = GetParameters();
        //const int inside_color = parameters["color"].GetInt();
        //const int outside_color = parameters["outside_color"].GetInt();
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
                mColors.AddGeometry(r_geometry, false);
            }
        } 

        mColors.CalculateNodalRayColors(min_ray_position, number_of_divisions, inside_color, outside_color);
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
                }
            }
        }
    }

    std::vector<SurrogateBoundaryNode>& GetSurrogateBoundaryData() 
    {
        return mSurrogateBoundaryData;
    }

    SurrogateBoundaryNode& GetSurrogateBoundaryNode(std::size_t I, std::size_t J, std::size_t K) 
    {
        return mSurrogateBoundaryData[mMeshingData.GetNodeIndex(I,J,K)];
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

    /// The data of every node in the voxel mesh, with its corresponging SB information
    std::vector<SurrogateBoundaryNode> mSurrogateBoundaryData;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

}; /* Class SurrogateBoundaryModeler */

}  /* namespace Kratos.*/