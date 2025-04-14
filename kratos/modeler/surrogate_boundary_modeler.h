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
 * @brief Generates the voxel mesh plus stores the vector and signed distance from every node to the skin of the model.
 */
class KRATOS_API(KRATOS_CORE) SurrogateBoundaryModeler 
    : public VoxelMeshGeneratorModeler
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

    /// Pointer definition of SurrogateBoundaryModeler
    KRATOS_CLASS_POINTER_DEFINITION( SurrogateBoundaryModeler );

    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * @brief Default constructor
     */
    SurrogateBoundaryModeler(){}

    //Constructor
    SurrogateBoundaryModeler(
        Model& rModel, Parameters rParameters = Parameters()) 
        : VoxelMeshGeneratorModeler(rModel, rParameters ) {
        }

    /// Destructor
    virtual ~SurrogateBoundaryModeler(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Data stored for every voxel in the mesh. Nodes that are not created by the Modeler
     * are marked as non active and pointer to node will be null. Active nodes will have 
     * a vector distance to skin and a signed distance (magnitude of vector distance with
     * positive sign if inside the volume and negative sign otherwise)
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

    /**
     * Computes the distance (vector and signed) for every node in the voxel mesher.
     * @note This method should only be called after calling SetupGeometryModel() and SetupModelPart().
     */
    void ComputeSurrogateBoundary() 
    {
        int inside_color = 14;
        int outside_color = 1;

        GeometricalObjectsBins triangle_bin(mpInputModelPart->ElementsBegin(),mpInputModelPart->ElementsEnd());

        array_1d<std::size_t, 3> number_of_divisions = mMeshingData.GetNumberOfDivisions();
        mSurrogateBoundaryData.resize(number_of_divisions[0]*number_of_divisions[1]*number_of_divisions[2]);

        array_1d< std::size_t, 3 > min_ray_position{0,0,0};
        mColors.InitializeRays(min_ray_position, number_of_divisions, "nodes");
        for(auto& r_element : mpInputModelPart->Elements()) {
            Element::GeometryType& r_geometry = r_element.GetGeometry();
            mColors.AddGeometry(r_geometry, false);
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
                        if (mColors.GetNodalColor(i,j,k) == inside_color) 
                        {
                            mSurrogateBoundaryData[node_index].SetSignedDistance( d);
                        } else {
                            mSurrogateBoundaryData[node_index].SetSignedDistance(-d);
                        }
                    }
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