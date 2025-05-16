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

#pragma once

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/geometry.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"
#include "voxel_mesher_operation.h"

namespace Kratos {

class ComputeSurrogateBoundaryData: public VoxelMesherOperation {

    typedef Node NodeType;
    typedef Node::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::PointsArrayType PointsArrayType;

public:
    ComputeSurrogateBoundaryData(VoxelMeshGeneratorModeler& rModeler, Parameters OperationParameters):
        VoxelMesherOperation(rModeler, OperationParameters)
    {}

    ~ComputeSurrogateBoundaryData() override = default ;

    Parameters GetDefaultParameters() const override;

    void ValidateParameters() override;

    void Apply() const override;

private:
     /**
     * @class SurrogateBoundaryNode
     * @brief Data for every voxel in the mesh. Nodes that are not created by the Modeler
     * are marked as non active and pointer to node will be null. Active nodes will have 
     * a vector distance to skin and a signed distance (magnitude of vector distance with
     * positive sign if inside the volume and negative sign otherwise).
     */
    class SurrogateBoundaryNode 
    {
        bool mIsActive = false; 
        bool mIsInside = false;
        double mSignedDistance = 0;
        array_1d<double,3> mDistanceToSkin{0,0,0};
        Node::Pointer node_pointer = nullptr;

     public: 
        bool& IsActive() { return mIsActive; }
        bool& IsInside() { return mIsInside; }
        double& GetSignedDistance() { return mSignedDistance; }
        array_1d<double,3>& GetVectorDistance() { return mDistanceToSkin; }
        
        // Pointer to original node in VoxelMesher
        Node::Pointer GetNodePtr(){ return node_pointer; }
        void SetNodePointer(Node::Pointer pNode) { node_pointer = pNode; }
    };
   

    // std::vector<SurrogateBoundaryNode>& GetSurrogateBoundaryNodes() 
    // {
    //     return mSurrogateBoundaryData;
    // }

    // SurrogateBoundaryNode& GetSurrogateBoundaryNode(std::size_t I, std::size_t J, std::size_t K) 
    // {
    //     return mSurrogateBoundaryData[GetNodeIndex(I,J,K)];
    // }

    std::string PrintSBData(std::vector<SurrogateBoundaryNode>& sbdata) const;

};
}