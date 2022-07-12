//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//
//

#if !defined(KRATOS_VOLUME_INSIDE_VOXEL)
#define  KRATOS_VOLUME_INSIDE_VOXEL

// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/hexahedra_3d_8.h"
#include "intersection_utilities.h"

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
 * @class VolumeInsideVoxel
 * @ingroup KratosCore
 * @brief Utilities to compute the real volume inside a voxel
 * @details This class provides static methods to compute (using different approximations) the portion of a 
 * voxel that is actually filled with volume.
 * @author Ariadna Cortés
 */
class VolumeInsideVoxelUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;


    /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( VolumeInsideVoxelUtility );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    VolumeInsideVoxelUtility(){}

    /// Destructor
    virtual ~VolumeInsideVoxelUtility(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @return Approximated volume 
     * @note This approximation assigns a fraction of volume (1/8) to each node of the
     * voxel, and counts it as volume if the node is inside the object (NodeDistance > 0)
     */  
    template<class TGeometryType>
    static double NodesApproximation(
        const TGeometryType& rVoxel        
    ) {
        double volume = 0;
        PointsArrayType nodes = rVoxel.Points();
        for (int i = 0; i < 8; i++) {
            if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=0.125; //heyho
            } 
        }
        return volume;
    }

    /*This method is completly useless since it does the same calculation as the previous 
    one but in a different way. Helps to illustrate use of edges
    */
    template<class TGeometryType>
    static double EdgesApproximation(
        const TGeometryType& rVoxel        
    ) {
        double volume = 0;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        PointsArrayType nodes = rVoxel.Points();
        for (int i = 0; i < 12; i++) {
            PointsArrayType ends = edges[i].Points();
            if(ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=1.0/12;
            } else if(
                ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) < 0 || 
                ends[0].GetSolutionStepValue(DISTANCE) < 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0 ) {
                volume+=1.0/24;
            }
        }
        return volume;
    }

    template<class TGeometryType, class TGeometryArrayType>
    static double EdgesPortionApproximation(
        const TGeometryType& rVoxel,  
        const TGeometryArrayType& rTriangles     
    ) {
        //still not implemented
        double volume = 0;
        return volume;
    }
    
private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

}; /* Class VoxelInsideVolume */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */