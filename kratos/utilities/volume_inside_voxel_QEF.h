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

#if !defined(KRATOS_VOLUME_INSIDE_VOXEL_QEF)
#define  KRATOS_VOLUME_INSIDE_VOXEL_QEF

// System includes

// External includes

// Project includes
#include "qef_utility.h"
#include "voxel_utilities.h"

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
class VolumeInsideVoxelQEF : public VoxelUtilities
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
    KRATOS_CLASS_POINTER_DEFINITION( VolumeInsideVoxelQEF );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    VolumeInsideVoxelQEF(){}

    /// Destructor
    virtual ~VolumeInsideVoxelQEF(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     * @note This approximation finds the portion of each edge that is part of the volume (using
     * intersection point with triangles of the mesh). Even if this class is templated for both 
     * parameters, it will only work with intersecting TRIANGLES, since the utility used to compute
     * the intersection does not allow templating.
     */ 
    static double SimpleNodesQEFApproximation(const GeometryType& rVoxel, const GeometryArrayType& rTriangles);

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     * @note This approximation finds the portion of each edge that is part of the volume (using
     * intersection point with triangles of the mesh).
     */ 
    static double FacesPortionQEFApproximation(const GeometryType& rVoxel, const GeometryArrayType& rTriangles);

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     */
    static double VoxelVolumeQEFApproximation(const GeometryType& rVoxel,const GeometryArrayType& rTriangles);

    /**
     * @brief Aproximates the actual volume inside the hexa (can have shifted angles!) 
     * @param rVoxel references to the hexa whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the hexa at some edge.
     * @return Approximated volume 
     */
    static double HexaVolumeQEFApproximation(const GeometryType& rVoxel, const GeometryArrayType& rTriangles);

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
    
    /**
     * @brief Approximates the portion of a face that actually corresponds to area (assigning each node 
     * 1/numberOfNodes portion if it is inside the volume)
     * @param Nodes The nodes of the geometry
     * @param NodesInFaces matrix containing the index of the nodes of the geometry that belong to each faces
     * @return 
     * */
    static double GetPortion(PointsArrayType& Nodes,const DenseMatrix<unsigned int>& NodesInFaces, int& face);

    /**
     * @brief Calculates the distance from a face (given by its nodes) and the given point, normalized 
     * to the size of the face side.
     * @param Nodes The nodes of the geometry
     * @param NodesInFaces matrix containing the index of the nodes of the geometry that belong to each faces
     * @return Distance from the face to the specified point
     * */
    static double NormalizedDistanceToQEF(
        PointsArrayType& Nodes,
        const DenseMatrix<unsigned int>& NodesInFaces, 
        const array_1d<double,3>& Point, int& face);

    /**
     * @brief Calculates the distance from a face (given by its nodes) and the given point, normalized 
     * to the size of the face side.
     * @param Nodes The nodes of the geometry
     * @param NodesInFaces matrix containing the index of the nodes of the geometry that belong to each faces
     * @return Distance from the face to the specified point
     * */
    static double NormalizedDistanceToQEF(
        GeometryType& rFace, 
        const array_1d<double,3>& Point, int& face);

}; /* Class VoxelInsideVolumeQEF */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */