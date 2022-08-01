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
#include "volume_inside_voxel_utility.h"

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
class VolumeInsideVoxelQEF : public VolumeInsideVoxelUtility
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
    template<class TGeometryType, class TGeometryArrayType>
    static double SimpleNodesQEFApproximation(
        const TGeometryType& rVoxel,  
        const TGeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        PointsArrayType Nodes = rVoxel.Points();

        array_1d<double,3> qef = QEF::QEFPoint(rVoxel,rTriangles); 
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 
        
        DenseMatrix<unsigned int> NodesInFaces(4,6);
        rVoxel.NodesInFaces(NodesInFaces);

        for(int i = 0; i < NodesInFaces.size2(); i++) {
            double Portion = GetPortion(Nodes,NodesInFaces,i);
            double dist = NormalizedDistanceToQEF(Nodes, NodesInFaces, qef, i);
            
            double PartialVolume = Portion*abs(dist)/3.0;   //Volume of a piramid
            volume += PartialVolume;
        }
        //if (volume == 0) return EdgesPortionApproximation(rVoxel,rTriangles);
        
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
    
    /**
     * @brief Approximates the portion of a face that actually corresponds to volume (assigning each node 
     * 1/numberOfNodes portion if it is inside the volume)
     * @param rFace Reference to the face whos portion must be calculated
     * @return 
     * */
    static double GetPortion(PointsArrayType& Nodes,const DenseMatrix<unsigned int>& NodesInFaces, int& face) 
    {
        double Portion = 0;
        for(int i = 0; i < NodesInFaces.size1(); i++) {
            if( Nodes[NodesInFaces(i,face)].GetSolutionStepValue(DISTANCE) > 0) Portion += 1.0/4;
        }
        return Portion;
    }

    static double NormalizedDistanceToQEF(
        PointsArrayType& Nodes,
        const DenseMatrix<unsigned int>& NodesInFaces, 
        const array_1d<double,3>& Point, int& face) 
    {
        array_1d<double, 3> edge1 = Nodes[NodesInFaces(1,face)] - Nodes[NodesInFaces(0,face)];
        array_1d<double, 3> edge2 = Nodes[NodesInFaces(2,face)] - Nodes[NodesInFaces(0,face)];
        array_1d<double, 3> mNormal;
        MathUtils<double>::UnitCrossProduct(mNormal, edge1, edge2);

        const double mConstant =  inner_prod(mNormal, Nodes[NodesInFaces(0,face)]);
        double Side = Distance(Nodes[NodesInFaces(1,face)],Nodes[NodesInFaces(0,face)]);
        double Distance = inner_prod(mNormal,Point) - mConstant;

        //std:: cout << "Plane Vector 1: " << edge1 <<std::endl << "Plain Vector 2: " << edge2 << std::endl;
        //std::cout << "Normal vector: " << mNormal <<std::endl;

        return Distance/Side;
    }

}; /* Class VoxelInsideVolumeQEF */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */