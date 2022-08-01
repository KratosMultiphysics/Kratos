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

#if !defined(KRATOS_QEF)
#define  KRATOS_QEF

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/ublas_interface.h"
#include "includes/geometrical_object.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_3d_3.h"
#include "intersection_utilities.h"
#include "../external_libraries/a_matrix/include/matrix.h"

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
 * @class QEF (quadratic error function)
 * @ingroup KratosCore
 * @brief Utilities to compute the minimum error point in a 3D voxel intersected by a triangle mesh
 * @note the methods in this class is templated for both parameters (voxel and triangle mesh) but can't 
 * actually be used for any type of mesh since the methods in Intersection utilities are not templated. 
 * This methods should be reimplemented in order to use this class as template
 * @author Ariadna Cortés
 */
class KRATOS_API(KRATOS_CORE) QEF
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
    typedef Matrix MatrixType;
    typedef Vector VectorType;


    /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( QEF );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    QEF(){}

    /// Destructor
    virtual ~QEF(){}

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose x-point will be calculated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return The Point (coordinates) of the x-point of the voxel
     */  
    template<class TGeometryType, class TGeometryArrayType>
    static array_1d<double,3> QEFPoint (
        const TGeometryType& rVoxel,  
        const TGeometryArrayType& rTriangles     
    ) {
        array_1d<double,3> xPoint;
        array_1d<double,3> center = CalculateCenter(rVoxel);
        VectorType vCenter = center; 
        MatrixType mCenter(3,1);
        column(mCenter,0) = center;
        GeometryArrayType edges = rVoxel.GenerateEdges();

        //Initialize the corresponding matrixes
        MatrixType AtA(3,3,0);  //3x3 matrix initialized to 0
        MatrixType AtB(3,1,0);  //3x1 matrix
        MatrixType BtB(1,1,0);  //1x1 matrix

        for (int i = 0; i < rTriangles.size(); i++) {
            array_1d<double,3> normal = CalculateNormal(rTriangles[i]);
            VectorType vNormal = normal; 
            MatrixType mNormal(3,1);
            column(mNormal,0) = normal;

            //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
            //but possible), only one will be taken into account to create the matrixes.
            int result = 0; 
            array_1d<double,3> intersection;
            VectorType vIntersection;
            int j = 0;
            while(!result && j < edges.size()) { 
                PointsArrayType ends = edges[j++].Points();
                result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],intersection);
                vIntersection = intersection;
            }
            if (result) {
                //Fill the matrixes with the corresponding information from the intersection and normal
                MatrixType mNormalTrans = trans(mNormal);
                MatrixType help = prod(mNormal,mNormalTrans);
                AtA = AtA + help;
                double aux = MathUtils<double>::Dot(normal,intersection);
                MatrixType mAux(1,1,aux);
                AtB = AtB + prod(mNormal,mAux);
                BtB = BtB + prod(mAux,mAux); 
            }
        }

        //Find the eigenvalues and eigenvectors to AtA
        MatrixType mEigenvectors;
        MatrixType mEigenvalues;
        const bool converged = MathUtils<double>::GaussSeidelEigenSystem(AtA, mEigenvectors, mEigenvalues);

        MatrixType D(3,3,0);
        for (int i : {0,1,2}) mEigenvalues(i,i) < 1e-12 ? D(i,i) = 0 : D(i,i) = check(1.0/mEigenvalues(i,i), 1e-12);

        MatrixType AtAInverse;  
        MathUtils<double>::BDBtProductOperation(AtAInverse, D, mEigenvectors);
       
        MatrixType AtAc = prod(AtA,mCenter);
        MatrixType solution = prod(AtAInverse, AtB - AtAc) + mCenter;        
        xPoint = column(solution,0);

        return xPoint;
    }

    /**
     * @brief Aproximates the portion of the edge that represents volume
     * @param rDistances references to a sorted vector containing the distances of each intersecting point with the edge
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated volume 
     */  
    static array_1d<double,3> CalculateCenter(const GeometryType& rVoxel);

    /**
     * @brief Aproximates the portion of the edge that represents volume
     * @param rDistances references to a sorted vector containing the distances of each intersecting point with the edge
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated volume 
     */  
    static array_1d<double,3> CalculateNormal(const GeometryType& triangle);
    
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

    static void write(MatrixType& A) {
        for (int i = 0; i < A.size1(); i++) {
            for (int j = 0; j < A.size2(); j++) {
                std::cout << A(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static double check(const double& d, const double& epsilon) {
        return d > epsilon ? d : 0;
    }

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */