//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
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
 * @author Ariadna Cortes
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
     * @brief Finds the QEF point of a voxel 
     * @param rVoxel references to the voxel whose x-point will be calculated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return The QEF point (x,y,z) 
     */

    static array_1d<double,3> QEFPoint (
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    );

    /**
     * @brief Calculates the center of a geometry
     * @param rVoxel reference to the geometry
     * @return Center (x,y,z)
     */  
    static array_1d<double,3> CalculateCenter(const GeometryType& rVoxel);

    /**
     * @brief Calculates the normal vector to the surface of a 3D triangle 
     * @param rTriangle reference to the triangle
     * @return Normal vector (x,y,z)
     */  
    static array_1d<double,3> CalculateNormal(const GeometryType& rTriangle);
    
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

    static double Check(
        const double D, 
        const double Epsilon
    );

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */