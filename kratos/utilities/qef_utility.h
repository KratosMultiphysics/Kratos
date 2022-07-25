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
#include "includes/element.h"
#include "includes/geometrical_object.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_3d_3.h"
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
 * @class QEF (quadratic error function)
 * @ingroup KratosCore
 * @brief Utilities to compute the minimum error point in a 3D voxel intersected by a triangle mesh
 * @note the methods in this class is templated for both parameters (voxel and triangle mesh) but can't 
 * actually be used for any type of mesh since the methods in Intersection utilities are not templated. 
 * This methods should be reimplemented in order to use this class as template
 * @author Ariadna Cortés
 */
class QEF
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
    static array_1d<double,3> QEF_point (
        const TGeometryType& rVoxel,  
        const TGeometryArrayType& rTriangles     
    ) {
        array_1d<double,3> xPoint;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        //Initialize the corresponding matrixes
        int n = edges.size();

        for (int i = 0; i < n; i++) {
            Distances.push_back(0);
            PointsArrayType ends = edges[i].Points();
            //std::cout << "Edge " << i << " has nodes " << ends[0] << " " << ends[1] << std::endl;

            for (int i = 0; i < rTriangles.size(); i++) {
                array_1d<double,3> intersection;
                int result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],intersection);
                
                if(result == 1) {
                    array_1d<double,3> normal = CalculateNormal<Triangle3D3<NodeType>>(rTriangles[i]);
                }  
            } 
                        
        }
        return xPoint;
    }

    /**
     * @brief Aproximates the portion of the edge that represents volume
     * @param rDistances references to a sorted vector containing the distances of each intersecting point with the edge
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated volume 
     */  
    static array_1d<double,3> CalculateCenter(GeometryType& rVoxel) {
        PointsArrayType nodes = rVoxel.Points();
        double x = (nodes[0].X() + nodes[1].X())/2.0;
        double y = (nodes[1].Y() + nodes[2].Y())/2.0;
        double z = (nodes[0].Z() + nodes[4].Z())/2.0;
        array_1d<double,3> center({x,y,z});
        return center;
    }

    /**
     * @brief Aproximates the portion of the edge that represents volume
     * @param rDistances references to a sorted vector containing the distances of each intersecting point with the edge
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated volume 
     */  
    static array_1d<double,3> CalculateNormal(GeometryType& triangle) {
        PointsArrayType nodes = triangle.Points();
        array_1d<double,3> u = nodes[1] - nodes[0];
        array_1d<double,3> v = nodes[2] - nodes[0];
        array_1d<double,3> normal;
        MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(normal,u,v);
        return normal;
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
     * @brief Returns the distance between two 3D points.
     * @param rPoint0 reference to the first point
     * @param rPoint1 reference to the second point
     * @return Distance 
     */  
    static double Distance(const NodeType& Point0, const NodeType& Point1) {
        const double lx = Point0.X() - Point1.X();
        const double ly = Point0.Y() - Point1.Y();
        const double lz = Point0.Z() - Point1.Z();

        const double length = lx * lx + ly * ly + lz * lz;

        return std::sqrt( length );
    }
    /**
     * @brief Returns the distance between two 3D points.
     * @param rPoint0 reference to the first point
     * @param rPoint1 reference an array of 3 coordinates representing the second point
     * @return Distance 
     */  
    static double Distance(const NodeType& Point0, const array_1d<double,3>& Point1) {
        const double lx = Point0.X() - Point1[0];
        const double ly = Point0.Y() - Point1[1];
        const double lz = Point0.Z() - Point1[2];

        const double length = lx * lx + ly * ly + lz * lz;

        return std::sqrt( length );
    }

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */