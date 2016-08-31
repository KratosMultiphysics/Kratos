/*
 * File:   AdvancedNMPointsMapper.hpp
 * Author: jcotela
 * Co-author: vmataix, rzorrilla
 *
 * Created on 19 January 2010, 10:20
 * Last update on 28 August 2016, 10:28
 */

#if !defined(KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED )
#define  KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED

// System includes
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "fsi_application.h"
#include "includes/model_part.h"
#include "containers/array_1d.h"
#include "spatial_containers/spatial_containers.h" // kd-tree
#include "utilities/math_utils.h"                  // Cross Product

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

/** @brief Custom Gauss Point container to be used by the mapper
 */
class GaussPointItem: public Point<3>
{
public:

    ///@name Type Definitions
    ///@{
    /// Auxiliar matrix 3x3 employed for the 3D cases
    typedef boost::numeric::ublas::bounded_matrix<double,3,3> MatrixVar;

    /// Counted pointer of GaussPointItem
    KRATOS_CLASS_POINTER_DEFINITION( GaussPointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    GaussPointItem():
        Point<3>(),
        mArea(0),
        mProjStatus(0)
    {
        mNormal = ZeroVector(3);
    }

    GaussPointItem(
            array_1d<double, 3> Coords,
            double Area,
            array_1d<double, 3> Normal
            ):
        Point<3>(Coords),
        mArea(Area),
        mNormal(Normal),
        mProjStatus(0)
    {}

    ///Copy constructor  (not really required)
    GaussPointItem(const GaussPointItem& rhs):
        Point<3>(rhs),
        mArea(rhs.mArea),
        mNormal(rhs.mNormal),
        mProjStatus(rhs.mProjStatus),
        mpOriginCond(rhs.mpOriginCond),
        mpOriginNode(rhs.mpOriginNode)
    {
        mOriginCoords[0] = rhs.mOriginCoords[0];
        mOriginCoords[1] = rhs.mOriginCoords[1];
    }

    /// Destructor.
    // ~GaussPointItem();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Returns the area of the condtition
     * @return The area of the condition
     */
    void GetArea(double& Area)
    {
        Area = mArea;
    }

    /**
     * Returns the normal of the condtition
     * @return The normal of the condition
     */
    void GetNormal(array_1d<double, 3>& Normal)
    {
        Normal = mNormal;
    }

    /**
     * It returns the distance along normal from Gauss point to a condition
     * @return
     */

    void GetDist(double& Dist)
    {
        Dist = mDist;
    }

    /**
     * It returns the projection status
     * @return Proj: The projection status
     */

    void GetProjStatus(int& Proj)
    {
        Proj = mProjStatus;
    }

    /**
     * Test function
     */

    boost::weak_ptr<Condition> GetOriginCond()
    {
        return mpOriginCond;
    }

    /**
     * It sets a projection for a condition
     */

    void SetProjection(
            Condition::WeakPointer Cond,
            array_1d<double,2> Coords,
            double Dist
            )
    {
        mpOriginCond = Cond;
        mOriginCoords = Coords;
        mDist = Dist;
        mProjStatus = 1;
    }

    /**
     * It sets a projection for a node
     */

    void SetProjection(
            Node<3>::WeakPointer pNode,
            const double SqDist
            )
    {
        mpOriginNode = pNode;
        mDist = SqDist;
        mProjStatus = 2;
        mOriginCoords[0] = 0.0;
        mOriginCoords[1] = 0.0;
    }

    /**
     * It projects in 2D/3D for a line/triangle a returns the local coordinates and distance
     */

    void Project(
            Condition::Pointer pOriginCond,
            array_1d<double,2> & Coords,
            double & Dist,
            const int dimension
            );
            
    /**
     * Project a point over a plane
     * @param PointInPlane: A point in the plane
     * @param PointToBeProjected: The point to be projected
     * @param Normal: The normal of the plane
     * @return PointProjected: The point pojected over the plane
     * @return dist: The distance between the point and the plane
     */

    void ProjectPointToPlane(
            const Point<3> & PointInPlane,
            const Point<3> & PointToBeProjected,
            Point<3> & PointProjected,
            double & dist,
            const array_1d<double,3> & Normal
            );

    /**
     * It gets the projected value for scalar variables
     * @param rOriginVar: The variable (scalar) in the original condition
     * @return Value: The projected value (scalar)
     * @param dimension: 2D/3D case
     */

    void GetProjectedValue(
            const Variable<double> & rOriginVar,
            double& Value,
            const int dimension
            );

    /**
     * It gets the projected value for vector variables
     * @param rOriginVar: The variable (vector) in the original condition
     * @return Value: The projected value (vector)
     * @param dimension: 2D/3D case
     */

    void GetProjectedValue(
            const Variable< array_1d<double,3> > & rOriginVar,
            array_1d<double,3>& Value,
            const int dimension
            );

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    double mArea;                        // Destinationn condition's area
    array_1d<double,3> mNormal;          // Destinationn condition's normal
    double mDist;                        // For GP projected to a Condition, Distance along Normal from Gauss Point to Condition
    int mProjStatus;                     // For GP projected to a Node, SQUARED Distance to Node
                                            // 0: Not Projected
                                            // 1: Projected to a condition
                                            // 2: Couldn't be projected, but a value can be obtained from a nearby node
    Condition::WeakPointer mpOriginCond; // Condition pointer
    array_1d<double,2> mOriginCoords;    // For GP projected to a condition
    Node<3>::WeakPointer mpOriginNode;   // For GP projected to a node

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class GaussPointItem

/** @brief Mapper Class whic allows to project values from a domain to other
 */

class AdvancedNMPointsMapper
{
    ///@name Type Definitions
    ///@{
    // Type definitions for the tree
    typedef GaussPointItem                              PointType;
    typedef GaussPointItem::Pointer                     PointTypePointer;
    typedef std::vector<PointType::Pointer>             GaussPointVector;
    typedef std::vector<PointType::Pointer>::iterator   GaussPointIterator;
    typedef std::vector<double>                         DistanceVector;
    typedef std::vector<double>::iterator               DistanceIterator;

    // KDtree definitions
    typedef Bucket< 3ul, PointType, GaussPointVector, PointTypePointer, GaussPointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > tree;

    // Auxiliar matrix 3x3 for 3D cases
    typedef boost::numeric::ublas::bounded_matrix<double, 3, 3> MatrixVar;

public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    /**
     * Class Constructor
     * WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
     * Use an InterfacePreprocess object to create such a model part from a regular one:
     * InterfaceMapper = InterfacePreprocess()
     * InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     * @param rOriginModelPart: The original model part
     * @return rDestinationModelPart: The destination model part
     */
    AdvancedNMPointsMapper(
            const ModelPart & rOriginModelPart,
            ModelPart & rDestinationModelPart
            );

    /// Destructor.
    //~AdvancedNMPointsMapper();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * It maps a scalar variable to a normal vector from a model part to other
     * @param rOriginVar: The original value (scalar) of the variable
     * @return rDestVar: The variable (normal vector) in the destiny modelpart
     * @return MaxIter: Maximum number of iteration allowed
     * @return TolIter: Tolerance accepted in the iteration
     * @return sign_pos: Positive or negative projection
     */

    void ScalarToNormalVectorMap(
            const Variable<double> & rOriginVar,
            Variable<array_1d<double,3> >& rDestVar,
            const int MaxIter,
            const double TolIter,
            const bool sign_pos
            );

    /**
     * It maps a normal vector variable to a scalar from a model part to other
     * @param rOriginVar: The original value (normal vector) of the variable
     * @return rDestVar: The variable (scalar) in the destiny modelpart
     * @return MaxIter: Maximum number of iteration allowed
     * @return TolIter: Tolerance accepted in the iteration
     * @return sign_pos: Positive or negative projection
     */

    void NormalVectorToScalarMap(
            const Variable<array_1d<double,3> >& rOriginVar,
            Variable<double> & rDestVar,
            const int MaxIter,
            const double TolIter,
            const bool sign_pos
            );

    /**
     * It maps a variable (scalar) from a model part to other
     * @param rOriginVar: The original value of the variable
     * @return rDestVar: The variable in the destiny modelpart
     * @return MaxIter: Maximum number of iteration allowed
     * @return TolIter: Tolerance accepted in the iteration
     * @return sign_pos: Positive or negative projection
     */

    void ScalarMap(
            const Variable<double> & rOriginVar,
            Variable<double> & rDestVar,
            const int MaxIter,
            const double TolIter,
            const bool sign_pos
            );

    /**
     * It maps a variable (vector) from a model part to other
     * @param rOriginVar: The original value of the variable
     * @return rDestVar: The variable in the destiny modelpart
     * @return MaxIter: Maximum number of iteration allowed
     * @return TolIter: Tolerance accepted in the iteration
     * @return sign_pos: Positive or negative projection
     */

    void VectorMap(
            const Variable< array_1d<double,3> > & rOriginVar,
            Variable< array_1d<double,3> > & rDestVar,
            const int MaxIter,
            const double TolIter,
            const bool sign_pos,
            const bool distributed
            );

    /**
     * It searches neighbours nodes in a specific radius
     * @param SearchRadiusFactor: The radius of search
     * @return A value of a close Gauss point, or alternative
     */

    void FindNeighbours(double SearchRadiusFactor);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    const ModelPart& mrOriginModelPart; // The original model part
    ModelPart& mrDestinationModelPart;  // The destination model part
    unsigned int mBucketSize;           // Bucket size for kd-tree
    GaussPointVector mGaussPointList;   // The list of Gauss points

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * It calculates the normal and area of a condition
     * @param pCond: The pointer to the condition
     * @return Normal: The normal of the condition
     * @return Area: The area of the condition
     * @param dimension: 2D/3D case
     */

    void CalcNormalAndArea(
            const Condition::Pointer pCond,
            array_1d<double,3>& Normal,
            double& Area,
            const int dimension
            );

    /**
     * It calculates the the center and raidus of a line
     * @param pCond: The pointer to the condition
     * @return Center: Center point (3D)
     * @return Radius: The radius of the line (half lenght)
     */

    void LineCenterAndRadius(
            const Condition::Pointer pCond,
            Point<3>& Center,
            double& Radius
            );

    /**
     * It calculates the the center and radius of a triangle
     * @param pCond: The pointer to the condition
     * @return Center: Center point (3D)
     * @return Radius: The radius of the line (half lenght)
     */

    void TriangleCenterAndRadius(
            const Condition::Pointer pCond,
            Point<3>& Center,
            double& Radius
            );

    /**
     * Desired outcome: It sets the projectioon of a Gauss node to a condition
     * @param GaussPoint: The origin Gauss Point
     * @return pCandidateCond: The candidate condition
     * @param Dist: The distance between the node and the Gauss Point
     * @param dimension: 2D/3D case
     */

    void SetProjectionToCond(
            GaussPointItem& GaussPoint,
            Condition::Pointer pCandidateCond,
            const int dimension
            );

    /**
     * Alternative when no condition is available: It sets the projection of a Gauss point to a node
     * @param GaussPoint: The origin Gauss Point
     * @return pCandidateNode: The candidate node
     * @param Dist: The distance between the node and the Gauss Point
     */

    void SetProjectionToNode(
            GaussPointItem& GaussPoint,
            Node<3>::Pointer pCandidateNode,
            const double& Dist
            );

    /**
     *  Test function, stores the distance between a Gauss Point and its projection
     */

    void DistanceCheck();

    /**
     *  Auxiliar function to compute the nodal length/area of each node in both origin and destiny model parts.
     */

    void ComputeNodalLengthArea();
    
    /**
     *  Auxiliar function to compute the equivalent nodal tractions to point loads minimizing the L2 norm of the error.
     */
    
    void ComputeEquivalentTractions(
            const Variable<array_1d<double,3> >& rOriginVar,
            const int MaxIter,
            const double TolIter
            );
    
    /**
     *  Auxiliar function to compute the equivalent nodal tractions to point loads minimizing the L2 norm of the error.
     */
    
    void ComputeNodalLoadsFromTractions(
            const Variable<array_1d<double,3> >& rDestVar
            );

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Un accessible methods
    ///@{
    ///@}


};  // Class AdvancedNMPointsMapper

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

}  // namespace Kratos.

#endif // KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED  defined
