// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED )
#define  KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include <map>
#include <set>
#include <math.h> 
#include "contact_structural_mechanics_application_variables.h"

// The geometry of the triangle for the "tessellation"
/* TRIANGLES */
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
/* QUADRILATERALS */
#include "geometries/quadrilateral_2d_4.h"

// /* The integration points (we clip triangles in 3D, so with line and triangle is enought)*/
// #include "integration/line_gauss_legendre_integration_points.h"
#include "integration/triangle_gauss_legendre_integration_points.h"

/* Utilities */
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/qr_utility.h"            //QR decomposition utility used in matrix inversion.

namespace Kratos
{
///@name Kratos Globals
///@{
    
///@}
///@name Type Definitions
///@{
    
    typedef Point<3>                                             PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                            GeometryNodeType;
    typedef Geometry<PointType>                          GeometryPointType;
    
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod              IntegrationMethod;
    typedef IntegrationPoint<2>                       IntegrationPointType;
    typedef GeometryNodeType::IntegrationPointsArrayType IntegrationPointsType;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

/** \brief ExactMortarIntegrationUtility 
 * This utility calculates the exact integration necessary for the Mortar Conditions
 */
template< unsigned int TDim, unsigned int TNumNodes>
class ExactMortarIntegrationUtility
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ExactMortarIntegrationUtility);
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Constructor
    
    /**
     * @param SlaveGeometry: The geometry of the slave condition
     * @param SlaveNormal: The normal of the slave condition
     * @param IntegrationOrder: The integration order to consider
     */
    
    ExactMortarIntegrationUtility(
        const unsigned int IntegrationOrder = 0
        )
    :mIntegrationOrder(IntegrationOrder)
    {
        GetIntegrationMethod();
    }
    
    /// Destructor.
    virtual ~ExactMortarIntegrationUtility() {}
    
    ///@}
    ///@name Operators
    ///@{
    
    ///@}
    ///@name Operations
    ///@{    
    
    /**
     * Get the integration method to consider
     */
        
    void GetIntegrationMethod()
    {
        // Setting the auxiliar integration points
        if (mIntegrationOrder == 1)
        {
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if (mIntegrationOrder == 2)
        {
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if (mIntegrationOrder == 3)
        {
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if (mIntegrationOrder == 4)
        {
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if (mIntegrationOrder == 5)
        {
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
        {
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
    }
    
    /**
     * Get the integration method to consider
     */
        
    GeometryNodeType::IntegrationPointsArrayType GetIntegrationTriangle()
    {
        // Setting the auxiliar integration points
        if (mIntegrationOrder == 1)
        {
            return Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        }
        else if (mIntegrationOrder == 2)
        {
            return Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        }
        else if (mIntegrationOrder == 3)
        {
            return Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        }
        else if (mIntegrationOrder == 4)
        {
            return Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        }
        else if (mIntegrationOrder == 5)
        {
            return Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        }
        else
        {
            return Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        }
    }
    
    /**
     * This utility computes the exact integration of the mortar condition
     * @param MasterGeometry: The geometry of the master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    inline bool GetExactIntegration(    
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        );
    
    /**
     * This utility computes the exact integration of the mortar condition
     * @param SlaveCond: The slave condition
     * @param MasterCond: The master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    bool TestGetExactIntegration(     
        Condition::Pointer& SlaveCond,
        Condition::Pointer& MasterCond,
        Matrix& CustomSolution
        )
    {
        IntegrationPointsType IntegrationPointsSlave;
        
        const bool solution = GetExactIntegration(SlaveCond->GetGeometry(), SlaveCond->GetValue(NORMAL), MasterCond->GetGeometry(), IntegrationPointsSlave);
        
        CustomSolution.resize(IntegrationPointsSlave.size(), TDim, false);
        
//         std::cout << "The Gauss Points obtained are: " << std::endl;
        for (unsigned int GP = 0; GP < IntegrationPointsSlave.size(); GP++)
        {
//             // For debug
//             KRATOS_WATCH(IntegrationPointsSlave[GP]);
            
            // Solution save:
            CustomSolution(GP, 0) = IntegrationPointsSlave[GP].Coordinate(1);
            if (TDim == 2)
            {
                CustomSolution(GP, 1) = IntegrationPointsSlave[GP].Weight();
            }
            else
            {
                CustomSolution(GP, 1) = IntegrationPointsSlave[GP].Coordinate(2);
                CustomSolution(GP, 2) = IntegrationPointsSlave[GP].Weight();
            }
        }
        
        return solution;
    }
    
    /**
     * This function computes the moment fitting matrix to redistribute the Gauss Points in the new arbitrary geometry 
     * @param PointToRotate: The points from the origin geometry
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @param SlaveGeometry: The geometry to consider
     * @return IntegrationPointsSlave: The integrations points that belong to the slave
     */
    
//     void MomentFittingIntegrationPoints( 
//         IntegrationPointsType& IntegrationPointsSlave,
//         GeometryNodeType& SlaveGeometry
//         )
//     {
//         // Initial values
//         const unsigned int IntegrationPointsSlaveSize = IntegrationPointsSlave.size();
//         const unsigned int SlaveGeometrySize = SlaveGeometry.size();
//         
//         QR<double, row_major> QRDecomposition; // QR decomposition object
//         
//         // We calculate the Moment Fitting matrix
//         Matrix A (SlaveGeometrySize, IntegrationPointsSlaveSize);
//         Vector M = ZeroVector(SlaveGeometrySize);
//         
//         for (unsigned int i = 0; i < IntegrationPointsSlaveSize; i++)
//         {
//             Point<3> GP;
//             GP.Coordinates() = IntegrationPointsSlave[i].Coordinates();
//             
//             Vector N( SlaveGeometrySize );
//             SlaveGeometry.ShapeFunctionsValues( N, GP );
//             
//             for (unsigned int j = 0; j < SlaveGeometrySize; j++)
//             {          
//                 M[j] += N[j] * IntegrationPointsSlave[i].Weight();
//                 
//                 A(j, i) = N[j];
//             }
//         }
//         
//         // Now we apply the QR decomposition
//         Vector w( IntegrationPointsSlaveSize );
//         
//         QRDecomposition.compute(SlaveGeometrySize, IntegrationPointsSlaveSize, &A(0, 0));
//         QRDecomposition.solve( &M[0], &w[0] );
//         
//         // Finally the calculate the weights for each Gauss Point
//         for (unsigned int i = 0; i < IntegrationPointsSlaveSize; i++)
//         {            
//             IntegrationPointsSlave[i].Weight() = w[i];
//         }
//     }
    
    void MomentFittingIntegrationPoints( 
        IntegrationPointsType& IntegrationPointsSlave,
        GeometryNodeType& SlaveGeometry
        )
    {
        // Initial values
        const unsigned int IntegrationPointsSlaveSize = IntegrationPointsSlave.size();
        const unsigned int SlaveGeometrySize = SlaveGeometry.size();
        
        QR<double, row_major> QRDecomposition; // QR decomposition object
        
        // We calculate the Moment Fitting matrix
        Matrix A (IntegrationPointsSlaveSize, SlaveGeometrySize);
        Vector M (IntegrationPointsSlaveSize);
        
        for (unsigned int i = 0; i < IntegrationPointsSlaveSize; i++)
        {
            M[i] = IntegrationPointsSlave[i].Weight();
            
            Point<3> GP;
            GP.Coordinates() = IntegrationPointsSlave[i].Coordinates();
            
            Vector N( SlaveGeometrySize );
            SlaveGeometry.ShapeFunctionsValues( N, GP );
            
            for (unsigned int j = 0; j < SlaveGeometrySize; j++)
            {                
                A(i, j) = N[j];
            }
        }
        
        // Now we apply the QR decomposition
        Vector w( SlaveGeometrySize );
        
        QRDecomposition.compute(IntegrationPointsSlaveSize, SlaveGeometrySize, &A(0, 0));
        QRDecomposition.solve( &M[0], &w[0] );
        
        // Finally the calculate the weights for each Gauss Point
        for (unsigned int i = 0; i < IntegrationPointsSlaveSize; i++)
        {
            Point<3> GP;
            GP.Coordinates() = IntegrationPointsSlave[i].Coordinates();
            
            Vector N( SlaveGeometrySize );
            SlaveGeometry.ShapeFunctionsValues( N, GP );
            
            IntegrationPointsSlave[i].Weight() = inner_prod(N, w);
        }
    }
    
    /**
     * This function rotates to align the projected points to a parallel plane to XY
     * @param PointToRotate: The points from the origin geometry
     * @param PointReferenceRotation: The center point used as reference to rotate
     * @param SlaveNormal: The normal vector of the slave condition
     * @param SlaveTangentXi: The first tangent vector of the slave condition
     * @param SlaveTangentEta: The second tangent vector of the slave condition
     * @param Inversed: If we rotate to the XY or we recover from XY
     * @return PointRotated: The point rotated 
     */
    
    void RotatePoint( 
        Point<3>& PointToRotate,
        const Point<3> PointReferenceRotation,
        const array_1d<double, 3> SlaveTangentXi,
        const array_1d<double, 3> SlaveTangentEta,
        const bool Inversed
        )
    {                
        // We move to the (0,0,0)
        Point<3> AuxPointToRotate;
        AuxPointToRotate.Coordinates() = PointToRotate.Coordinates() - PointReferenceRotation.Coordinates();
        
        boost::numeric::ublas::bounded_matrix<double, 3, 3> RotationMatrix = ZeroMatrix(3, 3);
        
        if (Inversed == false)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                RotationMatrix(0, i) = SlaveTangentXi[i];
                RotationMatrix(1, i) = SlaveTangentEta[i];
            }
        }
        else
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                RotationMatrix(i, 0) = SlaveTangentXi[i];
                RotationMatrix(i, 1) = SlaveTangentEta[i];
            }
        }
        
        PointToRotate.Coordinates() = prod(RotationMatrix, AuxPointToRotate) + PointReferenceRotation.Coordinates();
    }
    
    /**
     * This function intersects two lines in a 2D plane
     * @param PointOrig: The points from the origin geometry
     * @param PointDest: The points in the destination geometry
     * @return PointIntersection: The intersection point if there is any
     * @return True if there is a intersection point, false otherwise
     */
    
    bool Clipping2D(
        Point<3>& PointIntersection, 
        const Point<3> PointOrig1,
        const Point<3> PointOrig2,
        const Point<3> PointDest1,
        const Point<3> PointDest2
        )
    {
        const double SOrig1Orig2X = PointOrig2.Coordinate(1) - PointOrig1.Coordinate(1);
        const double SOrig1Orig2Y = PointOrig2.Coordinate(2) - PointOrig1.Coordinate(2);
        const double SDest1Dest2X = PointDest2.Coordinate(1) - PointDest1.Coordinate(1);
        const double SDest1Dest2Y = PointDest2.Coordinate(2) - PointDest1.Coordinate(2);
        
        const double Denom = SOrig1Orig2X * SDest1Dest2Y - SDest1Dest2X * SOrig1Orig2Y;
    
        const double Tolerance = 1.0e-12;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        if (std::abs(Denom) < Tolerance) // NOTE: Collinear
        {
            return false;
        }
        
        const double SOrig1Dest1X = PointOrig1.Coordinate(1) - PointDest1.Coordinate(1);
        const double SOrig1Dest1Y = PointOrig1.Coordinate(2) - PointDest1.Coordinate(2);
        
        const double S = (SOrig1Orig2X * SOrig1Dest1Y - SOrig1Orig2Y * SOrig1Dest1X)/Denom;
        
        const double T = (SDest1Dest2X * SOrig1Dest1Y - SDest1Dest2Y * SOrig1Dest1X)/Denom;
        
        if (S >= -Tolerance && S <= (1.0 + Tolerance) && T >= -Tolerance && T <= (1.0 + Tolerance))
        {
            PointIntersection.Coordinate(1) = PointOrig1.Coordinate(1) + T * SOrig1Orig2X; 
            PointIntersection.Coordinate(2) = PointOrig1.Coordinate(2) + T * SOrig1Orig2Y; 
            
            return true;
        }
        else
        {
            return false;
        }
    }
    
    /**
     * This function calculates in 2D the normal vector to a given one
     * @param v: The vector to compute the normal 
     * @return n: The normal vector
     */
    
    array_1d<double, 3> GetNormalVector2D(const array_1d<double, 3> v)
    {
        array_1d<double, 3> n;

        n[0] = - v[1];
        n[1] =   v[0];
        n[2] =    0.0;
    
        return n;
    }
    
    /**
     * This function calculates in 2D the angle between two points
     * @param PointOrig1: The points from the origin geometry
     * @param PointOrig2: The points in the destination geometry
     * @param axis_1: The axis respect the angle is calculated
     * @param axis_2: The normal to the previous axis
     * @return angle: The angle formed
     */
    
    double AnglePoints(
        const Point<3> PointOrig1,
        const Point<3> PointOrig2,
        const array_1d<double, 3> axis_1,
        const array_1d<double, 3> axis_2
        )
    {
        array_1d<double, 3> local_edge = PointOrig2.Coordinates() - PointOrig1.Coordinates();
        if (norm_2(local_edge) > 0.0)
        {
            local_edge /= norm_2(local_edge);
        }
        
        const double xi  = inner_prod(axis_1, local_edge);
        const double eta = inner_prod(axis_2, local_edge);
        
//         // Debug
//         KRATOS_WATCH(local_edge);
//         KRATOS_WATCH(axis_1);
//         KRATOS_WATCH(axis_2);
//         std::cout << std::atan2(eta, xi) << " " << xi << " "<< eta << std::endl;
        
        return (std::atan2(eta, xi));
    }

    /**
     * This function checks in 2D if two nodes are the same one
     * @param PointOrig: The points from the origin geometry
     * @param PointDest: The points in the destination geometry
     * @return check: The check done
     */
    
    bool CheckPoints2D(
        const Point<3> PointOrig1,
        const Point<3> PointOrig2
        )
    {
        const double Tolerance = 1.0e-6;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        const bool x = (std::abs(PointOrig2.X() - PointOrig1.X()) < Tolerance) ? true : false;
        const bool y = (std::abs(PointOrig2.Y() - PointOrig1.Y()) < Tolerance) ? true : false;
        
        return (x&&y);
    }
    
    /**
     * This function checks in 2D if two nodes are the same one
     * @param PointOrig: The points from the origin geometry
     * @param PointDest: The points in the destination geometry
     * @return check: The check done
     */
    
    bool CheckPoints3D(
        const Point<3> PointOrig1,
        const Point<3> PointOrig2
        )
    {
        const double Tolerance = 1.0e-6;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        return (norm_2(PointOrig2.Coordinates() - PointOrig1.Coordinates()) < Tolerance) ? true : false;
    }
    
    /**
     * This functions recovers if the nodes are inside a Quadrilateral conformed with Point<3>, instead of Nodes<3> (The input ones)
     * @param SlaveGeometry: The quadrilateral of interest
     * @param rPoint: The point to study if is inside
     * @return True if the point is inside the geometry
     */
    
    bool FasIsInsideQuadrilateral2D(
        Quadrilateral2D4 <Point<3>>& SlaveGeometry,
        const GeometryNodeType::CoordinatesArrayType& rPoint
        )
    {
        const double Tolerance = 1.0e-6;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        GeometryNodeType::CoordinatesArrayType rResult;
        
        Matrix J;

        rResult.clear();
        
        Vector DeltaXi = ZeroVector( 3 );
        
        GeometryNodeType::CoordinatesArrayType CurrentGlobalCoords( ZeroVector( 3 ) );

        //Newton iteration:

        const unsigned int maxiter = 30;

        for ( unsigned int k = 0; k < maxiter; k++ )
        {
            CurrentGlobalCoords = ZeroVector( 3 );
            SlaveGeometry.GlobalCoordinates( CurrentGlobalCoords, rResult );
            noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;
            SlaveGeometry.InverseOfJacobian( J, rResult );
            Matrix Jaux = ZeroMatrix(3, 3);
            subrange(Jaux, 0, 2, 0, 2) = J;
            noalias( DeltaXi ) = prod( Jaux, CurrentGlobalCoords );
            noalias( rResult ) += DeltaXi;

            if ( norm_2( DeltaXi ) > 30 )
            {
                break;
            }

            if ( norm_2( DeltaXi ) <  Tolerance )
            {
                break;
            }
        }
        
        if ( std::abs(rResult[0]) <= (1.0 + Tolerance) )
        {
            if ( std::abs(rResult[1]) <= (1.0 + Tolerance) )
            {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     * This functions calculates the determinant of a 2D triangle (using points) to check if invert the order
     * @param PointOrig1: First point
     * @param PointOrig2: Second point
     * @param PointOrig3: Third point
     * @return The DetJ
     */
    
    double FasTriagleCheck2D(
        const Point<3> PointOrig1,
        const Point<3> PointOrig2,
        const Point<3> PointOrig3
        )
    {
        const double x10 = PointOrig2.X() - PointOrig1.X();
        const double y10 = PointOrig2.Y() - PointOrig1.Y();

        const double x20 = PointOrig3.X() - PointOrig1.X();
        const double y20 = PointOrig3.Y() - PointOrig1.Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|	        |=	|	      |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|
        
        return x10 * y20 - y10 * x20;
    }
    
    /**
     * This functions calculates the determinant of a 3D condition (local dimension 2)
     * @param geom: The geometry where calculate the jacobian
     * @param GP: The gauss point
     * @return The DetJ
     */
    
    double DetJNonSquare(
        GeometryNodeType& geom,
        Point<3> GP
        )
    {
        Matrix J;
        geom.Jacobian(J, GP);
        const Matrix JTJ = prod( trans(J), J );
        
        return std::sqrt( JTJ(0,0) * JTJ(1,1) - JTJ(1,0) * JTJ(0,1) );
    }
    
    double DetJNonSquare(
        GeometryPointType& geom,
        Point<3> GP
        )
    {
        Matrix J;
        geom.Jacobian(J, GP);
        const Matrix JTJ = prod( trans(J), J );
        
        return std::sqrt( JTJ(0,0) * JTJ(1,1) - JTJ(1,0) * JTJ(0,1) );
    }
    
    /**
     * This function gives you the indexes needed to order a vector 
     * @param vect: The vector to order
     * @return idx: The vector of indexes
     */
    
    template <typename TType>
    std::vector<size_t> SortIndexes(const std::vector<TType> &vect) 
    {
        // Initialize original index locations
        std::vector<size_t> idx(vect.size());
        iota(idx.begin(), idx.end(), 0);

        // Sort indexes based on comparing values in vect
        std::sort(idx.begin(), idx.end(),
            [&vect](size_t i1, size_t i2) {return vect[i1] < vect[i2];});

        return idx;
    }
    
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

    unsigned int mIntegrationOrder;          // The integration order to consider
    IntegrationMethod mAuxIntegrationMethod; // The auxiliar list of Gauss Points taken from the geometry
    
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
}; // Class ExactMortarIntegrationUtility

///@name Explicit Specializations
///@{

    template<>  
    inline bool ExactMortarIntegrationUtility<2,2>::GetExactIntegration(         
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {             
        // We take the geometry GP from the core 
        const double tol = 1.0e-4; 
        
        double total_weight = 0.0;
        array_1d<double,2> coor_aux = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType projected_gp_global;
        GeometryNodeType::CoordinatesArrayType projected_gp_local;
        double aux_dist = 0.0;
        
        // Declare the boolean of full integral
        bool full_int = true;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = SlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::ProjectDirection(MasterGeometry, SlaveGeometry[i_slave].Coordinates(), projected_gp_global, aux_dist, -normal ); // The opposite direction
            
            const bool inside = MasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local );
            
            if (inside == false)
            {
                full_int = false;
            }
            else
            {
                if (i_slave == 0)
                {
                    coor_aux[0] = - 1.0;
                }
                else if (i_slave == 1)
                {
                    coor_aux[1] =   1.0;
                }
            }
        }
        
        if (full_int == true)
        {
            total_weight = 2.0;
        }
        else
        {
            std::vector<double> aux_xi;
            for (unsigned int i_master = 0; i_master < 2; i_master++)
            {
                const bool inside = ContactUtilities::ProjectIterativeLine2D(SlaveGeometry, MasterGeometry[i_master].Coordinates(), projected_gp_global, SlaveNormal);
                
                if (inside == true)
                {
                    aux_xi.push_back(projected_gp_local[0]);
                }
            }
            
            if (aux_xi.size() == 1)
            {
                if (coor_aux[0] == - 1.0)
                {
                    coor_aux[1] = aux_xi[0];
                }
                else if (coor_aux[1] == 1.0)
                {
                    coor_aux[0] = aux_xi[0];
                }
                else
                {
                    std::cout << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!!" << std::endl;
                }
            }
            else  if (aux_xi.size() == 2)
            {
                if (aux_xi[0] < aux_xi[1])
                {
                    coor_aux[0] = aux_xi[0];
                    coor_aux[1] = aux_xi[1];
                }
                else
                {
                    coor_aux[1] = aux_xi[0];
                    coor_aux[0] = aux_xi[1];
                }
            }
            
            total_weight = coor_aux[1] - coor_aux[0];
        }
        
        if(total_weight < 0.0)
        {
            KRATOS_ERROR << "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates: "<< coor_aux << std::endl;
        }
        
        if (total_weight > tol)
//             if (total_weight > 0.0)
        {
            // With the proportion of the weigth you recalculate the integration weight. Change the coordinates of the integration to accomodate
            const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
            IntegrationPointsSlave.resize(IntegrationPoints.size());
            for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
            {
                const double weight = IntegrationPoints[PointNumber].Weight() * total_weight/2.0;
                const double xi = 0.5 * (1.0 - IntegrationPoints[PointNumber].Coordinate(1)) * coor_aux[0] 
                                + 0.5 * (1.0 + IntegrationPoints[PointNumber].Coordinate(1)) * coor_aux[1];
                
                IntegrationPointsSlave[PointNumber] = IntegrationPoint<2>( xi, weight );
            }

//             // Debug
//             std::cout <<  SlaveGeometry[0].X() << " " << SlaveGeometry[0].Y() << " " << SlaveGeometry[1].X() << " " << SlaveGeometry[1].Y() << std::endl;
//             std::cout <<  MasterGeometry[0].X() << " " << MasterGeometry[0].Y() << " " << MasterGeometry[1].X() << " " << MasterGeometry[1].Y() << std::endl;
//             KRATOS_WATCH(coor_aux);
// 
//             std::cout << "IntegrationPoints : " << IntegrationPointsSlave.size( ) << std::endl;
//             for ( unsigned int i_vec = 0; i_vec < IntegrationPointsSlave.size( ); ++i_vec )
//             {
//                 KRATOS_WATCH( IntegrationPointsSlave[i_vec] );
//             }
            
            return true;
        }
        else
        {
            IntegrationPointsSlave.clear();
//             IntegrationPointsSlave.resize(0, false);
            return false;
        }
    
        IntegrationPointsSlave.clear();
        return false;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    inline bool ExactMortarIntegrationUtility<3,3>::GetExactIntegration(    
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        
        // We define the tolerance
        const double Tolerance = std::numeric_limits<double>::epsilon();;
//         const double Tolerance = 1.0e-8;
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const Point<3> SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<Point<3>, 3> MasterProjectedPoint;
        array_1d<bool, 3> AllInside;
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            MasterProjectedPoint[i_node] = ContactUtilities::FastProject(SlaveCenter, MasterGeometry[i_node], SlaveNormal);
            
            GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
        
            AllInside[i_node] = SlaveGeometry.IsInside( MasterProjectedPoint[i_node].Coordinates( ), ProjectedGPLocal ) ;
        }
        
        // We create the pointlist
        std::vector<Point<3>> PointList;
        
        // No point inside
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false))
        {
            IntegrationPointsSlave.clear();
//             IntegrationPointsSlave.resize(0, false);
            return false;
        }
        // All the points inside
        else if ((AllInside[0] == true) &&
                 (AllInside[1] == true) &&
                 (AllInside[2] == true))
        {
            std::vector<Point<3>::Pointer> AllPointsArray (3);
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                AllPointsArray[i_node] = boost::make_shared<Point<3>>(MasterProjectedPoint[i_node]);
            }
            
            Triangle3D3 <Point<3>> AllTriangle( AllPointsArray );
        
            const double LocalArea = AllTriangle.Area();
            
            if (LocalArea > Tolerance) // NOTE: Just in case we are not getting a real area
            {
                // We initialize our auxiliar integration point vector
                const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
                const size_t LocalIntegrationSize = IntegrationPoints.size();
                
                IntegrationPointsSlave.resize(LocalIntegrationSize, false);
                
                // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                {                    
                    // We convert the local coordinates to global coordinates
                    Point<3> gp_local;
                    gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                    Point<3> gp_global;
                    AllTriangle.GlobalCoordinates(gp_global, gp_local);
                    
                    // We recover this point to the triangle plane
                    RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                    
                    // Now we are supposed to project to the slave surface, but like the point it is already in the slave surface (with a triangle we work in his plane) we just calculate the local coordinates
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global);
                    
                    // We can construct now the integration local triangle
                    const double DetJ1 = DetJNonSquare(SlaveGeometry, gp_local);
                    const double DetJ2 = DetJNonSquare(AllTriangle, gp_local);
                    IntegrationPointsSlave[PointNumber] = IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ2 / DetJ1);
                }
                
//                 // We correct the weights using the moment fitting 
//                 MomentFittingIntegrationPoints(IntegrationPointsSlave, SlaveGeometry);
                
                return true;
            }
            else
            {
                IntegrationPointsSlave.clear();
                return false;
            }
        }
        else
        {            
            // Before clipping we rotate to a XY plane
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                RotatePoint(SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                RotatePoint(MasterProjectedPoint[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                
                if (AllInside[i_node] == true)
                {
                    PointList.push_back(MasterProjectedPoint[i_node]);
                }
            }
        
            // We consider the Z coordinate constant
            const double ZRef = MasterProjectedPoint[0].Coordinate(3);
            
            // We find the intersection in each side
            std::map<unsigned int, unsigned int> MapEdges;
            for (unsigned int i_edge = 0; i_edge < 3; i_edge++)
            {
                MapEdges.insert(std::make_pair(i_edge, 0));
//                 MapEdges [i_edge] = 0;
                
                const unsigned int ip_edge = (i_edge == 2) ? 0 : i_edge + 1;
                for (unsigned int j_edge = 0; j_edge < 3; j_edge++)
                {
                    const unsigned int jp_edge = (j_edge == 2) ? 0 : j_edge + 1;
                    
                    Point<3> IntersectedPoint;
                    const bool intersected = Clipping2D(
                        IntersectedPoint,
                        SlaveGeometry[i_edge],
                        SlaveGeometry[ip_edge],
                        MasterProjectedPoint[j_edge],
                        MasterProjectedPoint[jp_edge]
                        );
                    
                    if (intersected == true)
                    {
                        bool AddPoint = true;
                        for (unsigned int iter = 0; iter < PointList.size(); iter++)
                        {
                            if (CheckPoints2D(IntersectedPoint, PointList[iter]) == true)
                            {
                                AddPoint = false;
                            }
                        }
                        
                        if (AddPoint == true) 
                        {
                            IntersectedPoint.Coordinate(3) = ZRef;
                            PointList.push_back(IntersectedPoint);
                            MapEdges[i_edge] += 1;
                        }
                    }
                }
            }
            
            // No we check with edges are split just one time (which means that the corner belongs to the intersection)
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                unsigned int il_node = (i_node == 0) ? 2 : i_node - 1; // The first node is in edge 1 and 3
                
                if ((MapEdges[i_node]  == 1) && (MapEdges[il_node] == 1))
                {
                    bool AddPoint = true;
                    for (unsigned int iter = 0; iter < PointList.size(); iter++)
                    {
                        if (CheckPoints2D(SlaveGeometry[i_node], PointList[iter]) == true)
                        {
                            AddPoint = false;
                        }
                    }
                    
                    if (AddPoint == true) 
                    {
                        PointList.push_back(SlaveGeometry[i_node]);
                    }
                }
            }
            
            // We compose the triangles 
            const unsigned int ListSize = PointList.size();
            if (ListSize > 2) // Technically the minimum is three, just in case I consider 2
            {
                // We reorder the nodes according with the angle they form with the first node
                std::vector<double> Angles (ListSize - 1);
                array_1d<double, 3> v = PointList[1].Coordinates() - PointList[0].Coordinates();
                v /= norm_2(v);
                array_1d<double, 3> n = GetNormalVector2D(v);
                
                for (unsigned int elem = 1; elem < ListSize; elem++)
                {
                    Angles[elem - 1] = AnglePoints(PointList[0], PointList[elem], v, n);
                    if (Angles[elem - 1] < 0.0)
                    {
                        v = MasterProjectedPoint[elem].Coordinates() - MasterProjectedPoint[0].Coordinates();
                        v /= norm_2(v);
                        n = GetNormalVector2D(v);
                        for (unsigned int auxelem = 0; auxelem <= (elem - 1); auxelem++)
                        {
                            Angles[auxelem] -= Angles[elem - 1];
                        }
                    }
                }
                
                const std::vector<size_t> IndexVector = SortIndexes<double>(Angles);
                
                std::vector<Point<3>::Pointer> PointsArray (3);
                
                // We initialize our auxiliar integration point vector
                const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
                const size_t LocalIntegrationSize = IntegrationPoints.size();
                
//                 IntegrationPointsSlave.resize((ListSize - 2) * LocalIntegrationSize, false);
                IntegrationPointsSlave.clear();
            
                for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
                {
                    // NOTE: We add 1 because we removed from the list the fisrt point
                    if (FasTriagleCheck2D(PointList[0], PointList[IndexVector[elem] + 1], PointList[IndexVector[elem + 1] + 1]) > 0.0)
                    {
                        PointsArray[0] = boost::make_shared<Point<3>>(PointList[0]);
                        PointsArray[1] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 1] + 1]);
                    }
                    else
                    {
                        PointsArray[0] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 1] + 1]);
                        PointsArray[1] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<Point<3>>(PointList[0]);
                    }
                    
                    // We create the triangle
                    Triangle2D3 <Point<3>> triangle( PointsArray );
                    
                    // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                    const double LocalArea = triangle.Area();

                    if (LocalArea > Tolerance) // NOTE: Just in case we are not getting a real area
                    {
                        // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                        for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                        {                    
                            // We convert the local coordinates to global coordinates
                            Point<3> gp_local;
                            gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                            Point<3> gp_global;
                            triangle.GlobalCoordinates(gp_global, gp_local);
                            
                            // We recover this point to the triangle plane
                            RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                            
                            // Now we are supposed to project to the slave surface, but like the point it is already in the slave surface (with a triangle we work in his plane) we just calculate the local coordinates
                            SlaveGeometry.PointLocalCoordinates(gp_local, gp_global);
                            
                            // We can construct now the integration local triangle
                            const double DetJ1 = DetJNonSquare(SlaveGeometry, gp_local);
                            const double DetJ2 = triangle.DeterminantOfJacobian(gp_local);
                            IntegrationPointsSlave.push_back( IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ2 / DetJ1 ));
//                             IntegrationPointsSlave[elem * LocalIntegrationSize + PointNumber] = IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ );
                        }
                    }
                }
                
                if (IntegrationPointsSlave.size() > 0)
                {
//                     // We correct the weights using the moment fitting
//                     MomentFittingIntegrationPoints(IntegrationPointsSlave, SlaveGeometry);
                    
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else // No intersection
            {
                IntegrationPointsSlave.clear();
//                 IntegrationPointsSlave.resize(0, false);
                return false;
            }
        }
        
        IntegrationPointsSlave.clear();
        return false;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    inline bool ExactMortarIntegrationUtility<3,4>::GetExactIntegration(   
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {        
        // We define the tolerance
        const double Tolerance = 1.0e-8;
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const Point<3> SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<Point<3>, 4> SlaveProjectedPoint;
        array_1d<Point<3>, 4> MasterProjectedPoint;
        array_1d<bool, 4> AllInside;
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            SlaveProjectedPoint[i_node]  = ContactUtilities::FastProject( SlaveCenter,  SlaveGeometry[i_node], SlaveNormal);
            MasterProjectedPoint[i_node] = ContactUtilities::FastProject( SlaveCenter, MasterGeometry[i_node], SlaveNormal);
        }
        
        // Before clipping we rotate to a XY plane
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            RotatePoint( SlaveProjectedPoint[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
            RotatePoint(MasterProjectedPoint[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
        }
        
        std::vector<Point<3>::Pointer> DummyPointsArray (4);
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            DummyPointsArray[i_node] = boost::make_shared<Point<3>>(SlaveProjectedPoint[i_node]);
        }
        
        Quadrilateral2D4 <Point<3>> DummyQuadrilateral( DummyPointsArray );
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            GeometryNodeType::CoordinatesArrayType rResult;
            AllInside[i_node] = FasIsInsideQuadrilateral2D( DummyQuadrilateral, MasterProjectedPoint[i_node].Coordinates( ) ) ;
        }
        
        // We create the pointlist
        std::vector<Point<3>> PointList;
        
        // No point inside
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false) &&
            (AllInside[3] == false))
        {
            IntegrationPointsSlave.clear();
//             IntegrationPointsSlave.resize(0, false);
            return false;
        }
        // All the points inside
        else if ((AllInside[0] == true) &&
                 (AllInside[1] == true) &&
                 (AllInside[2] == true) &&
                 (AllInside[3] == true))
        {            
            // We reorder the nodes according with the angle they form with the first node
            std::vector<double> Angles (3);
            array_1d<double, 3> v = MasterProjectedPoint[1].Coordinates() - MasterProjectedPoint[0].Coordinates();
            v /= norm_2(v);
            array_1d<double, 3> n = GetNormalVector2D(v);
            
            for (unsigned int elem = 1; elem < 4; elem++)
            {
                Angles[elem - 1] = AnglePoints(MasterProjectedPoint[0], MasterProjectedPoint[elem], v, n);
                if (Angles[elem - 1] < 0.0)
                {
                    v = MasterProjectedPoint[elem].Coordinates() - MasterProjectedPoint[0].Coordinates();
                    v /= norm_2(v);
                    n = GetNormalVector2D(v);
                    for (unsigned int auxelem = 0; auxelem <= (elem - 1); auxelem++)
                    {
                        Angles[auxelem] -= Angles[elem - 1];
                    }
                }
            }
            
            const std::vector<size_t> IndexVector = SortIndexes<double>(Angles);
            
            std::vector<Point<3>::Pointer> PointsArray (3);
//             
            // We initialize our auxiliar integration point vector
            const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = GetIntegrationTriangle();
            const size_t LocalIntegrationSize = IntegrationPoints.size();
            
//             IntegrationPointsSlave.resize(2 * LocalIntegrationSize, false);
            IntegrationPointsSlave.clear();
            for (unsigned int elem = 0; elem < 2; elem++) // NOTE: We always have two points less that the number of nodes
            {
                // NOTE: We add 1 because we removed from the list the fisrt point
                if (FasTriagleCheck2D(MasterProjectedPoint[0], MasterProjectedPoint[IndexVector[elem] + 1], MasterProjectedPoint[IndexVector[elem + 1] + 1]) > 0.0)
                {
                    PointsArray[0] = boost::make_shared<Point<3>>(MasterProjectedPoint[0]);
                    PointsArray[1] = boost::make_shared<Point<3>>(MasterProjectedPoint[IndexVector[elem + 0] + 1]); 
                    PointsArray[2] = boost::make_shared<Point<3>>(MasterProjectedPoint[IndexVector[elem + 1] + 1]);
                }
                else
                {
                    PointsArray[0] = boost::make_shared<Point<3>>(MasterProjectedPoint[IndexVector[elem + 1] + 1]);
                    PointsArray[1] = boost::make_shared<Point<3>>(MasterProjectedPoint[IndexVector[elem + 0] + 1]); 
                    PointsArray[2] = boost::make_shared<Point<3>>(MasterProjectedPoint[0]);
                }
                
                // We create the triangle
                Triangle2D3 <Point<3>> triangle( PointsArray );
                
                // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                const double LocalArea = triangle.Area();
                
                if (LocalArea > Tolerance)
                {                       
//                     // Debug
//                     Point<3> aux1;
//                     aux1.Coordinates() = MasterProjectedPoint[0].Coordinates();
//                     RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     Point<3> aux2;
//                     aux2.Coordinates() = MasterProjectedPoint[IndexVector[elem + 0] + 1].Coordinates();
//                     RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     Point<3> aux3;
//                     aux3.Coordinates() = MasterProjectedPoint[IndexVector[elem + 1] + 1].Coordinates();
//                     RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                    
//                     std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                    
                    // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                    for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                    {                    
                        // We convert the local coordinates to global coordinates
                        Point<3> gp_local;
                        gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                        Point<3> gp_global;
                        triangle.GlobalCoordinates(gp_global, gp_local);
                        
                        // We recover this point to the triangle plane
                        RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                        
                        // Now we project to the slave surface
                        Point<3> gp_global_proj = ContactUtilities::FastProject(SlaveCenter, gp_global, - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        
                        // We can construct now the integration local triangle
                        const double DetJ1 = DetJNonSquare(SlaveGeometry, gp_local);
                        const double DetJ2 = triangle.DeterminantOfJacobian(gp_local); 
                        IntegrationPointsSlave.push_back(IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ2 / DetJ1 )); 
//                         IntegrationPointsSlave.push_back(IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), 4.0 * IntegrationPoints[PointNumber].Weight() * DetJ )); // NOTE:  The total weigh of a quadrilateral is 4
//                         IntegrationPointsSlave[elem * LocalIntegrationSize + PointNumber] = IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ ); // NOTE:  The total weigh of a quadrilateral is 4
                    }
                }
//                 // Debug
//                 else
//                 {
//                     std::cout << std::endl;
//                     KRATOS_WATCH(LocalArea);
//                     Point<3> aux1;
//                     aux1.Coordinates() = MasterProjectedPoint[0].Coordinates();
//                     RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     KRATOS_WATCH(aux1);
//                     
//                     Point<3> aux2;
//                     aux2.Coordinates() = MasterProjectedPoint[IndexVector[elem + 0] + 1].Coordinates();
//                     RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     KRATOS_WATCH(aux2);
//                     
//                     Point<3> aux3;
//                     aux3.Coordinates() = MasterProjectedPoint[IndexVector[elem + 1] + 1].Coordinates();
//                     RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     KRATOS_WATCH(aux3);
//                     
//                     KRATOS_ERROR;
//                 }
            }
            
            if (IntegrationPointsSlave.size() > 0)
            {
//                 // We correct the integration weights using the moment fitting
//                 MomentFittingIntegrationPoints(IntegrationPointsSlave, SlaveGeometry);
                
//                 // Debug
//                 for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++)
//                 {
//                     Point<3> gp_local;
//                     gp_local.Coordinates() = IntegrationPointsSlave[PointNumber].Coordinates();
//                     Point<3> gp_global;
//                     SlaveGeometry.GlobalCoordinates(gp_global, gp_local);
//                                                 
// //                         std::cout << IntegrationPoints[PointNumber].Weight() * DetJ2 / DetJ1 << "+";// << std::endl;
//                     std::cout << "Graphics3D[Sphere[{" << gp_global.X() << "," << gp_global.Y() << "," << gp_global.Z()  << "}," << 0.05 * IntegrationPointsSlave[PointNumber].Weight()<< "]],";// << std::endl;
//                 }
                
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            // We add the internal nodes
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                if (AllInside[i_node] == true)
                {
                    PointList.push_back(MasterProjectedPoint[i_node]);
                }
            }
            
            // We consider the Z coordinate constant
            const double ZRef = MasterProjectedPoint[0].Coordinate(3);
            
            // We find the intersection in each side
            std::map<unsigned int, unsigned int> MapEdges;
            for (unsigned int i_edge = 0; i_edge < 4; i_edge++)
            {
                MapEdges.insert(std::make_pair(i_edge, 0));
    //             MapEdges [i_edge] = 0;
                
                const unsigned int ip_edge = (i_edge == 3) ? 0 : i_edge + 1;
                for (unsigned int j_edge = 0; j_edge < 4; j_edge++)
                {
                    const unsigned int jp_edge = (j_edge == 3) ? 0 : j_edge + 1;
                    
                    Point<3> IntersectedPoint;
                    const bool intersected = Clipping2D(
                        IntersectedPoint,
                        SlaveProjectedPoint[i_edge],
                        SlaveProjectedPoint[ip_edge],
                        MasterProjectedPoint[j_edge],
                        MasterProjectedPoint[jp_edge]
                        );
                    
//                     // Debug
//                     std::cout << std::endl;
//                     Point<3> aux;
//                     aux.Coordinates() = SlaveProjectedPoint[i_edge].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     KRATOS_WATCH(aux);
//                     aux.Coordinates() = SlaveProjectedPoint[ip_edge].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     KRATOS_WATCH(aux);
//                     aux.Coordinates() = MasterProjectedPoint[j_edge].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     KRATOS_WATCH(aux);
//                     aux.Coordinates() = MasterProjectedPoint[jp_edge].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     KRATOS_WATCH(aux);
//                     std::cout << std::endl;
                    
                    if (intersected == true)
                    {
//                         // Debug
//                         std::cout << "Intersected" << std::endl;
//                         std::cout << std::endl;
//                         Point<3> aux;
//                         aux.Coordinates() = IntersectedPoint.Coordinates();
//                         RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         KRATOS_WATCH(aux);
//                         aux.Coordinates() = SlaveProjectedPoint[i_edge].Coordinates();
//                         RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         KRATOS_WATCH(aux);
//                         aux.Coordinates() = SlaveProjectedPoint[ip_edge].Coordinates();
//                         RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         KRATOS_WATCH(aux);
//                         aux.Coordinates() = MasterProjectedPoint[j_edge].Coordinates();
//                         RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         KRATOS_WATCH(aux);
//                         aux.Coordinates() = MasterProjectedPoint[jp_edge].Coordinates();
//                         RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         KRATOS_WATCH(aux);
                        
                        bool AddPoint = true;
                        for (unsigned int iter = 0; iter < PointList.size(); iter++)
                        {
                            if (CheckPoints2D(IntersectedPoint, PointList[iter]) == true)
                            {
                                AddPoint = false;
                            }
                        }
                        
                        if (AddPoint == true) 
                        {                            
                            IntersectedPoint.Coordinate(3) = ZRef;
                            PointList.push_back(IntersectedPoint);
                            MapEdges[i_edge] += 1;
                        }
                    }
                }
            }
            
//             // Debug
//             KRATOS_ERROR;
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                unsigned int il_node = (i_node == 0) ? 3 : i_node - 1; // The first node is in edge 1 and 4
                
                if ((MapEdges[i_node] == 1) && (MapEdges[il_node] == 1))
                {
                    bool AddPoint = true;
                    for (unsigned int iter = 0; iter < PointList.size(); iter++)
                    {
                        if (CheckPoints2D(SlaveProjectedPoint[i_node], PointList[iter]) == true) 
                        {
                            AddPoint = false;
                        }
                    }
                    
                    if (AddPoint == true) 
                    {
                        PointList.push_back(SlaveProjectedPoint[i_node]);
                    }
                }
            }
            
//             // Debug 
// //             KRATOS_WATCH(PointList.size());
//             for (unsigned int i_list = 0; i_list < PointList.size(); i_list++)
//             {
// //                 KRATOS_WATCH(PointList[i_list]);
//                 Point<3> aux;
//                 aux.Coordinates() = PointList[i_list].Coordinates();
//                 RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                 std::cout << aux.X() << "\t" << aux.Y() << "\t" << aux.Z() << std::endl;
//             }
//             std::cout << std::endl;

            // We compose the triangles
            const unsigned int ListSize = PointList.size();
            
//             // Debug
//             if (ListSize > 0 && ListSize < 3)
//             {
//                 for (unsigned int i = 0; i < ListSize; i++)
//                 {
//                     Point<3> aux;
//                     aux.Coordinates() = PointList[i].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);            
//                     std::cout << aux.X() << "\t" << aux.Y() << "\t" << aux.Z()  << std::endl;
//                 }
//                 std::cout << std::endl;
//                 for (unsigned int i_node = 0; i_node < 4; i_node++)
//                 {
//                     Point<3> aux;
//                     aux.Coordinates() = SlaveProjectedPoint[i_node].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);   
//                     
//                     std::cout << aux.X() << "\t" << aux.Y() << "\t" << aux.Z()  << std::endl;
//                 }
//                 std::cout << std::endl;
//                 for (unsigned int i_node = 0; i_node < 4; i_node++)
//                 {
//                     Point<3> aux;
//                     aux.Coordinates() = MasterProjectedPoint[i_node].Coordinates();
//                     RotatePoint(aux, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);   
//                     
//                     std::cout << aux.X() << "\t" << aux.Y() << "\t" << aux.Z()  << std::endl;
//                 }
//             }
            
            if (ListSize > 2) // Technically the minimum is three, just in case I consider 2
            {
                // We reorder the nodes according with the angle they form with the first node
                std::vector<double> Angles (ListSize - 1);
                array_1d<double, 3> v = PointList[1].Coordinates() - PointList[0].Coordinates();
                v /= norm_2(v);
                array_1d<double, 3> n = GetNormalVector2D(v);
                
                for (unsigned int elem = 1; elem < ListSize; elem++)
                {
                    Angles[elem - 1] = AnglePoints(PointList[0], PointList[elem], v, n);
                    if (Angles[elem - 1] < 0.0)
                    {
                        v = PointList[elem].Coordinates() - PointList[0].Coordinates();
                        v /= norm_2(v);
                        n = GetNormalVector2D(v);
                        for (unsigned int auxelem = 0; auxelem <= (elem - 1); auxelem++)
                        {
                            Angles[auxelem] -= Angles[elem - 1];
                        }
                    }
                }
                    
                const std::vector<size_t> IndexVector = SortIndexes<double>(Angles);
                
                std::vector<Point<3>::Pointer> PointsArray (3);
                
                // We initialize our auxiliar integration point vector
                const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = GetIntegrationTriangle();
                const size_t LocalIntegrationSize = IntegrationPoints.size();
                
//                 IntegrationPointsSlave.resize((ListSize - 2) * LocalIntegrationSize, false);
                IntegrationPointsSlave.clear();
        
                for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
                {
                    // NOTE: We add 1 because we removed from the list the fisrt point
                    if (FasTriagleCheck2D(PointList[0], PointList[IndexVector[elem] + 1], PointList[IndexVector[elem + 1] + 1]) > 0.0)
                    {
                        PointsArray[0] = boost::make_shared<Point<3>>(PointList[0]);
                        PointsArray[1] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 1] + 1]);
                    }
                    else
                    {
                        PointsArray[0] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 1] + 1]);
                        PointsArray[1] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<Point<3>>(PointList[0]);
                    }
                    
                    // We create the triangle
                    Triangle2D3 <Point<3>> triangle( PointsArray );
                    
                    // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                    const double LocalArea = triangle.Area(); // FIXME: Probably this is affecting to to the calculation of the local coordinates
                    
                    if (LocalArea > Tolerance)
                    {                                         
//                         // Debug
//                         Point<3> aux1;
//                         aux1.Coordinates() = PointList[0].Coordinates();
//                         RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         Point<3> aux2;
//                         aux2.Coordinates() = PointList[IndexVector[elem + 0] + 1].Coordinates();
//                         RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         Point<3> aux3;
//                         aux3.Coordinates() = PointList[IndexVector[elem + 1] + 1].Coordinates();
//                         RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                        
                        // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                        for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                        {                            
                            // We convert the local coordinates to global coordinates
                            Point<3> gp_local;
                            gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                            Point<3> gp_global;
                            triangle.GlobalCoordinates(gp_global, gp_local);
                            
                            // We recover this point to the triangle plane
                            RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                            
                            // Now we project to the slave surface
                            Point<3> gp_global_proj = ContactUtilities::FastProject(SlaveCenter, gp_global, - SlaveNormal); // We came back 
                            SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                            
                            // We can construct now the integration local triangle // FIXME: The weights I am getting are constant, probably you did something wrong
                            const double DetJ1 = DetJNonSquare(SlaveGeometry, gp_local);
                            const double DetJ2 = triangle.DeterminantOfJacobian(gp_local);
                            IntegrationPointsSlave.push_back(IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ2 / DetJ1 )); 
//                             IntegrationPointsSlave.push_back(IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), 4.0 * IntegrationPoints[PointNumber].Weight() * DetJ )); // NOTE: The total weight of a quadrilateral is 4
//                             IntegrationPointsSlave[elem * LocalIntegrationSize + PointNumber] = IntegrationPoint<3>( gp_local.Coordinate(1), gp_local.Coordinate(2), 4.0 * IntegrationPoints[PointNumber].Weight() * DetJ ); // NOTE: The total weight of a quadrilateral is 4
                        }
                    }
//                     // Debug
//                     else // FIXME: This is because the nodes are not in correct order, look how to do that (or zero area!!!)
//                     {
//                         std::cout << std::endl;
//                         KRATOS_WATCH(LocalArea);
//                         std::cout << CheckPoints2D(PointList[0], PointList[IndexVector[elem] + 1]) << "\t" << CheckPoints2D(PointList[IndexVector[elem + 1] + 1], PointList[IndexVector[elem] + 1]) << "\t" << CheckPoints2D(PointList[0], PointList[IndexVector[elem + 1] + 1]) << std::endl;
//                         Point<3> aux1;
//                         aux1.Coordinates() = PointList[0].Coordinates();
//                         RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         KRATOS_WATCH(aux1);
//                         
//                         Point<3> aux2;
//                         aux2.Coordinates() = PointList[IndexVector[elem + 0] + 1].Coordinates();
//                         RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         KRATOS_WATCH(aux2);
//                         
//                         Point<3> aux3;
//                         aux3.Coordinates() = PointList[IndexVector[elem + 1] + 1].Coordinates();
//                         RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         KRATOS_WATCH(aux3);
//                         
//                         KRATOS_ERROR;
//                     }
                }
                
                if (IntegrationPointsSlave.size() > 0)
                {
//                     // We correct the integration weights using the moment fitting
//                     MomentFittingIntegrationPoints(IntegrationPointsSlave, SlaveGeometry);
                    
//                     // Debug
//                     for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++)
//                     {
//                         Point<3> gp_local;
//                         gp_local.Coordinates() = IntegrationPointsSlave[PointNumber].Coordinates();
//                         Point<3> gp_global;
//                         SlaveGeometry.GlobalCoordinates(gp_global, gp_local);
//                                        
// //                         std::cout << IntegrationPoints[PointNumber].Weight() * DetJ2 / DetJ1 << "+";// << std::endl;
//                         std::cout << "Graphics3D[Sphere[{" << gp_global.X() << "," << gp_global.Y() << "," << gp_global.Z()  << "}," << 0.05 * IntegrationPointsSlave[PointNumber].Weight()<< "]],";// << std::endl;
//                     }
                    
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else // No intersection
            {
                IntegrationPointsSlave.clear();
//                 IntegrationPointsSlave.resize(0, false);
                return false;
            }
        }
        
        IntegrationPointsSlave.clear();
        return false;
    }
}

#endif  /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
