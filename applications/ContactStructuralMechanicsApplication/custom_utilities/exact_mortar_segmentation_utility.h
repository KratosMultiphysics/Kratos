// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
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

namespace Kratos
{
///@name Kratos Globals
///@{
    
///@}
///@name Type Definitions
///@{
    
    typedef Point<3>                                                 PointType;
    typedef Node<3>                                                   NodeType;
    typedef Geometry<NodeType>                                GeometryNodeType;
    typedef Geometry<PointType>                              GeometryPointType;
    
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod                  IntegrationMethod;
    typedef IntegrationPoint<2>                           IntegrationPointType;
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
    
    typedef std::vector<array_1d<PointType,TDim>>        ConditionArrayListType;
    
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
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave
        );
    
    /**
     * This utility computes the exact integration of the mortar condition (just the points, not the whole integration points)
     * @param MasterGeometry: The geometry of the master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    inline bool GetExactIntegration(    
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
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
        
        const bool solution = GetExactIntegration(SlaveCond->GetGeometry(), SlaveCond->GetValue(NORMAL), MasterCond->GetGeometry(), MasterCond->GetValue(NORMAL), IntegrationPointsSlave);
        
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
        PointType& PointToRotate,
        const PointType PointReferenceRotation,
        const array_1d<double, 3> SlaveTangentXi,
        const array_1d<double, 3> SlaveTangentEta,
        const bool Inversed
        )
    {                
        // We move to the (0,0,0)
        PointType AuxPointToRotate;
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
        PointType& PointIntersection, 
        const PointType PointOrig1,
        const PointType PointOrig2,
        const PointType PointDest1,
        const PointType PointDest2
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
        const PointType PointOrig1,
        const PointType PointOrig2,
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
        
        return (std::atan2(eta, xi));
    }

    /**
     * This function checks in 2D if two nodes are the same one
     * @param PointOrig: The points from the origin geometry
     * @param PointDest: The points in the destination geometry
     * @return check: The check done
     */
    
    bool CheckPoints2D(
        const PointType PointOrig1,
        const PointType PointOrig2
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
        const PointType PointOrig1,
        const PointType PointOrig2
        )
    {
        const double Tolerance = 1.0e-6;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        return (norm_2(PointOrig2.Coordinates() - PointOrig1.Coordinates()) < Tolerance) ? true : false;
    }
    
    /**
     * This functions recovers if the nodes are inside a Quadrilateral conformed with PointType, instead of Nodes<3> (The input ones)
     * @param SlaveGeometry: The quadrilateral of interest
     * @param rPoint: The point to study if is inside
     * @return True if the point is inside the geometry
     */
    
    bool FasIsInsideQuadrilateral2D(
        Quadrilateral2D4 <PointType>& SlaveGeometry,
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
        const PointType PointOrig1,
        const PointType PointOrig2,
        const PointType PointOrig3
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
     * This functions calculates the determinant of a 3D triangle (using points) to check if invert the order
     * @param PointOrig1: First point
     * @param PointOrig2: Second point
     * @param PointOrig3: Third point
     * @return The DetJ
     */
    
    double FasTriagleCheck3D(
        const PointType PointOrig1,
        const PointType PointOrig2,
        const PointType PointOrig3
        )
    {
        Matrix jacobian( 3, 2 );
        jacobian( 0, 0 ) = -( PointOrig1.X() ) + ( PointOrig2.X() ); //on the Gauss points (J is constant at each element)
        jacobian( 1, 0 ) = -( PointOrig1.Y() ) + ( PointOrig2.Y() );
        jacobian( 2, 0 ) = -( PointOrig1.Z() ) + ( PointOrig2.Z() );
        jacobian( 0, 1 ) = -( PointOrig1.X() ) + ( PointOrig3.X() );
        jacobian( 1, 1 ) = -( PointOrig1.Y() ) + ( PointOrig3.Y() );
        jacobian( 2, 1 ) = -( PointOrig1.Z() ) + ( PointOrig3.Z() );
        
        return MathUtils<double>::GeneralizedDet(jacobian);
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
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {             
        // We take the geometry GP from the core 
        const double Tolerance = 1.0e-4;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        double total_weight = 0.0;
        array_1d<double,2> coor_aux = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType projected_gp_global;
        GeometryNodeType::CoordinatesArrayType projected_gp_local;
        
        // Declare the boolean of full integral
        bool full_int = true;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = SlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::FastProjectDirection(MasterGeometry, SlaveGeometry[i_slave].Coordinates(), projected_gp_global, MasterNormal, -normal ); // The opposite direction
            
            const bool inside = MasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, Tolerance );
            
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
        
        if (total_weight > Tolerance)
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
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        
        // We define the tolerance
        const double Tolerance = 1.0e-8;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<PointType, 3> MasterProjectedPoint;
        array_1d<bool, 3> AllInside;
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            MasterProjectedPoint[i_node] = ContactUtilities::FastProject(SlaveCenter, MasterGeometry[i_node], SlaveNormal);
            
            GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
        
            AllInside[i_node] = SlaveGeometry.IsInside( MasterProjectedPoint[i_node].Coordinates( ), ProjectedGPLocal, Tolerance) ;
        }
        
        // We create the pointlist
        std::vector<PointType> PointList;
        
        // No point from the master is inside the slave
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false))
        {            
            // We check if all the nodes are inside the master element
            array_1d<PointType, 3> SlaveProjectedPoint;
            
            const PointType MasterCenter = MasterGeometry.Center();
            
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                SlaveProjectedPoint[i_node] = ContactUtilities::FastProject(MasterCenter, SlaveGeometry[i_node], MasterNormal);
                
                GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
            
                AllInside[i_node] = MasterGeometry.IsInside( SlaveProjectedPoint[i_node].Coordinates( ), ProjectedGPLocal, Tolerance) ;
            }
            
            // The whole slave is inside the master
            if ((AllInside[0] == true) &&
                (AllInside[1] == true) &&
                (AllInside[2] == true))
            {
                IntegrationPointsSlave = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
            }
            else
            {
                IntegrationPointsSlave.clear();
//                 IntegrationPointsSlave.resize(0, false);
                return false;
            }
        }
        // All the points inside
        else if ((AllInside[0] == true) &&
                 (AllInside[1] == true) &&
                 (AllInside[2] == true))
        {
            std::vector<PointType::Pointer> AllPointsArray (3);
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                AllPointsArray[i_node] = boost::make_shared<PointType>(MasterProjectedPoint[i_node]);
            }
            
            Triangle3D3 <PointType> AllTriangle( AllPointsArray );
        
            const double LocalArea = AllTriangle.Area();
            
            if (LocalArea > Tolerance) // NOTE: Just in case we are not getting a real area
            {
                // We initialize our auxiliar integration point vector
                const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
                const size_t LocalIntegrationSize = IntegrationPoints.size();
                
                IntegrationPointsSlave.resize(LocalIntegrationSize, false);
                
//                 // Debug
//                 PointType aux1;
//                 aux1.Coordinates() = MasterProjectedPoint[0].Coordinates();
//                 RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                 
//                 PointType aux2;
//                 aux2.Coordinates() = MasterProjectedPoint[1].Coordinates();
//                 RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                 
//                 PointType aux3;
//                 aux3.Coordinates() = MasterProjectedPoint[2].Coordinates();
//                 RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                 
//                 std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                
                // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                {                    
                    // We convert the local coordinates to global coordinates
                    PointType gp_local;
                    gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                    PointType gp_global;
                    AllTriangle.GlobalCoordinates(gp_global, gp_local);
                    
                    // Now we are supposed to project to the slave surface, but like the point it is already in the slave surface (with a triangle we work in his plane) we just calculate the local coordinates
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global);
                    
                    // We can construct now the integration local triangle
                    const double DetJ = AllTriangle.DeterminantOfJacobian(gp_local);
                    IntegrationPointsSlave[PointNumber] = IntegrationPointType( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ );
                }
                
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
                RotatePoint(       SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                RotatePoint(MasterProjectedPoint[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                
                if (AllInside[i_node] == true)
                {
                    PointList.push_back(MasterProjectedPoint[i_node]);
                }
            }
        
            // We consider the Z coordinate constant
            const double ZRef = SlaveCenter.Coordinate(3);
            
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
                    
                    PointType IntersectedPoint;
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
                
                std::vector<PointType::Pointer> PointsArray (3);
                
                // We initialize our auxiliar integration point vector
                const GeometryNodeType::IntegrationPointsArrayType& IntegrationPoints = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
                const size_t LocalIntegrationSize = IntegrationPoints.size();
            
                // We recover the rotation
                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    RotatePoint( SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                }
                
//                 IntegrationPointsSlave.resize((ListSize - 2) * LocalIntegrationSize, false);
                IntegrationPointsSlave.clear();
                
                for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
                {
                    // NOTE: We add 1 because we removed from the list the fisrt point
                    if (FasTriagleCheck2D(PointList[0], PointList[IndexVector[elem] + 1], PointList[IndexVector[elem + 1] + 1]) > 0.0)
                    {
                        PointsArray[0] = boost::make_shared<PointType>(PointList[0]);
                        PointsArray[1] = boost::make_shared<PointType>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<PointType>(PointList[IndexVector[elem + 1] + 1]);
                    }
                    else
                    {
                        PointsArray[0] = boost::make_shared<PointType>(PointList[IndexVector[elem + 1] + 1]);
                        PointsArray[1] = boost::make_shared<PointType>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<PointType>(PointList[0]);
                    }
                    
                    // We create the triangle
                    Triangle2D3 <PointType> triangle( PointsArray );
                    
                    // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                    const double LocalArea = triangle.Area();

                    if (LocalArea > Tolerance) // NOTE: Just in case we are not getting a real area
                    {
//                         // Debug
//                         PointType aux1;
//                         aux1.Coordinates() = MasterProjectedPoint[0].Coordinates();
//                         RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         PointType aux2;
//                         aux2.Coordinates() = MasterProjectedPoint[IndexVector[elem + 0] + 1].Coordinates();
//                         RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         PointType aux3;
//                         aux3.Coordinates() = MasterProjectedPoint[IndexVector[elem + 1] + 1].Coordinates();
//                         RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                            
                        // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                        for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                        {                    
                            // We convert the local coordinates to global coordinates
                            PointType gp_local;
                            gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                            PointType gp_global;
                            triangle.GlobalCoordinates(gp_global, gp_local);
                            
                            // We recover this point to the triangle plane
                            RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                            
                            // Now we are supposed to project to the slave surface, but like the point it is already in the slave surface (with a triangle we work in his plane) we just calculate the local coordinates
                            SlaveGeometry.PointLocalCoordinates(gp_local, gp_global);
                            
                            // We can construct now the integration local triangle
                            const double DetJ = triangle.DeterminantOfJacobian(gp_local);
                            IntegrationPointsSlave.push_back( IntegrationPointType( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ ));
//                             IntegrationPointsSlave[elem * LocalIntegrationSize + PointNumber] = IntegrationPointType( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ );
                        }
                    }
                }

                if (IntegrationPointsSlave.size() > 0)
                {                    
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else // No intersection
            {
                // We recover the rotation
                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    RotatePoint( SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                }
                
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
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {        
        // We define the tolerance
        const double Tolerance = 1.0e-8;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<PointType, 4> SlaveProjectedPoint;
        array_1d<PointType, 4> MasterProjectedPoint;
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
        
        std::vector<PointType::Pointer> DummyPointsArray (4);
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            DummyPointsArray[i_node] = boost::make_shared<PointType>(SlaveProjectedPoint[i_node]);
        }
        
        Quadrilateral2D4 <PointType> DummyQuadrilateral( DummyPointsArray );
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            GeometryNodeType::CoordinatesArrayType rResult;
            AllInside[i_node] = FasIsInsideQuadrilateral2D( DummyQuadrilateral, MasterProjectedPoint[i_node].Coordinates( ) ) ;
        }
        
        // We create the pointlist
        std::vector<PointType> PointList;
        
        // No point from the master is inside the slave
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false) &&
            (AllInside[3] == false))
        {
            // We check if all the nodes are inside the master element            
            const PointType MasterCenter = MasterGeometry.Center();
            
            const array_1d<double, 3> MasterTangentXi  = (MasterGeometry[2].Coordinates() - MasterGeometry[0].Coordinates())/norm_2(MasterGeometry[2].Coordinates() - MasterGeometry[0].Coordinates());
            const array_1d<double, 3> MasterTangentEta = MathUtils<double>::UnitCrossProduct(MasterTangentXi, MasterNormal);
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                SlaveProjectedPoint[i_node]  = ContactUtilities::FastProject( MasterCenter,  SlaveGeometry[i_node], MasterNormal);
                MasterProjectedPoint[i_node] = ContactUtilities::FastProject( MasterCenter, MasterGeometry[i_node], MasterNormal);
            }
            
            // Before clipping we rotate to a XY plane
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                RotatePoint( SlaveProjectedPoint[i_node], MasterCenter, MasterTangentXi, MasterTangentEta, false);
                RotatePoint(MasterProjectedPoint[i_node], MasterCenter, MasterTangentXi, MasterTangentEta, false);
            }
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                DummyPointsArray[i_node] = boost::make_shared<PointType>(MasterProjectedPoint[i_node]);
            }
            
            Quadrilateral2D4 <PointType> AuxDummyQuadrilateral( DummyPointsArray );
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                GeometryNodeType::CoordinatesArrayType rResult;
                AllInside[i_node] = FasIsInsideQuadrilateral2D( AuxDummyQuadrilateral, SlaveProjectedPoint[i_node].Coordinates( ) ) ;
            }
            
            // The whole slave is inside the master
            if ((AllInside[0] == true) &&
                (AllInside[1] == true) &&
                (AllInside[2] == true) &&
                (AllInside[3] == true))
            {
                IntegrationPointsSlave = SlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
            }
            else
            {
                IntegrationPointsSlave.clear();
    //             IntegrationPointsSlave.resize(0, false);
                return false;
            }
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
            
            std::vector<PointType::Pointer> PointsArray (3);
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
                    PointsArray[0] = boost::make_shared<PointType>(MasterProjectedPoint[0]);
                    PointsArray[1] = boost::make_shared<PointType>(MasterProjectedPoint[IndexVector[elem + 0] + 1]); 
                    PointsArray[2] = boost::make_shared<PointType>(MasterProjectedPoint[IndexVector[elem + 1] + 1]);
                }
                else
                {
                    PointsArray[0] = boost::make_shared<PointType>(MasterProjectedPoint[IndexVector[elem + 1] + 1]);
                    PointsArray[1] = boost::make_shared<PointType>(MasterProjectedPoint[IndexVector[elem + 0] + 1]); 
                    PointsArray[2] = boost::make_shared<PointType>(MasterProjectedPoint[0]);
                }
                
                // We create the triangle
                Triangle2D3 <PointType> triangle( PointsArray );
                
                // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                const double LocalArea = triangle.Area();
                
                if (LocalArea > Tolerance)
                {                       
//                     // Debug
//                     PointType aux1;
//                     aux1.Coordinates() = MasterProjectedPoint[0].Coordinates();
//                     RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     PointType aux2;
//                     aux2.Coordinates() = MasterProjectedPoint[IndexVector[elem + 0] + 1].Coordinates();
//                     RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     PointType aux3;
//                     aux3.Coordinates() = MasterProjectedPoint[IndexVector[elem + 1] + 1].Coordinates();
//                     RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                     
//                     std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                    
                    // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                    for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                    {                    
                        // We convert the local coordinates to global coordinates
                        PointType gp_local;
                        gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                        PointType gp_global;
                        triangle.GlobalCoordinates(gp_global, gp_local);
                        
                        // We recover this point to the triangle plane
                        RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                        
                        // Now we project to the slave surface
                        PointType gp_global_proj = ContactUtilities::FastProject(SlaveCenter, gp_global, - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        
                        // We can construct now the integration local triangle
                        const double DetJ = triangle.DeterminantOfJacobian(gp_local); 
                        IntegrationPointsSlave.push_back(IntegrationPointType( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ )); 
                    }
                }
            }
            
            if (IntegrationPointsSlave.size() > 0)
            {                
//                 // Debug
//                 for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++)
//                 {
//                     PointType gp_local;
//                     gp_local.Coordinates() = IntegrationPointsSlave[PointNumber].Coordinates();
//                     PointType gp_global;
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
            const double ZRef = SlaveCenter.Coordinate(3);
            
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
                    
                    PointType IntersectedPoint;
                    const bool intersected = Clipping2D(
                        IntersectedPoint,
                        SlaveProjectedPoint[i_edge],
                        SlaveProjectedPoint[ip_edge],
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
                
                std::vector<PointType::Pointer> PointsArray (3);
                
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
                        PointsArray[0] = boost::make_shared<PointType>(PointList[0]);
                        PointsArray[1] = boost::make_shared<PointType>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<PointType>(PointList[IndexVector[elem + 1] + 1]);
                    }
                    else
                    {
                        PointsArray[0] = boost::make_shared<PointType>(PointList[IndexVector[elem + 1] + 1]);
                        PointsArray[1] = boost::make_shared<PointType>(PointList[IndexVector[elem + 0] + 1]); 
                        PointsArray[2] = boost::make_shared<PointType>(PointList[0]);
                    }
                    
                    // We create the triangle
                    Triangle2D3 <PointType> triangle( PointsArray );
                    
                    // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                    const double LocalArea = triangle.Area(); // FIXME: Probably this is affecting to to the calculation of the local coordinates
                    
                    if (LocalArea > Tolerance)
                    {                                         
//                         // Debug
//                         PointType aux1;
//                         aux1.Coordinates() = PointList[0].Coordinates();
//                         RotatePoint(aux1, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         PointType aux2;
//                         aux2.Coordinates() = PointList[IndexVector[elem + 0] + 1].Coordinates();
//                         RotatePoint(aux2, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         PointType aux3;
//                         aux3.Coordinates() = PointList[IndexVector[elem + 1] + 1].Coordinates();
//                         RotatePoint(aux3, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
//                         
//                         std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
//                         
                        // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates                 
                        for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                        {                            
                            // We convert the local coordinates to global coordinates
                            PointType gp_local;
                            gp_local.Coordinates() = IntegrationPoints[PointNumber].Coordinates();
                            PointType gp_global;
                            triangle.GlobalCoordinates(gp_global, gp_local);
                            
                            // We recover this point to the triangle plane
                            RotatePoint(gp_global, SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                            
                            // Now we project to the slave surface
                            PointType gp_global_proj = ContactUtilities::FastProject(SlaveCenter, gp_global, - SlaveNormal); // We came back 
                            SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                            
                            // We can construct now the integration local triangle 
                            const double DetJ = triangle.DeterminantOfJacobian(gp_local);
                            IntegrationPointsSlave.push_back(IntegrationPointType( gp_local.Coordinate(1), gp_local.Coordinate(2), IntegrationPoints[PointNumber].Weight() * DetJ ));
                        }
                    }
                }
                
                if (IntegrationPointsSlave.size() > 0)
                {                    
//                     // Debug
//                     for (unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++)
//                     {
//                         PointType gp_local;
//                         gp_local.Coordinates() = IntegrationPointsSlave[PointNumber].Coordinates();
//                         PointType gp_global;
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

    /***********************************************************************************/
    /***********************************************************************************/


    template<>  
    inline bool ExactMortarIntegrationUtility<2,2>::GetExactIntegration(         
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
        )
    {             
        // We take the geometry GP from the core 
        const double Tolerance = 1.0e-4;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        double total_weight = 0.0;
        array_1d<double,2> coor_aux = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType projected_gp_global;
        GeometryNodeType::CoordinatesArrayType projected_gp_local;
        
        // Declare the boolean of full integral
        bool full_int = true;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = SlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::FastProjectDirection(MasterGeometry, SlaveGeometry[i_slave].Coordinates(), projected_gp_global, MasterNormal, -normal ); // The opposite direction
            
            const bool inside = MasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, Tolerance );
            
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
        
        if (total_weight > Tolerance)
        {
            ConditionsPointsSlave.resize(1);
            array_1d<PointType, 2> ListPoints;
            ListPoints[0].Coordinate(1) = coor_aux[0];
            ListPoints[1].Coordinate(1) = coor_aux[1];
            ConditionsPointsSlave[0] = ListPoints;
            
            return true;
        }
        else
        {
            ConditionsPointsSlave.clear();
//             ConditionsPointsSlave.resize(0, false);
            return false;
        }
    
        ConditionsPointsSlave.clear();
        return false;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    inline bool ExactMortarIntegrationUtility<3,3>::GetExactIntegration(    
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
        )
    {
        // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        
        // We define the tolerance
        const double Tolerance = 1.0e-8;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<PointType, 3> MasterProjectedPoint;
        array_1d<bool, 3> AllInside;
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            MasterProjectedPoint[i_node] = ContactUtilities::FastProject(SlaveCenter, MasterGeometry[i_node], SlaveNormal);
            
            GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
        
            AllInside[i_node] = SlaveGeometry.IsInside( MasterProjectedPoint[i_node].Coordinates( ), ProjectedGPLocal, Tolerance) ;
        }
        
        // We create the pointlist
        std::vector<PointType> PointList;
        
        // No point from the master is inside the slave
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false))
        {            
            // We check if all the nodes are inside the master element
            array_1d<PointType, 3> SlaveProjectedPoint;
            
            const PointType MasterCenter = MasterGeometry.Center();
            
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                SlaveProjectedPoint[i_node] = ContactUtilities::FastProject(MasterCenter, SlaveGeometry[i_node], MasterNormal);
                
                GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
            
                AllInside[i_node] = MasterGeometry.IsInside( SlaveProjectedPoint[i_node].Coordinates( ), ProjectedGPLocal, Tolerance) ;
            }
            
            // The whole slave is inside the master
            if ((AllInside[0] == true) &&
                (AllInside[1] == true) &&
                (AllInside[2] == true))
            {
                ConditionsPointsSlave.resize(1);

                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    PointType Point;
                    SlaveGeometry.PointLocalCoordinates(Point, SlaveGeometry[i_node]);
                    ConditionsPointsSlave[0][i_node] = Point;
                }
            }
            else
            {
                ConditionsPointsSlave.clear();
//                 ConditionsPointsSlave.resize(0, false);
                return false;
            }
        }
        // All the points inside
        else if ((AllInside[0] == true) &&
                 (AllInside[1] == true) &&
                 (AllInside[2] == true))
        {
            ConditionsPointsSlave.resize(1);

            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                PointType Point;
                SlaveGeometry.PointLocalCoordinates(Point, MasterProjectedPoint[i_node]);
                ConditionsPointsSlave[0][i_node] = Point;
            }
            
            return true;
        }
        else
        {            
            // Before clipping we rotate to a XY plane
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                RotatePoint(       SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                RotatePoint(MasterProjectedPoint[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                
                if (AllInside[i_node] == true)
                {
                    PointList.push_back(MasterProjectedPoint[i_node]);
                }
            }
        
            // We consider the Z coordinate constant
            const double ZRef = SlaveCenter.Coordinate(3);
            
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
                    
                    PointType IntersectedPoint;
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

                ConditionsPointsSlave.resize((ListSize - 2));
            
                // We recover this point to the triangle plane
                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    RotatePoint( SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                }
                for (unsigned int i_master = 0; i_master < MasterProjectedPoint.size(); i_master++)
                {
                    RotatePoint(MasterProjectedPoint[i_master], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                }
                
                for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
                {
                    if (FasTriagleCheck2D(MasterProjectedPoint[0], MasterProjectedPoint[IndexVector[elem] + 1], MasterProjectedPoint[IndexVector[elem + 1] + 1]) > 0.0)
                    {
                        array_1d<PointType, 3> PointsLocals;
                        
                        // Now we project to the slave surface
                        PointType gp_global_proj; 
                        PointType gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[0], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[0] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 0] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[1] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 1] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[2] = gp_local;
                        
                        ConditionsPointsSlave[elem] = PointsLocals;
                    }
                    else
                    {
                        array_1d<PointType, 3> PointsLocals;
                        
                        // Now we project to the slave surface
                        PointType gp_global_proj; 
                        PointType gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 1] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[0] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 0] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[1] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[0], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[2] = gp_local;
                        
                        ConditionsPointsSlave[elem] = PointsLocals;
                    }
                }
                
                if (ConditionsPointsSlave.size() > 0)
                {                    
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else // No intersection
            {
                // We recover the rotation
                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    RotatePoint( SlaveGeometry[i_node], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                }
                
                ConditionsPointsSlave.clear();
//                 ConditionsPointsSlave.resize(0, false);
                return false;
            }
        }
        
        ConditionsPointsSlave.clear();
        return false;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    inline bool ExactMortarIntegrationUtility<3,4>::GetExactIntegration(   
        GeometryNodeType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& MasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
        )
    {        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<PointType, 4> SlaveProjectedPoint;
        array_1d<PointType, 4> MasterProjectedPoint;
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
        
        std::vector<PointType::Pointer> DummyPointsArray (4);
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            DummyPointsArray[i_node] = boost::make_shared<PointType>(SlaveProjectedPoint[i_node]);
        }
        
        Quadrilateral2D4 <PointType> DummyQuadrilateral( DummyPointsArray );
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            GeometryNodeType::CoordinatesArrayType rResult;
            AllInside[i_node] = FasIsInsideQuadrilateral2D( DummyQuadrilateral, MasterProjectedPoint[i_node].Coordinates( ) ) ;
        }
        
        // We create the pointlist
        std::vector<PointType> PointList;
        
        // No point from the master is inside the slave
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false) &&
            (AllInside[3] == false))
        {
            // We check if all the nodes are inside the master element            
            const PointType MasterCenter = MasterGeometry.Center();
            
            const array_1d<double, 3> MasterTangentXi  = (MasterGeometry[2].Coordinates() - MasterGeometry[0].Coordinates())/norm_2(MasterGeometry[2].Coordinates() - MasterGeometry[0].Coordinates());
            const array_1d<double, 3> MasterTangentEta = MathUtils<double>::UnitCrossProduct(MasterTangentXi, MasterNormal);
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                SlaveProjectedPoint[i_node]  = ContactUtilities::FastProject( MasterCenter,  SlaveGeometry[i_node], MasterNormal);
                MasterProjectedPoint[i_node] = ContactUtilities::FastProject( MasterCenter, MasterGeometry[i_node], MasterNormal);
            }
            
            // Before clipping we rotate to a XY plane
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                RotatePoint( SlaveProjectedPoint[i_node], MasterCenter, MasterTangentXi, MasterTangentEta, false);
                RotatePoint(MasterProjectedPoint[i_node], MasterCenter, MasterTangentXi, MasterTangentEta, false);
            }
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                DummyPointsArray[i_node] = boost::make_shared<PointType>(MasterProjectedPoint[i_node]);
            }
            
            Quadrilateral2D4 <PointType> AuxDummyQuadrilateral( DummyPointsArray );
            
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                GeometryNodeType::CoordinatesArrayType rResult;
                AllInside[i_node] = FasIsInsideQuadrilateral2D( AuxDummyQuadrilateral, SlaveProjectedPoint[i_node].Coordinates( ) ) ;
            }
            
            // The whole slave is inside the master
            if ((AllInside[0] == true) &&
                (AllInside[1] == true) &&
                (AllInside[2] == true) &&
                (AllInside[3] == true))
            {
                ConditionsPointsSlave.resize(2);
                for (unsigned int i_tri = 0; i_tri < 2; i_tri++)
                {
                    for (unsigned int i_node = 0; i_node < 3; i_node++)
                    {
                        PointType Point;
                        SlaveGeometry.PointLocalCoordinates(Point, MasterProjectedPoint[i_node + i_tri]);
                        ConditionsPointsSlave[i_tri][i_node] = Point;
                    }
                }
            }
            else
            {
                ConditionsPointsSlave.clear();
    //             ConditionsPointsSlave.resize(0, false);
                return false;
            }
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
            
            // We recover this point to the triangle plane
            for (unsigned int i_master = 0; i_master < MasterProjectedPoint.size(); i_master++)
            {
                RotatePoint(MasterProjectedPoint[i_master], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
            }
            
            ConditionsPointsSlave.resize(2);
            for (unsigned int elem = 0; elem < 2; elem++) 
            {
                if (FasTriagleCheck3D(MasterProjectedPoint[0], MasterProjectedPoint[IndexVector[elem] + 1], MasterProjectedPoint[IndexVector[elem + 1] + 1]) > 0.0)
                {
                    array_1d<PointType, 3> PointsLocals;
                    
                    // Now we project to the slave surface
                    PointType gp_global_proj; 
                    PointType gp_local;
                    gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[0], - SlaveNormal); // We came back 
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                    PointsLocals[0] = gp_local;
                    gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 0] + 1], - SlaveNormal); // We came back 
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                    PointsLocals[1] = gp_local;
                    gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 1] + 1], - SlaveNormal); // We came back 
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                    PointsLocals[2] = gp_local;
                    
                    ConditionsPointsSlave[elem] = PointsLocals;
                }
                else
                {
                    array_1d<PointType, 3> PointsLocals;
                    
                    // Now we project to the slave surface
                    PointType gp_global_proj; 
                    PointType gp_local;
                    gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 1] + 1], - SlaveNormal); // We came back 
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                    PointsLocals[0] = gp_local;
                    gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 0] + 1], - SlaveNormal); // We came back 
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                    PointsLocals[1] = gp_local;
                    gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[0], - SlaveNormal); // We came back 
                    SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                    PointsLocals[2] = gp_local;
                    
                    ConditionsPointsSlave[elem] = PointsLocals;
                }                
            }
            
            if (ConditionsPointsSlave.size() > 0)
            {                                
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
            const double ZRef = SlaveCenter.Coordinate(3);
            
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
                    
                    PointType IntersectedPoint;
                    const bool intersected = Clipping2D(
                        IntersectedPoint,
                        SlaveProjectedPoint[i_edge],
                        SlaveProjectedPoint[ip_edge],
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
                
                ConditionsPointsSlave.resize((ListSize - 2));
        
                // We recover this point to the triangle plane
                for (unsigned int i_master = 0; i_master < MasterProjectedPoint.size(); i_master++)
                {
                    RotatePoint(MasterProjectedPoint[i_master], SlaveCenter, SlaveTangentXi, SlaveTangentEta, true);
                }
                
                for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
                {
                    if (FasTriagleCheck3D(MasterProjectedPoint[0], MasterProjectedPoint[IndexVector[elem] + 1], MasterProjectedPoint[IndexVector[elem + 1] + 1]) > 0.0)
                    {
                        array_1d<PointType, 3> PointsLocals;
                    
                        // Now we project to the slave surface
                        PointType gp_global_proj; 
                        PointType gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[0], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[0] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 0] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[1] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 1] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[2] = gp_local;
                        
                        ConditionsPointsSlave[elem] = PointsLocals;
                    }
                    else
                    {
                        array_1d<PointType, 3> PointsLocals;
                        
                        // Now we project to the slave surface
                        PointType gp_global_proj; 
                        PointType gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 1] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[0] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[IndexVector[elem + 0] + 1], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[1] = gp_local;
                        gp_global_proj = ContactUtilities::FastProject(SlaveCenter, MasterProjectedPoint[0], - SlaveNormal); // We came back 
                        SlaveGeometry.PointLocalCoordinates(gp_local, gp_global_proj);
                        PointsLocals[2] = gp_local;
                        
                        ConditionsPointsSlave[elem] = PointsLocals;
                    }
                }
                
                if (ConditionsPointsSlave.size() > 0)
                {                                        
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else // No intersection
            {
                ConditionsPointsSlave.clear();
//                 ConditionsPointsSlave.resize(0, false);
                return false;
            }
        }
        
        ConditionsPointsSlave.clear();
        return false;
    }
}

#endif  /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
