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
#include "geometries/triangle_3d_3.h"
/* QUADRILATERALS */
#include "geometries/quadrilateral_3d_4.h"

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
    
    typedef Point<3>                                                                    PointType;
    typedef Node<3>                                                                      NodeType;
    typedef Geometry<NodeType>                                                   GeometryNodeType;
    typedef Geometry<PointType>                                                 GeometryPointType;
    
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod                                     IntegrationMethod;
    typedef IntegrationPoint<2>                                              IntegrationPointType;
    typedef GeometryNodeType::IntegrationPointsArrayType                    IntegrationPointsType;
    
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
    
    typedef std::vector<array_1d<PointType,TDim>>                          ConditionArrayListType;
    
    typedef Line2D2<Point<3>>                                                            LineType;
    
    typedef Triangle3D3<Point<3>>                                                    TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;
    
    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ExactMortarIntegrationUtility);
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Constructor
    
    /**
     * This is the default constructor
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
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        ConditionArrayListType ConditionsPointsSlave;
        
        const bool IsInside = GetExactIntegration(OriginalSlaveGeometry,SlaveNormal,OriginalMasterGeometry,MasterNormal,ConditionsPointsSlave);
        
        for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); i_geom++)
        {
            std::vector<PointType::Pointer> PointsArray (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (unsigned int inode = 0; inode < TDim; inode++)
            {
                PointType GlobalPoint;
                OriginalSlaveGeometry.GlobalCoordinates(GlobalPoint, ConditionsPointsSlave[i_geom][inode]);
                PointsArray[inode] = boost::make_shared<PointType>(GlobalPoint);
            }
            
            DecompositionType DecompGeom( PointsArray );
            
            const GeometryType::IntegrationPointsArrayType& LocalIntegrationPointsSlave = DecompGeom.IntegrationPoints( mAuxIntegrationMethod );
            
            // Integrating the mortar operators
            for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationPointsSlave.size(); PointNumber++ )
            {
                const double Weight = LocalIntegrationPointsSlave[PointNumber].Weight();
                const PointType LocalPointDecomp = LocalIntegrationPointsSlave[PointNumber].Coordinates();
                PointType LocalPointParent;
                PointType GPGlobal;
                DecompGeom.GlobalCoordinates(GPGlobal, LocalPointDecomp);
                OriginalSlaveGeometry.PointLocalCoordinates(LocalPointParent, GPGlobal);
                
                const double DetJ = DecompGeom.DeterminantOfJacobian( LocalPointDecomp ) * (TDim == 2 ? 2.0 : 1.0);
                
                IntegrationPointsSlave.push_back( IntegrationPointType( LocalPointParent.Coordinate(1), LocalPointParent.Coordinate(2), Weight * DetJ )); // TODO: Change push_back for a fic opoeration
            }
        }
        
        return IsInside;
    }
    
    /**
     * This utility computes the exact integration of the mortar condition (just the points, not the whole integration points)
     * @param MasterGeometry: The geometry of the master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    inline bool GetExactIntegration(    
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
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
     * This function checks in 2D if two points are the same one
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
     * This function checks in 2D if two points are the same one
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
     * This function calculates the triangles intersections (this is a module, that can be used directly in the respective function)
     * @param ConditionsPointsSlave: The final solution vector, containing all the nodes
     * @param PointList: The intersection points
     * @param AllInside: The nodes that are already known as inside the other geometry
     * @param Geometry1/Geometry2: The geometries studied (projected)
     * @param SlaveTangentXi/SlaveTangentEta: The vectors used as base to rotate
     * @param RefCenter: The reference point to rotate
     * @param IsAllInside: To simplify and consider the PointList directly
     * @return If there is intersection or not (true/false)
     */
    
    bool TriangleIntersections(
        ConditionArrayListType& ConditionsPointsSlave,
        std::vector<PointType>& PointList,
        const array_1d<bool, TNumNodes> AllInside,
        GeometryPointType& Geometry1,
        GeometryPointType& Geometry2,
        const array_1d<double, 3>& SlaveTangentXi,
        const array_1d<double, 3>& SlaveTangentEta,
        const PointType& RefCenter,
        const bool IsAllInside = false
        )
    {   
        if (IsAllInside == false)
        {
            // We consider the Z coordinate constant
            const double ZRef = RefCenter.Coordinate(3);
            
            // We find the intersection in each side
            std::map<unsigned int, unsigned int> MapEdges;
            for (unsigned int i_edge = 0; i_edge < TNumNodes; i_edge++)
            {
                MapEdges.insert(std::make_pair(i_edge, 0));
//                 MapEdges [i_edge] = 0;
                
                const unsigned int ip_edge = (i_edge == (TNumNodes - 1)) ? 0 : i_edge + 1;
                for (unsigned int j_edge = 0; j_edge < TNumNodes; j_edge++)
                {
                    const unsigned int jp_edge = (j_edge == (TNumNodes - 1)) ? 0 : j_edge + 1;
                    
                    PointType IntersectedPoint;
                    const bool intersected = Clipping2D(
                        IntersectedPoint,
                        Geometry1[i_edge],
                        Geometry1[ip_edge],
                        Geometry2[j_edge],
                        Geometry2[jp_edge]
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
            for (unsigned int inode = 0; inode < TNumNodes; inode++)
            {
                unsigned int il_node = (inode == 0) ? (TNumNodes - 1) : inode - 1; // The first node is in edge 1 and 3
                
                if ((MapEdges[inode]  == 1) && (MapEdges[il_node] == 1))
                {
                    bool AddPoint = true;
                    for (unsigned int iter = 0; iter < PointList.size(); iter++)
                    {
                        if (CheckPoints2D(Geometry1[inode], PointList[iter]) == true)
                        {
                            AddPoint = false;
                        }
                    }
                    
                    if (AddPoint == true) 
                    {
                        PointList.push_back(Geometry1[inode]);
                    }
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
            for (unsigned int inode = 0; inode < TNumNodes; inode++)
            {
                RotatePoint(Geometry1[inode], RefCenter, SlaveTangentXi, SlaveTangentEta, true);
                RotatePoint(Geometry2[inode], RefCenter, SlaveTangentXi, SlaveTangentEta, true);
            }
            for (unsigned int i_pointlist = 0; i_pointlist < PointList.size(); i_pointlist++)
            {
                RotatePoint(PointList[i_pointlist], RefCenter, SlaveTangentXi, SlaveTangentEta, true);
            }
            
            for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
            {
//                     // Debug
//                     PointType aux1;
//                     aux1.Coordinates() = PointList[0].Coordinates();
//                     
//                     PointType aux2;
//                     aux2.Coordinates() = PointList[IndexVector[elem + 0] + 1].Coordinates();
//                     
//                     PointType aux3;
//                     aux3.Coordinates() = PointList[IndexVector[elem + 1] + 1].Coordinates();
//                     
//                     std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                
                if (FasTriagleCheck2D(PointList[0], PointList[IndexVector[elem] + 1], PointList[IndexVector[elem + 1] + 1]) > 0.0)
                {
                    array_1d<PointType, 3> PointsLocals;
                    
                    // Now we project to the slave surface
                    PointType point_local;

                    Geometry1.PointLocalCoordinates(point_local, PointList[0]);
                    PointsLocals[0] = point_local;
                    
                    Geometry1.PointLocalCoordinates(point_local, PointList[IndexVector[elem + 0] + 1]);
                    PointsLocals[1] = point_local;
                    
                    Geometry1.PointLocalCoordinates(point_local, PointList[IndexVector[elem + 1] + 1]);
                    PointsLocals[2] = point_local;
                    
                    ConditionsPointsSlave[elem] = PointsLocals;
                }
                else
                {
                    array_1d<PointType, 3> PointsLocals;
                    
                    // Now we project to the slave surface
                    PointType point_local;

                    Geometry1.PointLocalCoordinates(point_local, PointList[IndexVector[elem + 1] + 1]);
                    PointsLocals[0] = point_local;

                    Geometry1.PointLocalCoordinates(point_local, PointList[IndexVector[elem + 0] + 1]);
                    PointsLocals[1] = point_local;
 
                    Geometry1.PointLocalCoordinates(point_local, PointList[0]);
                    PointsLocals[2] = point_local;
                    
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
        else if(ListSize == 2) // NOTE: If 1 could be because there is an intersection of just one node (concident node)
        {
            unsigned int AuxSum = 0;
            for (unsigned int iter = 0; iter < TNumNodes; iter++)
            {
                AuxSum += AllInside[iter];
            }
            if (AuxSum < 2)
            {
                KRATOS_WATCH(ListSize);
                KRATOS_WATCH(AuxSum);
//                 KRATOS_WATCH(Geometry1);
//                 KRATOS_WATCH(Geometry2);
                KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN" << std::endl; 
            }
        }
        else // No intersection
        {
            ConditionsPointsSlave.clear();
//                 ConditionsPointsSlave.resize(0, false);
            return false;
        }
        
        ConditionsPointsSlave.clear();
//         ConditionsPointsSlave.resize(0, false);
        return false;
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
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
        )
    {             
        // We take the geometry GP from the core 
        const double Tolerance = 1.0e-6;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        double total_weight = 0.0;
        array_1d<double,2> AuxiliarCoordinates = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType ProjectedGPGlobal;
        GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = OriginalSlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::FastProjectDirection(OriginalMasterGeometry, OriginalSlaveGeometry[i_slave].Coordinates(), ProjectedGPGlobal, MasterNormal, -normal ); // The opposite direction
            
            const bool IsInside = OriginalMasterGeometry.IsInside( ProjectedGPGlobal.Coordinates( ), ProjectedGPLocal, Tolerance );
            
            if (IsInside == true) 
            {
                if (i_slave == 0)
                {
                    AuxiliarCoordinates[0] = - 1.0;
                }
                else if (i_slave == 1)
                {
                    AuxiliarCoordinates[1] =   1.0;
                }
            }
        }
        
        if ((AuxiliarCoordinates[0] == - 1.0 && AuxiliarCoordinates[1] == 1.0) == true)
        {
            total_weight = 2.0;
        }
        else
        {
            std::vector<double> AuxiliarXi;
            for (unsigned int i_master = 0; i_master < 2; i_master++)
            {
                const bool IsInside = ContactUtilities::ProjectIterativeLine2D(OriginalSlaveGeometry, OriginalMasterGeometry[i_master].Coordinates(), ProjectedGPLocal, SlaveNormal);
                
                if (IsInside == true)
                {
                    AuxiliarXi.push_back(ProjectedGPLocal[0]);
                }
            }
            
            if (AuxiliarXi.size() == 1)
            {
                if (AuxiliarCoordinates[0] == - 1.0)
                {
                    AuxiliarCoordinates[1] = AuxiliarXi[0];
                }
                else if (AuxiliarCoordinates[1] == 1.0)
                {
                    AuxiliarCoordinates[0] = AuxiliarXi[0];
                }
                else
                {
                    KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!!" << std::endl;
                }
            }
            else  if (AuxiliarXi.size() == 2)
            {
                if (AuxiliarCoordinates[0] == - 1.0)
                {
                    AuxiliarCoordinates[1] = AuxiliarXi[0] < AuxiliarXi[1] ? AuxiliarXi[1] : AuxiliarXi[0];
                }
                else if (AuxiliarCoordinates[1] == 1.0)
                {
                    AuxiliarCoordinates[0] = AuxiliarXi[0] < AuxiliarXi[1] ? AuxiliarXi[0] : AuxiliarXi[1];
                }
                else
                {
                    if (AuxiliarXi[0] < AuxiliarXi[1])
                    {
                        AuxiliarCoordinates[0] = AuxiliarXi[0];
                        AuxiliarCoordinates[1] = AuxiliarXi[1];
                    }
                    else
                    {
                        AuxiliarCoordinates[1] = AuxiliarXi[0];
                        AuxiliarCoordinates[0] = AuxiliarXi[1];
                    }
                }
            }
            else
            {
                KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!!" << std::endl;
            }
            
            total_weight = AuxiliarCoordinates[1] - AuxiliarCoordinates[0];
        }
        
        if(total_weight < 0.0)
        {
            KRATOS_ERROR << "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates: "<< AuxiliarCoordinates << std::endl;
        }
        
//         // Debug
//         std::cout << "xi1 " << AuxiliarCoordinates[0] << " xi2 " << AuxiliarCoordinates[1] << std::endl;
        
        if (total_weight > std::numeric_limits<double>::epsilon())
        {
            ConditionsPointsSlave.resize(1);
            array_1d<PointType, 2> ListPoints;
            ListPoints[0].Coordinate(1) = AuxiliarCoordinates[0];
            ListPoints[1].Coordinate(1) = AuxiliarCoordinates[1];
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
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
        )
    {
        // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        
        // We define the tolerance
        const double Tolerance = 1.0e-8;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // We define the auxiliar geometry
        std::vector<PointType::Pointer> PointsArraySlave  (3);
        std::vector<PointType::Pointer> PointsArrayMaster (3);
        for (unsigned int inode = 0; inode < 3; inode++)
        {
            PointType AuxPoint;
            
            AuxPoint.Coordinates() = OriginalSlaveGeometry[inode].Coordinates();
            PointsArraySlave[inode] = boost::make_shared<PointType>(AuxPoint);
            
            AuxPoint.Coordinates() = OriginalMasterGeometry[inode].Coordinates();
            PointsArrayMaster[inode] = boost::make_shared<PointType>(AuxPoint);
        }
        
        Triangle3D3 <PointType> SlaveGeometry(  PointsArraySlave  );
        Triangle3D3 <PointType> MasterGeometry( PointsArrayMaster );
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[1].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<bool, 3> AllInside;
        
        for (unsigned int inode = 0; inode < 3; inode++)
        {
            MasterGeometry[inode] = ContactUtilities::FastProject(SlaveCenter, MasterGeometry[inode], SlaveNormal);
            
            GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
        
            AllInside[inode] = SlaveGeometry.IsInside( MasterGeometry[inode].Coordinates( ), ProjectedGPLocal, Tolerance) ;
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
            
            for (unsigned int inode = 0; inode < 3; inode++)
            {
                SlaveProjectedPoint[inode] = ContactUtilities::FastProject(MasterCenter, SlaveGeometry[inode], MasterNormal);
                
                GeometryNodeType::CoordinatesArrayType ProjectedGPLocal;
            
                AllInside[inode] = MasterGeometry.IsInside( SlaveProjectedPoint[inode].Coordinates( ), ProjectedGPLocal, Tolerance) ;
            }
            
            // The whole slave is inside the master
            if ((AllInside[0] == true) &&
                (AllInside[1] == true) &&
                (AllInside[2] == true))
            {
                ConditionsPointsSlave.resize(1);

                for (unsigned int inode = 0; inode < 3; inode++)
                {
                    PointType Point;
                    SlaveGeometry.PointLocalCoordinates(Point, SlaveGeometry[inode]);
                    ConditionsPointsSlave[0][inode] = Point;
                }
                
                return true;
            }
            else
            {
                // Before clipping we rotate to a XY plane
                for (unsigned int inode = 0; inode < 3; inode++)
                {
                    RotatePoint( SlaveGeometry[inode], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                    RotatePoint(MasterGeometry[inode], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                    
                    if (AllInside[inode] == true)
                    {
                        PointList.push_back(SlaveGeometry[inode]);
                    }
                }
                
                return TriangleIntersections(ConditionsPointsSlave, PointList, AllInside, SlaveGeometry, MasterGeometry, SlaveTangentXi, SlaveTangentEta, SlaveCenter);
            }
        }
        // All the points inside
        else if ((AllInside[0] == true) &&
                 (AllInside[1] == true) &&
                 (AllInside[2] == true))
        {
            ConditionsPointsSlave.resize(1);

            for (unsigned int inode = 0; inode < 3; inode++)
            {
                PointType Point;
                SlaveGeometry.PointLocalCoordinates(Point, MasterGeometry[inode]);
                ConditionsPointsSlave[0][inode] = Point;
            }
            
            return true;
        }
        else
        {            
            // Before clipping we rotate to a XY plane
            for (unsigned int inode = 0; inode < 3; inode++)
            {
                RotatePoint( SlaveGeometry[inode], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                RotatePoint(MasterGeometry[inode], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
                
                if (AllInside[inode] == true)
                {
                    PointList.push_back(MasterGeometry[inode]);
                }
            }
            
            return TriangleIntersections(ConditionsPointsSlave, PointList, AllInside, SlaveGeometry, MasterGeometry, SlaveTangentXi, SlaveTangentEta, SlaveCenter);
        }
        
        ConditionsPointsSlave.clear();
        return false;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    inline bool ExactMortarIntegrationUtility<3,4>::GetExactIntegration(   
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave
        )
    {        
        // We define the tolerance
        const double Tolerance = 1.0e-8;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // We define the auxiliar geometry
        std::vector<PointType::Pointer> PointsArraySlave  (4);
        std::vector<PointType::Pointer> PointsArrayMaster (4);
        for (unsigned int inode = 0; inode < 4; inode++)
        {
            PointType AuxPoint;
            
            AuxPoint.Coordinates() = OriginalSlaveGeometry[inode].Coordinates();
            PointsArraySlave[inode] = boost::make_shared<PointType>(AuxPoint);
            
            AuxPoint.Coordinates() = OriginalMasterGeometry[inode].Coordinates();
            PointsArrayMaster[inode] = boost::make_shared<PointType>(AuxPoint);
        }
        
        Quadrilateral3D4 <PointType> SlaveGeometry(  PointsArraySlave  );
        Quadrilateral3D4 <PointType> MasterGeometry( PointsArrayMaster );
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType SlaveCenter = SlaveGeometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> SlaveTangentXi  = (SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates())/norm_2(SlaveGeometry[2].Coordinates() - SlaveGeometry[0].Coordinates());
        const array_1d<double, 3> SlaveTangentEta = MathUtils<double>::UnitCrossProduct(SlaveTangentXi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<bool, 4> AllInside;
        
        for (unsigned int inode = 0; inode < 4; inode++)
        {
            SlaveGeometry[inode]  = ContactUtilities::FastProject( SlaveCenter,  SlaveGeometry[inode], SlaveNormal);
            MasterGeometry[inode] = ContactUtilities::FastProject( SlaveCenter, MasterGeometry[inode], SlaveNormal);
        }
        
        // Before clipping we rotate to a XY plane
        for (unsigned int inode = 0; inode < 4; inode++)
        {
            RotatePoint( SlaveGeometry[inode], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
            RotatePoint(MasterGeometry[inode], SlaveCenter, SlaveTangentXi, SlaveTangentEta, false);
        }
        
        for (unsigned int inode = 0; inode < 4; inode++)
        {
            GeometryNodeType::CoordinatesArrayType rResult;
            AllInside[inode] = SlaveGeometry.IsInside( MasterGeometry[inode].Coordinates( ), rResult, Tolerance ) ;
        }
        
        // We create the pointlist
        std::vector<PointType> PointList;
        
        // No point from the master is inside the slave
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false) &&
            (AllInside[3] == false))
        {
            for (unsigned int inode = 0; inode < 4; inode++)
            {
                GeometryNodeType::CoordinatesArrayType rResult;
                AllInside[inode] = MasterGeometry.IsInside( SlaveGeometry[inode].Coordinates( ), rResult, Tolerance ) ;
            }
            
            // The whole slave is inside the master
            if ((AllInside[0] == true) &&
                (AllInside[1] == true) &&
                (AllInside[2] == true) &&
                (AllInside[3] == true))
            {
                PointList.resize(4);
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    PointList[inode] = SlaveGeometry[inode];
                }
                
                return TriangleIntersections(ConditionsPointsSlave, PointList, AllInside, SlaveGeometry, MasterGeometry, SlaveTangentXi, SlaveTangentEta, SlaveCenter, true);
            }
            else
            {
                // We add the internal nodes
                for (unsigned int inode = 0; inode < 4; inode++)
                {
                    if (AllInside[inode] == true)
                    {
                        PointList.push_back(SlaveGeometry[inode]);
                    }
                }
                
                return TriangleIntersections(ConditionsPointsSlave, PointList, AllInside, SlaveGeometry, MasterGeometry, SlaveTangentXi, SlaveTangentEta, SlaveCenter);
            }
        }
        // All the points inside
        else if ((AllInside[0] == true) &&
                 (AllInside[1] == true) &&
                 (AllInside[2] == true) &&
                 (AllInside[3] == true))
        {
            PointList.resize(4);
            for (unsigned int inode = 0; inode < 4; inode++)
            {
                PointList[inode] = MasterGeometry[inode];
            }
            
            return TriangleIntersections(ConditionsPointsSlave, PointList, AllInside, SlaveGeometry, MasterGeometry, SlaveTangentXi, SlaveTangentEta, SlaveCenter, true);
        }
        else
        {
            // We add the internal nodes
            for (unsigned int inode = 0; inode < 4; inode++)
            {
                if (AllInside[inode] == true)
                {
                    PointList.push_back(MasterGeometry[inode]);
                }
            }
            
            return TriangleIntersections(ConditionsPointsSlave, PointList, AllInside, SlaveGeometry, MasterGeometry, SlaveTangentXi, SlaveTangentEta, SlaveCenter);
        }
        
        ConditionsPointsSlave.clear();
        return false;
    }
}

#endif  /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
