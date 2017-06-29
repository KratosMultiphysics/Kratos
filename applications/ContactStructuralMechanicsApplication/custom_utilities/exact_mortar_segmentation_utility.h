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
     * This utility computes the exact integration of the mortar condition (just the points, not the whole integration points)
     * @param OriginalSlaveGeometry: The geometry of the slave condition
     * @param SlaveNormal: The normal of the slave condition
     * @param OriginalMasterGeometry: The geometry of the master condition
     * @param MasterNormal: The normal of the master condition
     * @param ConditionsPointsSlave: The points that perform the exact integration
     * @param FilterFarGeometries: If the geometries are checked in first place if they are far away
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    inline bool GetExactIntegration(    
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave, 
        const bool FilterFarGeometries = true
        );
    
    /**
     * This utility computes the exact integration of the mortar condition
     * @param OriginalSlaveGeometry: The geometry of the slave condition
     * @param SlaveNormal: The normal of the slave condition
     * @param OriginalMasterGeometry: The geometry of the master condition
     * @param MasterNormal: The normal of the master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @param FilterFarGeometries: If the geometries are checked in first place if they are far away
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    inline bool GetExactIntegration(    
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave,
        const bool FilterFarGeometries = true
        )
    {
        ConditionArrayListType conditions_points_slave;
        
        const bool is_inside = GetExactIntegration(OriginalSlaveGeometry, SlaveNormal, OriginalMasterGeometry, MasterNormal, conditions_points_slave, FilterFarGeometries);
        
        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (unsigned int i_node = 0; i_node < TDim; i_node++)
            {
                PointType global_point;
                OriginalSlaveGeometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = boost::make_shared<PointType>(global_point);
            }
            
            DecompositionType decomp_geom( points_array );
            
            const GeometryType::IntegrationPointsArrayType& local_integration_slave = decomp_geom.IntegrationPoints( mAuxIntegrationMethod );
            
            // Integrating the mortar operators
            for ( unsigned int point_number = 0; point_number < local_integration_slave.size(); point_number++ )
            {
                const double weight = local_integration_slave[point_number].Weight();
                const PointType local_point_decomp = local_integration_slave[point_number].Coordinates();
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                OriginalSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);
                
                const double det_J = decomp_geom.DeterminantOfJacobian( local_point_decomp ) * (TDim == 2 ? 2.0 : 1.0);
                
                IntegrationPointsSlave.push_back( IntegrationPointType( local_point_parent.Coordinate(1), local_point_parent.Coordinate(2), weight * det_J )); // TODO: Change push_back for a fic opoeration
            }
        }
        
        return is_inside;
    }
    
    /**
     * This utility computes the exact integration of the mortar condition
     * @param SlaveCond: The slave condition
     * @param MasterCond: The master condition
     * @param CustomSolution: The matrix containing the integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    bool TestGetExactIntegration(     
        Condition::Pointer& SlaveCond,
        Condition::Pointer& MasterCond,
        Matrix& CustomSolution
        )
    {
        IntegrationPointsType integration_points_slave;
        
        const bool solution = GetExactIntegration(SlaveCond->GetGeometry(), SlaveCond->GetValue(NORMAL), MasterCond->GetGeometry(), MasterCond->GetValue(NORMAL), integration_points_slave, false);
        
        CustomSolution.resize(integration_points_slave.size(), TDim, false);
        
//         std::cout << "The Gauss Points obtained are: " << std::endl;
        for (unsigned int GP = 0; GP < integration_points_slave.size(); GP++)
        {
//             // For debug
//             KRATOS_WATCH(integration_points_slave[GP]);
            
            // Solution save:
            CustomSolution(GP, 0) = integration_points_slave[GP].Coordinate(1);
            if (TDim == 2)
            {
                CustomSolution(GP, 1) = integration_points_slave[GP].Weight();
            }
            else
            {
                CustomSolution(GP, 1) = integration_points_slave[GP].Coordinate(2);
                CustomSolution(GP, 2) = integration_points_slave[GP].Weight();
            }
        }
        
        return solution;
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
     * This function rotates to align the projected points to a parallel plane to XY
     * @param PointToRotate: The points from the origin geometry
     * @param PointReferenceRotation: The center point used as reference to rotate
     * @param SlaveNormal: The normal vector of the slave condition
     * @param slave_tangent_xi: The first tangent vector of the slave condition
     * @param slave_tangent_eta: The second tangent vector of the slave condition
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
        PointType aux_point_to_rotate;
        aux_point_to_rotate.Coordinates() = PointToRotate.Coordinates() - PointReferenceRotation.Coordinates();
        
        boost::numeric::ublas::bounded_matrix<double, 3, 3> rotation_matrix = ZeroMatrix(3, 3);
        
        if (Inversed == false)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                rotation_matrix(0, i) = SlaveTangentXi[i];
                rotation_matrix(1, i) = SlaveTangentEta[i];
            }
        }
        else
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                rotation_matrix(i, 0) = SlaveTangentXi[i];
                rotation_matrix(i, 1) = SlaveTangentEta[i];
            }
        }
        
        PointToRotate.Coordinates() = prod(rotation_matrix, aux_point_to_rotate) + PointReferenceRotation.Coordinates();
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
        const double s_orig1_orig2_x = PointOrig2.Coordinate(1) - PointOrig1.Coordinate(1);
        const double s_orig1_orig2_y = PointOrig2.Coordinate(2) - PointOrig1.Coordinate(2);
        const double s_dest1_dest2_x = PointDest2.Coordinate(1) - PointDest1.Coordinate(1);
        const double s_dest1_dest2_y = PointDest2.Coordinate(2) - PointDest1.Coordinate(2);
        
        const double denom = s_orig1_orig2_x * s_dest1_dest2_y - s_dest1_dest2_x * s_orig1_orig2_y;
    
        const double tolerance = 1.0e-12; // std::numeric_limits<double>::epsilon();
        if (std::abs(denom) < tolerance) // NOTE: Collinear
        {
            return false;
        }
        
        const double s_orig1_dest1x = PointOrig1.Coordinate(1) - PointDest1.Coordinate(1);
        const double s_orig1_dest1_y = PointOrig1.Coordinate(2) - PointDest1.Coordinate(2);
        
        const double s = (s_orig1_orig2_x * s_orig1_dest1_y - s_orig1_orig2_y * s_orig1_dest1x)/denom;
        
        const double t = (s_dest1_dest2_x * s_orig1_dest1_y - s_dest1_dest2_y * s_orig1_dest1x)/denom;
        
        if (s >= -tolerance && s <= (1.0 + tolerance) && t >= -tolerance && t <= (1.0 + tolerance))
        {
            PointIntersection.Coordinate(1) = PointOrig1.Coordinate(1) + t * s_orig1_orig2_x; 
            PointIntersection.Coordinate(2) = PointOrig1.Coordinate(2) + t * s_orig1_orig2_y; 
            
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
     * @param Axis1: The axis respect the angle is calculated
     * @param Axis2: The normal to the previous axis
     * @return angle: The angle formed
     */
    
    double AnglePoints(
        const PointType PointOrig1,
        const PointType PointOrig2,
        const array_1d<double, 3> Axis1,
        const array_1d<double, 3> Axis2
        )
    {
        array_1d<double, 3> local_edge = PointOrig2.Coordinates() - PointOrig1.Coordinates();
        if (norm_2(local_edge) > 0.0)
        {
            local_edge /= norm_2(local_edge);
        }
        
        const double xi  = inner_prod(Axis1, local_edge);
        const double eta = inner_prod(Axis2, local_edge);
        
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
        const double tolerance = 1.0e-6; // std::numeric_limits<double>::epsilon();
        
        const bool x = (std::abs(PointOrig2.X() - PointOrig1.X()) < tolerance) ? true : false;
        const bool y = (std::abs(PointOrig2.Y() - PointOrig1.Y()) < tolerance) ? true : false;
        
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
        const double tolerance = 1.0e-6; // std::numeric_limits<double>::epsilon();
        
        return (norm_2(PointOrig2.Coordinates() - PointOrig1.Coordinates()) < tolerance) ? true : false;
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
    
    double FastTriangleCheck(
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
     * @param point_list: The intersection points
     * @param all_inside: The nodes that are already known as inside the other geometry
     * @param Geometry1/Geometry2: The geometries studied (projected)
     * @param slave_tangent_xi/slave_tangent_eta: The vectors used as base to rotate
     * @param RefCenter: The reference point to rotate
     * @param IsAllInside: To simplify and consider the point_list directly
     * @return If there is intersection or not (true/false)
     */
    
    bool TriangleIntersections(
        ConditionArrayListType& ConditionsPointsSlave,
        std::vector<PointType>& point_list,
        const array_1d<bool, TNumNodes> all_inside,
        GeometryPointType& Geometry1,
        GeometryPointType& Geometry2,
        const array_1d<double, 3>& slave_tangent_xi,
        const array_1d<double, 3>& slave_tangent_eta,
        const PointType& RefCenter,
        const bool IsAllInside = false
        )
    {   
        if (IsAllInside == false)
        {
            // We consider the Z coordinate constant
            const double z_ref = RefCenter.Coordinate(3);
            
            // We find the intersection in each side
            std::map<unsigned int, unsigned int> map_edges;
            for (unsigned int i_edge = 0; i_edge < TNumNodes; i_edge++)
            {
                map_edges.insert(std::make_pair(i_edge, 0));
//                 map_edges [i_edge] = 0;
                
                const unsigned int ip_edge = (i_edge == (TNumNodes - 1)) ? 0 : i_edge + 1;
                for (unsigned int j_edge = 0; j_edge < TNumNodes; j_edge++)
                {
                    const unsigned int jp_edge = (j_edge == (TNumNodes - 1)) ? 0 : j_edge + 1;
                    
                    PointType intersected_point;
                    const bool intersected = Clipping2D(
                        intersected_point,
                        Geometry1[i_edge],
                        Geometry1[ip_edge],
                        Geometry2[j_edge],
                        Geometry2[jp_edge]
                        );
                    
                    if (intersected == true)
                    {                        
                        bool add_point = true;
                        for (unsigned int iter = 0; iter < point_list.size(); iter++)
                        {
                            if (CheckPoints2D(intersected_point, point_list[iter]) == true)
                            {
                                add_point = false;
                            }
                        }
                        
                        if (add_point == true) 
                        {
                            intersected_point.Coordinate(3) = z_ref;
                            point_list.push_back(intersected_point);
                            map_edges[i_edge] += 1;
                        }
                    }
                }
            }
            
            // No we check with edges are split just one time (which means that the corner belongs to the intersection)
            for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
            {
                unsigned int il_node = (i_node == 0) ? (TNumNodes - 1) : i_node - 1; // The first node is in edge 1 and 3
                
                if ((map_edges[i_node]  == 1) && (map_edges[il_node] == 1))
                {
                    bool add_point = true;
                    for (unsigned int iter = 0; iter < point_list.size(); iter++)
                    {
                        if (CheckPoints2D(Geometry1[i_node], point_list[iter]) == true)
                        {
                            add_point = false;
                        }
                    }
                    
                    if (add_point == true) 
                    {
                        point_list.push_back(Geometry1[i_node]);
                    }
                }
            }
        }
        
        // We compose the triangles 
        const unsigned int list_size = point_list.size();
        if (list_size > 2) // Technically the minimum is three, just in case I consider 2
        {
            // We reorder the nodes according with the angle they form with the first node
            std::vector<double> angles (list_size - 1);
            array_1d<double, 3> v = point_list[1].Coordinates() - point_list[0].Coordinates();
            v /= norm_2(v);
            array_1d<double, 3> n = GetNormalVector2D(v);
            
            for (unsigned int elem = 1; elem < list_size; elem++)
            {
                angles[elem - 1] = AnglePoints(point_list[0], point_list[elem], v, n);
                if (angles[elem - 1] < 0.0)
                {
                    v = point_list[elem].Coordinates() - point_list[0].Coordinates();
                    v /= norm_2(v);
                    n = GetNormalVector2D(v);
                    for (unsigned int aux_elem = 0; aux_elem <= (elem - 1); aux_elem++)
                    {
                        angles[aux_elem] -= angles[elem - 1];
                    }
                }
            }
            
            const std::vector<size_t> index_vector = SortIndexes<double>(angles);

            ConditionsPointsSlave.resize((list_size - 2));
        
            // We recover this point to the triangle plane
            for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
            {
                RotatePoint(Geometry1[i_node], RefCenter, slave_tangent_xi, slave_tangent_eta, true);
                RotatePoint(Geometry2[i_node], RefCenter, slave_tangent_xi, slave_tangent_eta, true);
            }
            for (unsigned int ipoint_list = 0; ipoint_list < point_list.size(); ipoint_list++)
            {
                RotatePoint(point_list[ipoint_list], RefCenter, slave_tangent_xi, slave_tangent_eta, true);
            }
            
            for (unsigned int elem = 0; elem < list_size - 2; elem++) // NOTE: We always have two points less that the number of nodes
            {
//                     // Debug
//                     PointType aux1;
//                     aux1.Coordinates() = point_list[0].Coordinates();
//                     
//                     PointType aux2;
//                     aux2.Coordinates() = point_list[index_vector[elem + 0] + 1].Coordinates();
//                     
//                     PointType aux3;
//                     aux3.Coordinates() = point_list[index_vector[elem + 1] + 1].Coordinates();
//                     
//                     std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << aux1.X() << "," << aux1.Y() << "," << aux1.Z()  << "},{" << aux2.X() << "," << aux2.Y() << "," << aux2.Z()  << "},{" << aux3.X() << "," << aux3.Y() << "," << aux3.Z()  << "}}]}],";// << std::endl;
                
                array_1d<PointType, 3> points_locals;
                    
                // Now we project to the slave surface
                PointType point_local;
                
                if (FastTriangleCheck(point_list[0], point_list[index_vector[elem] + 1], point_list[index_vector[elem + 1] + 1]) > 0.0)
                {
                    Geometry1.PointLocalCoordinates(point_local, point_list[0]);
                    points_locals[0] = point_local;
                    
                    Geometry1.PointLocalCoordinates(point_local, point_list[index_vector[elem + 0] + 1]);
                    points_locals[1] = point_local;
                    
                    Geometry1.PointLocalCoordinates(point_local, point_list[index_vector[elem + 1] + 1]);
                    points_locals[2] = point_local;
                }
                else
                {
                    Geometry1.PointLocalCoordinates(point_local, point_list[index_vector[elem + 1] + 1]);
                    points_locals[0] = point_local;

                    Geometry1.PointLocalCoordinates(point_local, point_list[index_vector[elem + 0] + 1]);
                    points_locals[1] = point_local;

                    Geometry1.PointLocalCoordinates(point_local, point_list[0]);
                    points_locals[2] = point_local;
                }
                
                ConditionsPointsSlave[elem] = points_locals;
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
//         else if(list_size == 1 || list_size == 2) // NOTE: Activate this in case you consider that your are missing something important
//         {
//             unsigned int AuxSum = 0;
//             for (unsigned int isum = 0; isum < all_inside.size(); isum++)
//             {
//                 AuxSum += all_inside[isum];
//             }
//             
//             if (AuxSum == list_size) // NOTE: One or two can be due to concident nodes on the edges
//             {
//                 ConditionsPointsSlave.clear();
// //                 ConditionsPointsSlave.resize(0, false);
//                 return false;
//             }
//             else
//             {
// //                 // Debug
// //                 KRATOS_WATCH(Geometry1);
// //                 KRATOS_WATCH(Geometry2);
// //                 for (unsigned int ipoint = 0; ipoint < list_size; ipoint++)
// //                 {
// //                     KRATOS_WATCH(point_list[ipoint]);
// //                 }
//                 
//                 // Debug (Mathematica plot!!!)
// //                 for (unsigned int isum = 0; isum < all_inside.size(); isum++)
// //                 {
// //                     KRATOS_WATCH(all_inside[isum]);
// //                 }
// //                 
// //                 PointType aux1;
// //                 aux1.Coordinates() = Geometry1[0].Coordinates();
// //                 
// //                 PointType aux2;
// //                 aux2.Coordinates() = Geometry1[1].Coordinates();
// //                 
// //                 PointType aux3;
// //                 aux3.Coordinates() = Geometry1[2].Coordinates();
// //                 
// //                 PointType aux4;
// //                 aux4.Coordinates() = Geometry2[0].Coordinates();
// //                 
// //                 PointType aux5;
// //                 aux5.Coordinates() = Geometry2[1].Coordinates();
// //                 
// //                 PointType aux6;
// //                 aux6.Coordinates() = Geometry2[2].Coordinates();
// //                 
// //                 std::cout << "Show[Graphics[{EdgeForm[Thick], Red ,Triangle[{{" << aux1.X() << "," << aux1.Y() << "},{" << aux2.X() << "," << aux2.Y() << "},{" << aux3.X() << "," << aux3.Y() << "}}]}],Graphics[{EdgeForm[Thick], Blue ,Triangle[{{" << aux4.X() << "," << aux4.Y() << "},{" << aux5.X() << "," << aux5.Y() << "},{" << aux6.X() << "," << aux6.Y() << "}}]}]";
// //                 
// //                 for (unsigned int ipoint = 0; ipoint < list_size; ipoint++)
// //                 {
// //                     std::cout << ",Graphics[{PointSize[Large],Point[{" << point_list[ipoint].X() << "," << point_list[ipoint].Y() << "}]}]";
// //                 }
// //                     
// //                 std::cout << "]" << std::endl;
//                 
//                 std::cout << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN (check if it is the edge)" << std::endl; 
// // //                 KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN" << std::endl; 
//             }
//         }
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

// TODO: Move this to a cpp file
///@name Explicit Specializations
///@{

    template<>  
    inline bool ExactMortarIntegrationUtility<2,2>::GetExactIntegration(         
        GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave,
        const bool FilterFarGeometries
        )
    {
        // First we check if the geometries are far away
        if (FilterFarGeometries == true)
        {
            if (ContactUtilities::DistanceCheck(OriginalSlaveGeometry, OriginalMasterGeometry) == true)
            {
                ConditionsPointsSlave.clear();
                return false;
            }
        }
        
        // We take the geometry GP from the core 
        const double tolerance = 1.0e-6; // std::numeric_limits<double>::epsilon();
        
        double total_weight = 0.0;
        array_1d<double,2> auxiliar_coordinates = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType projected_gp_global;
        GeometryNodeType::CoordinatesArrayType projected_gp_local;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = OriginalSlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::FastProjectDirection(OriginalMasterGeometry, OriginalSlaveGeometry[i_slave].Coordinates(), projected_gp_global, MasterNormal, -normal ); // The opposite direction
            
            const bool is_inside = OriginalMasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, tolerance );
            
            if (is_inside == true) 
            {
                if (i_slave == 0)
                {
                    auxiliar_coordinates[0] = - 1.0;
                }
                else if (i_slave == 1)
                {
                    auxiliar_coordinates[1] =   1.0;
                }
            }
        }
        
        if ((auxiliar_coordinates[0] == - 1.0 && auxiliar_coordinates[1] == 1.0) == true)
        {
            total_weight = 2.0;
        }
        else
        {
            std::vector<double> auxiliar_xi;
            for (unsigned int i_master = 0; i_master < 2; i_master++)
            {
                projected_gp_local[0] = (i_master == 0) ? -1.0 : 1.0;
                double delta_xi = (i_master == 0) ? 0.5 : -0.5;
                const bool is_inside = ContactUtilities::ProjectIterativeLine2D(OriginalSlaveGeometry, OriginalMasterGeometry[i_master].Coordinates(), projected_gp_local, SlaveNormal, tolerance, delta_xi);
                
                if (is_inside == true)
                {
                    auxiliar_xi.push_back(projected_gp_local[0]);
                }
            }
            
            if (auxiliar_xi.size() == 1 && ((auxiliar_coordinates[0] == - 1.0 || auxiliar_coordinates[1] == 1.0)))
            {
                if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) // NOTE: Equivalent to == -1.0
                {
                    auxiliar_coordinates[1] = auxiliar_xi[0];
                }
                else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) // NOTE: Equivalent to == 1.0
                {
                    auxiliar_coordinates[0] = auxiliar_xi[0];
                }
                else
                {
                    KRATOS_WATCH(auxiliar_xi[0]);
                    KRATOS_WATCH(auxiliar_coordinates[0]);
                    KRATOS_WATCH(auxiliar_coordinates[1]);
                    KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!! (TYPE 0)" << std::endl;
                }
            }
            else  if (auxiliar_xi.size() == 2)
            {
                if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) // NOTE: Equivalent to == -1.0
                {
                    auxiliar_coordinates[1] = auxiliar_xi[0] < auxiliar_xi[1] ? auxiliar_xi[1] : auxiliar_xi[0];
                }
                else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) // NOTE: Equivalent to == 1.0
                {
                    auxiliar_coordinates[0] = auxiliar_xi[0] < auxiliar_xi[1] ? auxiliar_xi[0] : auxiliar_xi[1];
                }
                else
                {
                    if (auxiliar_xi[0] < auxiliar_xi[1])
                    {
                        auxiliar_coordinates[0] = auxiliar_xi[0];
                        auxiliar_coordinates[1] = auxiliar_xi[1];
                    }
                    else
                    {
                        auxiliar_coordinates[1] = auxiliar_xi[0];
                        auxiliar_coordinates[0] = auxiliar_xi[1];
                    }
                }
            }
            else
            {
                return false; // NOTE: Giving problems
//                 KRATOS_WATCH(OriginalSlaveGeometry);
//                 KRATOS_WATCH(OriginalMasterGeometry);
//                 KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!! (TYPE 1)" << std::endl;
            }
            
            total_weight = auxiliar_coordinates[1] - auxiliar_coordinates[0];
        }
        
        if(total_weight < 0.0)
        {
            KRATOS_ERROR << "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates: "<< auxiliar_coordinates << std::endl;
        }
        else if(total_weight > 2.0)
        {
            KRATOS_ERROR << "WAAAAAAAAAAAAARNING!!!!!!!!, impossible, Weight higher than 2: "<< auxiliar_coordinates << std::endl;
        }
        
        if (total_weight > std::numeric_limits<double>::epsilon())
        {
            ConditionsPointsSlave.resize(1);
            array_1d<PointType, 2> list_points;
            list_points[0].Coordinate(1) = auxiliar_coordinates[0];
            list_points[1].Coordinate(1) = auxiliar_coordinates[1];
            ConditionsPointsSlave[0] = list_points;
            
            return true;
        }
        else
        {
            ConditionsPointsSlave.clear();
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
        ConditionArrayListType& ConditionsPointsSlave,
        const bool FilterFarGeometries
        )
    {
        // First we check if the geometries are far away
        if (FilterFarGeometries == true)
        {
            if (ContactUtilities::DistanceCheck(OriginalSlaveGeometry, OriginalMasterGeometry) == true)
            {
                ConditionsPointsSlave.clear();
                return false;
            }
        }
        
        // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        
        // We define the tolerance
        const double tolerance = 1.0e-8; // std::numeric_limits<double>::epsilon();
        
        // We define the auxiliar geometry
        std::vector<PointType::Pointer> points_array_slave  (3);
        std::vector<PointType::Pointer> points_array_master (3);
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            PointType aux_point;
            
            aux_point.Coordinates() = OriginalSlaveGeometry[i_node].Coordinates();
            points_array_slave[i_node] = boost::make_shared<PointType>(aux_point);
            
            aux_point.Coordinates() = OriginalMasterGeometry[i_node].Coordinates();
            points_array_master[i_node] = boost::make_shared<PointType>(aux_point);
        }
        
        Triangle3D3 <PointType> slave_geometry(  points_array_slave  );
        Triangle3D3 <PointType> master_geometry( points_array_master );
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType slave_center = slave_geometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> slave_tangent_xi  = (slave_geometry[1].Coordinates() - slave_geometry[0].Coordinates())/norm_2(slave_geometry[1].Coordinates() - slave_geometry[0].Coordinates());
        const array_1d<double, 3> slave_tangent_eta = MathUtils<double>::UnitCrossProduct(slave_tangent_xi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<bool, 3> all_inside;
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            master_geometry[i_node] = ContactUtilities::FastProject(slave_center, master_geometry[i_node], SlaveNormal);
        }
        
        // Before clipping we rotate to a XY plane
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            RotatePoint( slave_geometry[i_node], slave_center, slave_tangent_xi, slave_tangent_eta, false);
            RotatePoint(master_geometry[i_node], slave_center, slave_tangent_xi, slave_tangent_eta, false);
        }
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            GeometryNodeType::CoordinatesArrayType projected_gp_local;
        
            all_inside[i_node] = slave_geometry.IsInside( master_geometry[i_node].Coordinates( ), projected_gp_local, tolerance) ;
        }
        
        // We create the pointlist
        std::vector<PointType> point_list;
        
        // No point from the master is inside the slave
        if ((all_inside[0] == false) &&
            (all_inside[1] == false) &&
            (all_inside[2] == false))
        {            
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                GeometryNodeType::CoordinatesArrayType result;
                all_inside[i_node] = master_geometry.IsInside( slave_geometry[i_node].Coordinates( ), result, tolerance ) ;
            }
            
            // The whole slave is inside the master
            if ((all_inside[0] == true) &&
                (all_inside[1] == true) &&
                (all_inside[2] == true))
            {
                ConditionsPointsSlave.resize(1);

                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    RotatePoint( slave_geometry[i_node], slave_center, slave_tangent_xi, slave_tangent_eta, true);

                    PointType point;
                    slave_geometry.PointLocalCoordinates(point, slave_geometry[i_node]);
                    ConditionsPointsSlave[0][i_node] = point;
                }
                
                return true;
            }
            else
            {
                // Before clipping we rotate to a XY plane
                for (unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    if (all_inside[i_node] == true)
                    {
                        point_list.push_back(slave_geometry[i_node]);
                    }
                }
                
                return TriangleIntersections(ConditionsPointsSlave, point_list, all_inside, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
            }
        }
        // All the points inside
        else if ((all_inside[0] == true) &&
                 (all_inside[1] == true) &&
                 (all_inside[2] == true))
        {
            ConditionsPointsSlave.resize(1);
            
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                RotatePoint( master_geometry[i_node], slave_center, slave_tangent_xi, slave_tangent_eta, true);

                PointType point;
                slave_geometry.PointLocalCoordinates(point, master_geometry[i_node]);
                ConditionsPointsSlave[0][i_node] = point;
            }
            
            return true;
        }
        else
        {            
            // Before clipping we rotate to a XY plane
            for (unsigned int i_node = 0; i_node < 3; i_node++)
            {
                if (all_inside[i_node] == true)
                {
                    point_list.push_back(master_geometry[i_node]);
                }
            }
            
            return TriangleIntersections(ConditionsPointsSlave, point_list, all_inside, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
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
        ConditionArrayListType& ConditionsPointsSlave,
        const bool FilterFarGeometries
        )
    {        
        // First we check if the geometries are far away
        if (FilterFarGeometries == true)
        {
            if (ContactUtilities::DistanceCheck(OriginalSlaveGeometry, OriginalMasterGeometry) == true)
            {
                ConditionsPointsSlave.clear();
                return false;
            }
        }
        
        // We define the tolerance
        const double tolerance = 1.0e-8; // std::numeric_limits<double>::epsilon();
        
        // We define the auxiliar geometry
        std::vector<PointType::Pointer> points_array_slave  (4);
        std::vector<PointType::Pointer> points_array_master (4);
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            PointType aux_point;
            
            aux_point.Coordinates() = OriginalSlaveGeometry[i_node].Coordinates();
            points_array_slave[i_node] = boost::make_shared<PointType>(aux_point);
            
            aux_point.Coordinates() = OriginalMasterGeometry[i_node].Coordinates();
            points_array_master[i_node] = boost::make_shared<PointType>(aux_point);
        }
        
        Quadrilateral3D4 <PointType> slave_geometry(  points_array_slave  );
        Quadrilateral3D4 <PointType> master_geometry( points_array_master );
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const PointType slave_center = slave_geometry.Center();
        
        // We define the condition tangents
        const array_1d<double, 3> slave_tangent_xi  = (slave_geometry[2].Coordinates() - slave_geometry[0].Coordinates())/norm_2(slave_geometry[2].Coordinates() - slave_geometry[0].Coordinates());
        const array_1d<double, 3> slave_tangent_eta = MathUtils<double>::UnitCrossProduct(slave_tangent_xi, SlaveNormal);
        
        // No we project both nodes from the slave side and the master side
        array_1d<bool, 4> all_inside;
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            slave_geometry[i_node]  = ContactUtilities::FastProject( slave_center,  slave_geometry[i_node], SlaveNormal);
            master_geometry[i_node] = ContactUtilities::FastProject( slave_center, master_geometry[i_node], SlaveNormal);
        }
        
        // Before clipping we rotate to a XY plane
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            RotatePoint( slave_geometry[i_node], slave_center, slave_tangent_xi, slave_tangent_eta, false);
            RotatePoint(master_geometry[i_node], slave_center, slave_tangent_xi, slave_tangent_eta, false);
        }
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            GeometryNodeType::CoordinatesArrayType rResult;
            all_inside[i_node] = slave_geometry.IsInside( master_geometry[i_node].Coordinates( ), rResult, tolerance ) ;
        }
        
        // We create the pointlist
        std::vector<PointType> point_list;
        
        // No point from the master is inside the slave
        if ((all_inside[0] == false) &&
            (all_inside[1] == false) &&
            (all_inside[2] == false) &&
            (all_inside[3] == false))
        {
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                GeometryNodeType::CoordinatesArrayType rResult;
                all_inside[i_node] = master_geometry.IsInside( slave_geometry[i_node].Coordinates( ), rResult, tolerance ) ;
            }
            
            // The whole slave is inside the master
            if ((all_inside[0] == true) &&
                (all_inside[1] == true) &&
                (all_inside[2] == true) &&
                (all_inside[3] == true))
            {
                point_list.resize(4);
                for (unsigned int i_node = 0; i_node < 4; i_node++)
                {
                    point_list[i_node] = slave_geometry[i_node];
                }
                
                return TriangleIntersections(ConditionsPointsSlave, point_list, all_inside, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center, true);
            }
            else
            {
                // We add the internal nodes
                for (unsigned int i_node = 0; i_node < 4; i_node++)
                {
                    if (all_inside[i_node] == true)
                    {
                        point_list.push_back(slave_geometry[i_node]);
                    }
                }
                
                return TriangleIntersections(ConditionsPointsSlave, point_list, all_inside, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
            }
        }
        // All the points inside
        else if ((all_inside[0] == true) &&
                 (all_inside[1] == true) &&
                 (all_inside[2] == true) &&
                 (all_inside[3] == true))
        {
            point_list.resize(4);
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                point_list[i_node] = master_geometry[i_node];
            }
            
            return TriangleIntersections(ConditionsPointsSlave, point_list, all_inside, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center, true);
        }
        else
        {
            // We add the internal nodes
            for (unsigned int i_node = 0; i_node < 4; i_node++)
            {
                if (all_inside[i_node] == true)
                {
                    point_list.push_back(master_geometry[i_node]);
                }
            }
            
            return TriangleIntersections(ConditionsPointsSlave, point_list, all_inside, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
        }
        
        ConditionsPointsSlave.clear();
        return false;
    }
}

#endif  /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
