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
#include <math.h> 
#include "utilities/openmp_utils.h"
#include "contact_structural_mechanics_application_variables.h"

// The geometry of the triangle for the "tessellation"
#include "geometries/triangle_3d_3.h"

/* The integration points (we clip triangles in 3D, so with line and triangle is enought)*/
#include "integration/line_gauss_legendre_IntegrationPoints.h"
#include "integration/triangle_gauss_legendre_IntegrationPoints.h"

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
    
    typedef Point<3>                                             PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;
    
    ///Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsType;
    
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
     * @param IntegrationOrder: The integration order to consider
     */
    ExactMortarIntegrationUtility(
        const GeometryType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        unsigned int IntegrationOrder = 0
    )
    :mSlaveGeometry(SlaveGeometry),
    mSlaveNormal(SlaveNormal),
    mIntegrationOrder(IntegrationOrder)
    {
        // Setting the auxiliar integration points
        if (mIntegrationOrder == 1)
        {
            if (TDim == 2)
            {
                mAuxIntegrationMethod = Quadrature<LineGaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
            else
            {
                mAuxIntegrationMethod = Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
        }
        else if (mIntegrationOrder == 2)
        {
            if (TDim == 2)
            {
                mAuxIntegrationMethod = Quadrature<LineGaussLegendreIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
            else
            {
                mAuxIntegrationMethod = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
        }
        else if (mIntegrationOrder == 3)
        {
            if (TDim == 2)
            {
                mAuxIntegrationMethod = Quadrature<LineGaussLegendreIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
            else
            {
                mAuxIntegrationMethod = Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
        }
        else if (mIntegrationOrder == 4)
        {
            if (TDim == 2)
            {
                mAuxIntegrationMethod = Quadrature<LineGaussLegendreIntegrationPoints4, 1, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
            else
            {
                mAuxIntegrationMethod = Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
        }
        else if (mIntegrationOrder == 5)
        {
            if (TDim == 2)
            {
                mAuxIntegrationMethod = Quadrature<LineGaussLegendreIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
            else
            {
                mAuxIntegrationMethod = Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
        }
        else
        {
            if (TDim == 2)
            {
                mAuxIntegrationMethod = Quadrature<LineGaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
            else
            {
                mAuxIntegrationMethod = Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            }
        }
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
     * This utility computes the exact integration of the mortar condition
     * @param MasterGeometry: The geometry of the master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    
    bool GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        );
    
    /**
     * This function rotates to align the projected points to a parallel plane to XY
     * @param PointToRotate: The points from the origin geometry
     * @param PointReferenceRotation: The center point used as reference to rotate
     * @param Inversed: If we rotate to the XY or we recover from XY
     * @return PointRotated: The point rotated 
     */
    
    void RotatePoint(
        Point<3>& PointToRotate,
        const Point<3> PointReferenceRotation,
        const bool Inversed
        )
    {
        const double Tolerance = 1.0e-8;
        
        // We move to the (0,0,0)
        Point<3> AuxPointToRotate;
        AuxPointToRotate.Coordinates() = PointToRotate.Coordinates() - PointReferenceRotation.Coordinates();
        
        // We calculate the normal angle
        double Phi0 = std::atan(mSlaveNormal[0]/mSlaveNormal[2]) + Tolerance;
        double Phi1 = M_PI/2.0;
        
        if (Inversed == true)
        {
            Phi0 = Phi1;
            Phi1 = Phi0;
        }
        
        AuxPointToRotate.Coordinate(1) = AuxPointToRotate.Coordinate(1) * std::sin(Phi1)/std::sin(Phi0);
        AuxPointToRotate.Coordinate(2) = AuxPointToRotate.Coordinate(2) * std::sin(Phi1)/std::sin(Phi0);
        if (Inversed == false)
        {
            AuxPointToRotate.Coordinate(3) = 0.0;
        }
        else
        {
            const double Radio = std::sqrt(AuxPointToRotate.Coordinate(1) * AuxPointToRotate.Coordinate(1) + AuxPointToRotate.Coordinate(2) * AuxPointToRotate.Coordinate(2));
            AuxPointToRotate.Coordinate(3) = Radio * std::cos(Phi1);
        }
        
        PointToRotate.Coordinates() = AuxPointToRotate.Coordinates() + PointReferenceRotation.Coordinates();
    }
    
    /**
     * This function intersects two lines in a 2D plane
     * @param PointOrig: The points from the origin geometry
     * @param PointDest: The points in the destination geometry
     * @return PointIntersection: The intersection point if there is any
     * @return True if there is a intersection point, false otherwise
     */
    
    bool Clipping(
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
        
        if (Denom == 0.0) // NOTE: Collinear
        {
            return false;
        }
        
        const bool DenomPositive = (Denom > 0.0);
        
        const double SOrig1Dest1X = PointOrig2.Coordinate(1) - PointDest1.Coordinate(1);
        const double SOrig1Dest1Y = PointOrig2.Coordinate(2) - PointDest1.Coordinate(2);
        
        const double SNumer = SOrig1Orig2X * SOrig1Dest1Y - SOrig1Orig2Y * SOrig1Dest1X;
        
        if ((SNumer < 0.0) == DenomPositive) // NOTE: No collision
        {
            return false;
        }
        
        const double TNumer = SDest1Dest2X * SOrig1Dest1Y - SDest1Dest2Y * SOrig1Dest1X;
        
        if ((TNumer < 0.0) == DenomPositive) // NOTE: No collision
        {
            return false;
        }
        
        if (((SNumer > Denom) == DenomPositive) || ((TNumer > Denom) == DenomPositive))
        {
            return false;
        }
        
        // Collision detected
        const double T = TNumer/Denom;
        
        PointIntersection.Coordinate(1) = PointOrig1.Coordinate(1) + T * SOrig1Orig2X; 
        PointIntersection.Coordinate(2) = PointOrig1.Coordinate(2) + T * SOrig1Orig2Y; 
        
        return true;
    }
    
    /**
     * This function calculates in 2D the angle between two points
     * @param PointOrig: The points from the origin geometry
     * @param PointDest: The points in the destination geometry
     * @return angle: The angle formed
     */
    
    double AnglePoints(
        const Point<3> PointOrig1,
        const Point<3> PointOrig2
        )
    {
        const double Tolerance = 1.0e-8; // Avoid division by zero
        const double x = PointOrig2.Coordinate(1) - PointOrig1.Coordinate(1) + Tolerance;
        const double y = PointOrig2.Coordinate(2) - PointOrig1.Coordinate(2);
        
        const double angle = std::atan(y/x);
        
        return angle;
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

    GeometryType mSlaveGeometry;             // The geometry of the slave condition
    array_1d<double, 3> mSlaveNormal;        // The normal of the slave condition
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
    bool ExactMortarIntegrationUtility<2,2>::GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {             
        // We take the geometry GP from the core 
        const double tol = 1.0e-4; 
        
        double total_weight = 0.0;
        array_1d<double,2> coor_aux = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType projected_gp_global;
        GeometryType::CoordinatesArrayType projected_gp_local;
        double aux_dist = 0.0;
        
        // Declare the boolean of full integral
        bool full_int = true;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = mSlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::ProjectDirection(MasterGeometry, mSlaveGeometry[i_slave].Coordinates(), projected_gp_global, aux_dist, -normal ); // The opposite direction
            
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
                // NOTE: Check this!!!!
                const bool inside = ContactUtilities::ProjectIterative(mSlaveGeometry, MasterGeometry[i_master].Coordinates(), projected_gp_global);
                
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
                    KRATOS_WATCH("WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!!");
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
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = mSlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
            IntegrationPointsSlave.resize(IntegrationPoints.size());
            for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
            {
                const double weight = IntegrationPoints[PointNumber].Weight() * total_weight/2.0;
                const double xi = 0.5 * (1.0 - IntegrationPoints[PointNumber].Coordinate(1)) * coor_aux[0] 
                                + 0.5 * (1.0 + IntegrationPoints[PointNumber].Coordinate(1)) * coor_aux[1];
                
                IntegrationPointsSlave[PointNumber] = IntegrationPoint<2>( xi, weight );
            }
        }
        else
        {
//                 IntegrationPointsSlave.resize(0); // An empty std::vector
            IntegrationPointsSlave.clear(); // An empty std::vector
        }
        
//             if (IntegrationPointsSlave.size() > 0)
//             {
//                 std::cout <<  GetGeometry()[0].X() << " " << GetGeometry()[0].Y() << " " << GetGeometry()[1].X() << " " << GetGeometry()[1].Y() << std::endl;
//                 std::cout <<  MasterGeometry[0].X() << " " << MasterGeometry[0].Y() << " " << MasterGeometry[1].X() << " " << MasterGeometry[1].Y() << std::endl;
//                 KRATOS_WATCH(coor_aux);

//                 std::cout << "IntegrationPoints : " << IntegrationPointsSlave.size( ) << std::endl;
//                 for ( unsigned int i_vec = 0; i_vec < IntegrationPointsSlave.size( ); ++i_vec )
//                 {
//                     KRATOS_WATCH( IntegrationPointsSlave[i_vec] );
//                 }
//             }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    bool ExactMortarIntegrationUtility<3,3>::GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        
        // Firt we create an auxiliar plane based in the condition center and its normal
        const Point<3> SlaveCenter = mSlaveGeometry.Center();
        
        array_1d<Point<3>, 3> SlaveProjectedPoint;
        array_1d<Point<3>, 3> MasterProjectedPoint;
        array_1d<bool, 3> AllInside;
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            MasterProjectedPoint[i_node] = ContactUtilities::FastProject(MasterGeometry[i_node], SlaveCenter, mSlaveNormal);
            
            GeometryType::CoordinatesArrayType ProjectedGPLocal;
            
            AllInside[i_node] = mSlaveGeometry.IsInside( MasterProjectedPoint[i_node].Coordinates( ), ProjectedGPLocal ) ;
        }
        
        if ((AllInside[0] == false) &&
            (AllInside[1] == false) &&
            (AllInside[2] == false))
        {
            return false;
        }
        
        // Before clipping we rotate to a XY plane
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            RotatePoint(mSlaveGeometry[i_node], SlaveProjectedPoint[i_node], false);
            const Point<3> AuxPointMaster = MasterProjectedPoint[i_node];
            RotatePoint(MasterProjectedPoint[i_node], AuxPointMaster, false);
        }
        
        // We find the intersection in each side
        std::map<unsigned int, unsigned int> MapEdges;
        std::vector<Point<3>> PointList;
        for (unsigned int i_edge = 0; i_edge < 3; i_edge++)
        {
            MapEdges.insert(std::make_pair(i_edge, 0));
//             MapEdges [i_edge] = 0;
            
            unsigned int ip_edge = (i_edge == 2) ? 0 : i_edge + 1;
            for (unsigned int j_edge = 0; j_edge < 3; j_edge++)
            {
                unsigned int jp_edge = (j_edge == 2) ? 0 : j_edge + 1;
                
                Point<3> IntersectedPoint;
                const bool intersected = Clipping(
                    IntersectedPoint,
                    SlaveProjectedPoint[i_edge],
                    SlaveProjectedPoint[ip_edge],
                    MasterProjectedPoint[j_edge],
                    MasterProjectedPoint[jp_edge]
                    );
                
                if (intersected == true)
                {
                    PointList.push_back(IntersectedPoint);
                    MapEdges[i_edge] += 1;
                }
            }
        }
        
        for (unsigned int i_node = 0; i_node < 3; i_node++)
        {
            unsigned int il_node = (i_node == 0) ? 2 : i_node - 1; // The first node is in edge 1 and 3
            
            if ((MapEdges[i_node]  == 1) && (MapEdges[il_node] == 1))
            {
                PointList.push_back(SlaveProjectedPoint[i_node]);
            }
        }
        
        // We compose the triangles 
        const unsigned int ListSize = PointList.size();
        if (ListSize > 2) // Technically the minimum is three, just in case I consider 2
        {
            // We reorder the nodes according with the angle they form with the first node
            std::vector<double> Angles (ListSize - 1);
            for (unsigned int elem = 1; elem < ListSize; elem++)
            {
                Angles[elem - 1] = AnglePoints(PointList[0], PointList[elem]);
            }
            
            const std::vector<size_t> IndexVector = SortIndexes<double>(Angles);
            
            std::vector<Point<3>::Pointer> PointsArray (3);
            
            PointsArray[0] = boost::make_shared<Point<3>>(PointList[0]);
            
            // We initialize the total area
            double TotalArea = 0.0;
            
            // We initialize our auxiliar integration point vector
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = mSlaveGeometry.IntegrationPoints(mAuxIntegrationMethod);
            const size_t LocalIntegrationSize = IntegrationPoints.size();
            
            IntegrationPointsSlave.resize((ListSize - 2) * LocalIntegrationSize);
           
            for (unsigned int elem = 0; elem < ListSize - 2; elem++) // NOTE: We always have two points less that the number of nodes
            {
                // NOTE: We add 1 because we removed from the list the fisrt point
                PointsArray[1] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 1] + 1]); 
                PointsArray[2] = boost::make_shared<Point<3>>(PointList[IndexVector[elem + 2] + 1]);
                
                // We create the triangle
                Triangle3D3 <Point<3>> triangle( PointsArray );
                
                // Now we get the GP from this triangle (and weights, will be later with the total area summed)
                const double LocalArea = triangle.Area();
                TotalArea += LocalArea;
                
                // Local points should be calculated in the global space of the XY plane, then move to the plane, then invert the projection to the original geometry (that in the case of the triangle is not necessary), then we can calculate the local points which will be final coordinates 

                // TODO: Finish me
                
                for ( unsigned int PointNumber = 0; PointNumber < LocalIntegrationSize; PointNumber++ )
                {
                    const double weight = IntegrationPoints[PointNumber].Weight * LocalArea;
                    double xi  = IntegrationPoints[PointNumber].Coordinate(1);
                    double eta = IntegrationPoints[PointNumber].Coordinate(2);
                    
                    // TODO: Convert coordinates
                    
                    IntegrationPointsSlave[elem * LocalIntegrationSize + PointNumber] = IntegrationPoint<3>( xi, eta, weight );
                }
            }
            
            // We repeat the loop to ajust the weights
            for (unsigned int elem = 0; elem < ListSize - 2; elem++)
            {
                for ( unsigned int PointNumber = LocalIntegrationSize * elem; PointNumber < LocalIntegrationSize * (elem + 1); PointNumber++ )
                {
                    IntegrationPointsSlave[PointNumber].Weight() /= TotalArea;
                }
            }
        }
        else // No intersection
        {
            return false;
        }
        
        return true;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    bool ExactMortarIntegrationUtility<3,4>::GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        // Firt we create an auxiliar plane based in the condition center and its normal
        const Point<3> SlaveCenter = mSlaveGeometry.Center();
        
        // No we project both nodes from the slave side and the master side
        array_1d<Point<3>, 4> SlaveProjectedPoint;
        array_1d<Point<3>, 4> MasterProjectedPoint;
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            SlaveProjectedPoint[i_node] = ContactUtilities::FastProject(mSlaveGeometry[i_node], SlaveCenter, mSlaveNormal);
            MasterProjectedPoint[i_node] = ContactUtilities::FastProject(MasterGeometry[i_node], SlaveCenter, mSlaveNormal);
        }
        
        // Before clipping we rotate to a XY plane
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            const Point<3> AuxPointSlave = SlaveProjectedPoint[i_node];
            RotatePoint(SlaveProjectedPoint[i_node], SlaveProjectedPoint[i_node], false);
            const Point<3> AuxPointMaster = MasterProjectedPoint[i_node];
            RotatePoint(MasterProjectedPoint[i_node], AuxPointMaster, false);
        }
        
        // We find the intersection in each side
        std::map<unsigned int, unsigned int> MapEdges;
        std::vector<Point<3>> PointList;
        for (unsigned int i_edge = 0; i_edge < 4; i_edge++)
        {
            MapEdges.insert(std::make_pair(i_edge, 0));
//             MapEdges [i_edge] = 0;
            
            unsigned int ip_edge = (i_edge == 3) ? 0 : i_edge + 1;
            for (unsigned int j_edge = 0; j_edge < 4; j_edge++)
            {
                unsigned int jp_edge = (j_edge == 3) ? 0 : j_edge + 1;
                
                Point<3> IntersectedPoint;
                const bool intersected = Clipping(
                    IntersectedPoint,
                    SlaveProjectedPoint[i_edge],
                    SlaveProjectedPoint[ip_edge],
                    MasterProjectedPoint[j_edge],
                    MasterProjectedPoint[jp_edge]
                    );
                
                if (intersected == true)
                {
                    PointList.push_back(IntersectedPoint);
                    MapEdges[i_edge] += 1;
                }
            }
        }
        
        for (unsigned int i_node = 0; i_node < 4; i_node++)
        {
            unsigned int il_node = (i_node == 0) ? 3 : i_node - 1; // The first node is in edge 1 and 4
            
            if ((MapEdges[i_node]  == 1) && (MapEdges[il_node] == 1))
            {
                PointList.push_back(SlaveProjectedPoint[i_node]);
            }
        }
        
        // We compose the triangles 
        const unsigned int ListSize = PointList.size();
        if (ListSize > 2) // Technically the minimum is three, just in case I consider 2
        {
            
        }
        else // No intersection
        {
            return false;
        }
        
        return true;
    }
}

#endif  /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
