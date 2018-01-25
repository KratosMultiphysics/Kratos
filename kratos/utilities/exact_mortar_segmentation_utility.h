//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
#if !defined(KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED)
#define KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED

// System includes
#include <iostream>
#include <iomanip>
#include <math.h>

// External includes

// Project includes
#include "includes/mortar_classes.h"

// The geometry of the triangle for the "tessellation"
/* LINES */
#include "geometries/line_2d_2.h"
/* TRIANGLES */
#include "geometries/triangle_3d_3.h"
/* QUADRILATERALS */
#include "geometries/quadrilateral_3d_4.h"

// /* The integration points (we clip triangles in 3D, so with line and triangle is enought)*/
// #include "integration/line_gauss_legendre_integration_points.h"
#include "integration/triangle_gauss_legendre_integration_points.h"

/* Utilities */
#include "utilities/math_utils.h"
#include "utilities/mortar_utilities.h"

namespace Kratos {
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

typedef Point PointType;
typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryNodeType;
typedef Geometry<PointType> GeometryPointType;

///Type definition for integration methods
typedef GeometryData::IntegrationMethod IntegrationMethod;
typedef IntegrationPoint<2> IntegrationPointType;
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
template <unsigned int TDim, unsigned int TNumNodes, bool TBelong = false>
class KRATOS_API(KRATOS_CORE) ExactMortarIntegrationUtility {
   public:
    ///@name Type Definitions
    ///@{

    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N,
        typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N,
            PointBelongsQuadrilateral3D4N>::type>::type BelongType;

    typedef std::vector<array_1d<PointBelong<TNumNodes>, TDim>>
        VectorArrayPointsBelong;

    typedef std::vector<array_1d<PointType, TDim>> VectorArrayPoints;

    typedef typename std::conditional<TBelong, VectorArrayPointsBelong,
        VectorArrayPoints>::type ConditionArrayListType;

    typedef std::vector<PointBelong<TNumNodes>> VectorPointsBelong;

    typedef std::vector<PointType> VectorPoints;

    typedef typename std::conditional<TBelong, VectorPointsBelong,
        VectorPoints>::type PointListType;

    typedef array_1d<PointBelong<TNumNodes>, 3> ArrayPointsBelong;

    typedef array_1d<PointType, 3> ArrayPoints;

    typedef
        typename std::conditional<TBelong, ArrayPointsBelong, ArrayPoints>::type
            ArrayTriangleType;

    typedef Line2D2<Point> LineType;

    typedef Triangle3D3<Point> TriangleType;

    typedef typename std::conditional<TDim == 2, LineType, TriangleType>::type
        DecompositionType;

    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ExactMortarIntegrationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    /**
     * This is the default constructor
     * @param IntegrationOrder The integration order to consider
     */
    
    ExactMortarIntegrationUtility(
        const unsigned int IntegrationOrder = 0,
        const double DistanceThreshold = std::numeric_limits<double>::max()
        ) :mIntegrationOrder(IntegrationOrder),
           mDistanceThreshold(DistanceThreshold)
    {
        GetIntegrationMethod();
    }

    /// Destructor.
    virtual ~ExactMortarIntegrationUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This utility computes the exact integration of the mortar condition (just the points, not the whole integration points)
     * @param OriginalSlaveGeometry The geometry of the slave condition
     * @param SlaveNormal The normal of the slave condition
     * @param OriginalMasterGeometry The geometry of the master condition
     * @param MasterNormal The normal of the master condition
     * @param ConditionsPointsSlave The points that perform the exact integration
     * @return True if there is a common area (the geometries intersect), false otherwise
     */

    bool GetExactIntegration(GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        ConditionArrayListType& ConditionsPointsSlave);

    /**
     * This utility computes the exact integration of the mortar condition
     * @param OriginalSlaveGeometry The geometry of the slave condition
     * @param SlaveNormal The normal of the slave condition
     * @param OriginalMasterGeometry The geometry of the master condition
     * @param MasterNormal The normal of the master condition
     * @param IntegrationPointsSlave The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */

    bool GetExactIntegration(GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        IntegrationPointsType& IntegrationPointsSlave);

    /**
     * This utility computes the exact integration of the mortar condition and returns the area
     * @param OriginalSlaveGeometry The geometry of the slave condition
     * @param SlaveNormal The normal of the slave condition
     * @param OriginalMasterGeometry The geometry of the master condition
     * @param MasterNormal The normal of the master condition
     * @param rArea The total area integrated
     * @return True if there is a common area (the geometries intersect), false otherwise
     */

    bool GetExactAreaIntegration(GeometryNodeType& OriginalSlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        GeometryNodeType& OriginalMasterGeometry,
        const array_1d<double, 3>& MasterNormal, double& rArea);

    /**
     * It returns the total area inside the integration area
     * @param ConditionsPointsSlave The points that perform the exact integration
     * @param rArea The total area integrated
     */

    void GetTotalArea(GeometryNodeType& OriginalSlaveGeometry,
        ConditionArrayListType& ConditionsPointsSlave, double& rArea);

    /**
     * This utility computes the exact integration of the mortar condition
     * @param SlaveCond The slave condition
     * @param MasterCond The master condition
     * @param CustomSolution The matrix containing the integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */

    bool TestGetExactIntegration(Condition::Pointer& SlaveCond,
        Condition::Pointer& MasterCond, Matrix& CustomSolution);

    /**
     * This utility computes the exact integration of the mortar condition and returns the area
     * @param SlaveCond The slave condition
     * @return The total area integrated
     */
    
    double TestGetExactAreaIntegration(    
        ModelPart& rMainModelPart,
        Condition::Pointer& SlaveCond
        );
    
    /**
    * This method is used for debugging purposes
    * @param IndexSlave The index of the slave geometry
    * @param SlaveGeometry The slave geometry
    * @param IndexMaster The index of the master geometry
    * @param MasterGeometry The master geometry
    * @param ConditionsPointSlave The triangular decomposition
    */

    static inline void MathematicaDebug(const unsigned int IndexSlave,
        GeometryType& SlaveGeometry, const unsigned int IndexMaster,
        GeometryType& MasterGeometry,
        std::vector<array_1d<PointBelong<TNumNodes>, 3>>&
            ConditionsPointSlave) {
        typedef Triangle3D3<PointType> TriangleType;

        std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Red}],FaceForm[],"
                     "Polygon[{{";

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            std::cout << SlaveGeometry[i].X() << "," << SlaveGeometry[i].Y()
                      << "," << SlaveGeometry[i].Z();

            if (i < TNumNodes - 1)
                std::cout << "},{";
        }
        std::cout << "}}],Text[Style[" << IndexSlave << ", Tiny],{"
                  << SlaveGeometry.Center().X() << ","
                  << SlaveGeometry.Center().Y() << ","
                  << SlaveGeometry.Center().Z() << "}]}],";  // << std::endl;

        std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Blue}],FaceForm[],"
                     "Polygon[{{";

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            std::cout << MasterGeometry[i].X() << "," << MasterGeometry[i].Y()
                      << "," << MasterGeometry[i].Z();

            if (i < TNumNodes - 1)
                std::cout << "},{";
        }

        std::cout << "}}],Text[Style[" << IndexMaster << ", Tiny],{"
                  << MasterGeometry.Center().X() << ","
                  << MasterGeometry.Center().Y() << ","
                  << MasterGeometry.Center().Z() << "}]}],";  // << std::endl;

        for (unsigned int i_geom = 0; i_geom < ConditionsPointSlave.size();
             ++i_geom) {
            std::vector<PointType::Pointer> points_array(
                3);  // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (unsigned int i_node = 0; i_node < 3; ++i_node) {
                PointType global_point;
                SlaveGeometry.GlobalCoordinates(
                    global_point, ConditionsPointSlave[i_geom][i_node]);
                points_array[i_node] =
                    boost::make_shared<PointType>(global_point);
            }

            TriangleType decomp_geom(points_array);

            std::cout << "\nGraphics3D[{Opacity[.3],Triangle[{{";
            for (unsigned int i = 0; i < 3; ++i) {
                std::cout << std::setprecision(16) << decomp_geom[i].X() << ","
                          << decomp_geom[i].Y() << "," << decomp_geom[i].Z();

                if (i < 2)
                    std::cout << "},{";
            }
            std::cout << "}}]}],";  // << std::endl;
        }

        std::cout << std::endl;
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

    void GetIntegrationMethod();

    /**
     * Get the integration method to consider
     */

    GeometryNodeType::IntegrationPointsArrayType GetIntegrationTriangle();

    /**
     * This method checks if the whole array is true
     * @param AllInside The nodes that are inside or not the geometry
     * @return True if all the nodes are inside, false otherwise
     */

    static inline bool CheckAllInside(
        const array_1d<bool, TNumNodes>& AllInside) {
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            if (!AllInside[i_node])
                return false;

        return true;
    }

    /**
     * This function intersects two lines in a 2D plane
     * @param PointOrig1 The first point from the origin geometry
     * @param PointOrig2 The second point from the origin geometry
     * @param PointDest1 The first point in the destination geometry
     * @param PointDest2 The second point in the destination geometry
     * @param PointIntersection The intersection point if there is any
     * @return True if there is a intersection point, false otherwise
     */

    static inline bool Clipping2D(PointType& PointIntersection,
        const PointType& PointOrig1, const PointType& PointOrig2,
        const PointType& PointDest1, const PointType& PointDest2) {
        
        const array_1d<double, 3>& coord_point_orig1 = PointOrig1.Coordinates();
        const array_1d<double, 3>& coord_point_orig2 = PointOrig2.Coordinates();
        const array_1d<double, 3>& coord_point_dest1 = PointDest1.Coordinates();
        const array_1d<double, 3>& coord_point_dest2 = PointDest2.Coordinates();
        
        const double s_orig1_orig2_x =
            coord_point_orig2[0] - coord_point_orig1[0];
        const double s_orig1_orig2_y =
            coord_point_orig2[1] - coord_point_orig1[1];
        const double s_dest1_dest2_x =
            coord_point_dest2[0] - coord_point_dest1[0];
        const double s_dest1_dest2_y =
            coord_point_dest2[1] - coord_point_dest1[1];

        const double denom = s_orig1_orig2_x * s_dest1_dest2_y -
                             s_dest1_dest2_x * s_orig1_orig2_y;

        const double tolerance = 1.0e-12;
        //         const double tolerance = std::numeric_limits<double>::epsilon();

        if (std::abs(denom) < tolerance) // NOTE: Collinear
            return false;
        
        const double s_orig1_dest1_x = coord_point_orig1[0] - coord_point_dest1[0];
        const double s_orig1_dest1_y = coord_point_orig1[1] - coord_point_dest1[1];
        
        const double s = (s_orig1_orig2_x * s_orig1_dest1_y - s_orig1_orig2_y * s_orig1_dest1_x)/denom;
        
        const double t = (s_dest1_dest2_x * s_orig1_dest1_y - s_dest1_dest2_y * s_orig1_dest1_x)/denom;
        
        if (s >= -tolerance && s <= (1.0 + tolerance) && t >= -tolerance && t <= (1.0 + tolerance))
        {
            PointIntersection.Coordinates()[0] = coord_point_orig1[0] + t * s_orig1_orig2_x; 
            PointIntersection.Coordinates()[1] = coord_point_orig1[1] + t * s_orig1_orig2_y; 

            return true;
        }
        else
            return false;
    }

    /**
     * This function calculates in 2D the normal vector to a given one
     * @param v The vector to compute the normal 
     * @return n The normal vector
     */

    static inline array_1d<double, 3> GetNormalVector2D(
        const array_1d<double, 3>& v) {
        array_1d<double, 3> n;

        n[0] = -v[1];
        n[1] = v[0];
        n[2] = 0.0;

        return n;
    }

    /**
     * This function calculates in 2D the angle between two points
     * @param PointOrig1 The points from the origin geometry
     * @param PointOrig2 The points in the destination geometry
     * @param Axis1 The axis respect the angle is calculated
     * @param Axis2 The normal to the previous axis
     * @return angle The angle formed
     */
    
    static inline double AnglePoints(
        const PointType& PointOrig1,
        const PointType& PointOrig2,
        const array_1d<double, 3>& Axis1,
        const array_1d<double, 3>& Axis2
        )
    {
        array_1d<double, 3> local_edge = PointOrig2.Coordinates() - PointOrig1.Coordinates();
        if (norm_2(local_edge) > 0.0)
            local_edge /= norm_2(local_edge);
        
        const double xi  = inner_prod(Axis1, local_edge);
        const double eta = inner_prod(Axis2, local_edge);

        return (std::atan2(eta, xi));
    }

    /**
     * This function checks if two points are the same one
     * @param PointOrig The points from the origin geometry
     * @param PointDest The points in the destination geometry
     * @return check The check done
     */

    static inline bool CheckPoints(
        const PointType& PointOrig, const PointType& PointDest) {
        const double tolerance = std::numeric_limits<double>::epsilon();
        return (norm_2(PointDest.Coordinates() - PointOrig.Coordinates()) < tolerance) ? true : false;
    }

    /**
     * This functions calculates the determinant of a 2D triangle (using points) to check if invert the order
     * @param PointOrig1 First point
     * @param PointOrig2 Second point
     * @param PointOrig3 Third point
     * @return The DetJ
     */

    static inline double FastTriagleCheck2D(const PointType& PointOrig1,
        const PointType& PointOrig2, const PointType& PointOrig3) {
        const double x10 = PointOrig2.X() - PointOrig1.X();
        const double y10 = PointOrig2.Y() - PointOrig1.Y();

        const double x20 = PointOrig3.X() - PointOrig1.X();
        const double y20 = PointOrig3.Y() - PointOrig1.Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|	            |=	|	          |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|

        return x10 * y20 - y10 * x20;
    }

    /**
     * This function push backs the points that are inside
     * @param PointList The intersection points
     * @param AllInside The nodes that are already known as inside the other geometry
     * @param ThisGeometry The geometry considered
     */

    inline void PushBackPoints(VectorPoints& PointList,
        const array_1d<bool, TNumNodes>& AllInside,
        GeometryPointType& ThisGeometry);

    /**
     * This function push backs the points that are inside
     * @param PointList The intersection points
     * @param AllInside The nodes that are already known as inside the other geometry
     * @param ThisGeometry The geometry considered
     */

    inline void PushBackPoints(VectorPointsBelong& PointList,
        const array_1d<bool, TNumNodes>& AllInside,
        GeometryPointType& ThisGeometry, const PointBelongs& ThisBelongs);

    /**
     * This function checks if the points of Geometry2 are inside Geometry1
     * @param AllInside The nodes that are inside or not the geometry
     * @param Geometry1 The geometry where the points are checked
     * @param Geometry2 The geometry to check
     */

    inline void CheckInside(array_1d<bool, TNumNodes>& AllInside,
        GeometryPointType& Geometry1, GeometryPointType& Geometry2,
        const double Tolerance);

    /**
     * This function computes the angles indexes
     * @param PointList The intersection points
     */

    inline std::vector<std::size_t> ComputeAnglesIndexes(
        PointListType& PointList) const;

    /**
     * This function computes the angles indexes
     * @param PointList The intersection points
     * @param Geometry1 The first geometry studied (projected)
     * @param Geometry2 The second geometry studied (projected)
     * @param RefCenter The reference point to rotate
     */

    inline void ComputeClippingIntersections(PointListType& PointList,
        GeometryPointType& Geometry1, GeometryPointType& Geometry2,
        const PointType& RefCenter);

    /**
     * This function calculates the triangles intersections (this is a module, that can be used directly in the respective function)
     * @param ConditionsPointsSlave The final solution vector, containing all the nodes
     * @param PointList The intersection points
     * @param Geometry1 The first geometry studied (projected)
     * @param Geometry2 The second geometry studied (projected)
     * @param SlaveTangentXi The first vector used as base to rotate
     * @param SlaveTangentEta The second vector used as base to rotate
     * @param RefCenter The reference point to rotate
     * @param IsAllInside To simplify and consider the point_list directly
     * @return If there is intersection or not (true/false)
     */

    template <class TGeometryType = GeometryNodeType>
    inline bool TriangleIntersections(
        ConditionArrayListType& ConditionsPointsSlave, PointListType& PointList,
        TGeometryType& OriginalSlaveGeometry, GeometryPointType& Geometry1,
        GeometryPointType& Geometry2, const array_1d<double, 3>& SlaveTangentXi,
        const array_1d<double, 3>& SlaveTangentEta, const PointType& RefCenter,
        const bool IsAllInside = false);

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

    const unsigned int mIntegrationOrder;    // The integration order to consider
    const double mDistanceThreshold;         // The distance where we directly  consider out of integration limits
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
};  // Class ExactMortarIntegrationUtility
}
#endif /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
