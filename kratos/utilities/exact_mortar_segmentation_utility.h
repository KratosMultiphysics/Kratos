//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <iostream>
#include <iomanip>
#include <cmath>

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

// /* The integration points (we clip triangles in 3D, so with line and triangle is enough)*/
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

    /// The definition of the size type
    using SizeType = std::size_t;

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
 * @class ExactMortarIntegrationUtility
 * @ingroup KratosCore
 * @brief This utility calculates the exact integration necessary for the Mortar Conditions
 * @details The utility performs a mortar segmentation in order to obtain the exact integration of the geometry intersected
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TBelong If we consider belonging of nodes or not. When you do the intersections in order to get the directional derivatives you need to know where the intersections belongs to calculate the derivatives. This says between which nodes the intersection belongs
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template <SizeType TDim, SizeType TNumNodes, bool TBelong = false, SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(KRATOS_CORE) ExactMortarIntegrationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the class
    using ClassType = ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>;
    
    /// Geometric definitions
    using GeometryType = Geometry<Node>;
    using GeometryPointType = Geometry<Point>;

    /// Type definition for integration methods
    using IntegrationMethod = GeometryData::IntegrationMethod;
    using IntegrationPointType = IntegrationPoint<2>;
    using IntegrationPointsType = GeometryType::IntegrationPointsArrayType;

    /// The type of points to be considered
    using BelongType = typename std::conditional<TNumNodes == 2,
        PointBelongsLine2D2N,
        typename std::conditional<TNumNodes == 3,
            typename std::conditional<TNumNodesMaster == 3,
                PointBelongsTriangle3D3N,
                PointBelongsTriangle3D3NQuadrilateral3D4N>::type,
            typename std::conditional<TNumNodesMaster == 3,
                PointBelongsQuadrilateral3D4NTriangle3D3N,
                PointBelongsQuadrilateral3D4N>::type>::type>::type;

    /// The definition of the point with belonging
    using PointBelongType = PointBelong<TNumNodes, TNumNodesMaster>;

    /// An array of points with belonging
    using VectorArrayPointsBelong = std::vector<array_1d<PointBelongType, TDim>>;

    /// A vector of normal points
    using VectorArrayPoints = std::vector<array_1d<Point, TDim>>;

    /// The type of array of points to be considered depending if we are interested in derivatives or not
    using ConditionArrayListType = typename std::conditional<TBelong, VectorArrayPointsBelong, VectorArrayPoints>::type;

    /// A vector of points for derivatives
    using VectorPointsBelong = std::vector<PointBelongType>;

    /// A vector of normal points
    using VectorPoints = std::vector<Point>;

    /// The type of vector of points to be considered depending if we are interested in define derivatives or not
    using PointListType = typename std::conditional<TBelong, VectorPointsBelong, VectorPoints>::type;

    /// An array of points with belonging
    using ArrayPointsBelong = array_1d<PointBelongType, 3>;

    /// An array of normal points
    using ArrayPoints = array_1d<Point, 3>;

    /// The type of array of points to be used depending if we are interested in derivatives or not
    using ArrayTriangleType = typename std::conditional<TBelong, ArrayPointsBelong, ArrayPoints>::type;

    /// The points line geometry
    using LineType = Line2D2<Point>;

    /// The points triangle geometry
    using TriangleType = Triangle3D3<Point>;

    /// The geometry that will be considered for decomposition
    using DecompositionType = typename std::conditional<TDim == 2, LineType, TriangleType>::type;

    /// The definition of the index type
    using IndexType = std::size_t;

    /// Definition of epsilon
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ExactMortarIntegrationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param IntegrationOrder The integration order to consider
     * @param DistanceThreshold The maximum distance to be considered (if too far the integration will be skiped)
     */
    ExactMortarIntegrationUtility(
        const SizeType IntegrationOrder = 0,
        const double DistanceThreshold = std::numeric_limits<double>::max(),
        const SizeType EchoLevel = 0,
        const double ZeroToleranceFactor = 1.0,
        const bool ConsiderDelaunator = false
        ) :mIntegrationOrder(IntegrationOrder),
           mDistanceThreshold(DistanceThreshold),
           mEchoLevel(EchoLevel),
           mZeroToleranceFactor(ZeroToleranceFactor),
           mConsiderDelaunator(ConsiderDelaunator)
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
     * @brief This utility computes the exact integration of the mortar condition (just the points, not the whole integration points)
     * @param rOriginalSlaveGeometry The geometry of the slave condition
     * @param rSlaveNormal The normal of the slave condition
     * @param rOriginalMasterGeometry The geometry of the master condition
     * @param rMasterNormal The normal of the master condition
     * @param rConditionsPointsSlave The points that perform the exact integration
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    bool GetExactIntegration(
        const GeometryType& rOriginalSlaveGeometry,
        const array_1d<double, 3>& rSlaveNormal,
        const GeometryType& rOriginalMasterGeometry,
        const array_1d<double, 3>& rMasterNormal,
        ConditionArrayListType& rConditionsPointsSlave
        );

    /**
     * @brief This utility computes the exact integration of the mortar condition
     * @param rOriginalSlaveGeometry The geometry of the slave condition
     * @param rSlaveNormal The normal of the slave condition
     * @param rOriginalMasterGeometry The geometry of the master condition
     * @param rMasterNormal The normal of the master condition
     * @param rIntegrationPointsSlave The integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    bool GetExactIntegration(
        const GeometryType& rOriginalSlaveGeometry,
        const array_1d<double, 3>& rSlaveNormal,
        const GeometryType& rOriginalMasterGeometry,
        const array_1d<double, 3>& rMasterNormal,
        IntegrationPointsType& rIntegrationPointsSlave
        );

    /**
     * @brief This utility computes the exact integration of the mortar condition and returns the area
     * @param rOriginalSlaveGeometry The geometry of the slave condition
     * @param rSlaveNormal The normal of the slave condition
     * @param rOriginalMasterGeometry The geometry of the master condition
     * @param rMasterNormal The normal of the master condition
     * @param rArea The total area integrated
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    bool GetExactAreaIntegration(
        const GeometryType& rOriginalSlaveGeometry,
        const array_1d<double, 3>& rSlaveNormal,
        const GeometryType& rOriginalMasterGeometry,
        const array_1d<double, 3>& rMasterNormal,
        double& rArea
        );

    /**
     * @brief It returns the total area inside the integration area
     * @param rOriginalSlaveGeometry The geometry of the slave condition
     * @param rConditionsPointsSlave The points that perform the exact integration
     * @param rArea The total area integrated
     */
    void GetTotalArea(
        const GeometryType& rOriginalSlaveGeometry,
        ConditionArrayListType& rConditionsPointsSlave,
        double& rArea
        );

    /**
     * @brief This utility computes the exact integration of the mortar condition
     * @param pSlaveCond The slave condition
     * @param pMasterCond The master condition
     * @param rCustomSolution The matrix containing the integrations points that belong to the slave
     * @return True if there is a common area (the geometries intersect), false otherwise
     */
    bool TestGetExactIntegration(
        Condition::Pointer pSlaveCond,
        Condition::Pointer pMasterCond,
        Matrix& rCustomSolution
        );


    /**
     * @brief This utility computes the exact integration of the mortar condition and returns the area
     * @param pSlaveCond The slave condition
     * @param pMasterCond The master condition
     * @return The total area integrated
     */
    double TestGetExactAreaIntegration(
        Condition::Pointer pSlaveCond,
        Condition::Pointer pMasterCond
        );


    /**
     * @brief This utility computes the exact integration of the mortar condition and returns the area
     * @param rMainModelPart The main model part
     * @param pSlaveCond The slave condition
     * @return The total area integrated
     */
    double TestGetExactAreaIntegration(
        ModelPart& rMainModelPart,
        Condition::Pointer pSlaveCond
        );

    /**
     * @brief This method is used for debugging purposes. Generates a GiD mesh to check
     * @param rMainModelPart The main model part
     * @param IOConsidered The IO considered
     */
    void TestIODebug(
        ModelPart& rMainModelPart,
        const std::string IOConsidered = "GiD"
        );

    ///@}
    ///@name  Access
    ///@{

    /**
     * @brief This method gets the current mIntegrationOrder
     * @return The integration order considered
     */
    SizeType& GetIntegrationOrder()
    {
        return mIntegrationOrder;
    }

    /**
     * @brief This method sets the current mIntegrationOrder
     * @param IntegrationOrder The integration order considered
     */
    void SetIntegrationOrder(const SizeType IntegrationOrder)
    {
        mIntegrationOrder = IntegrationOrder;
        GetIntegrationMethod();
    }

    /**
     * @brief This method gets the current mDistanceThreshold
     * @return The distance threshold considered
     */
    double& GetDistanceThreshold()
    {
        return mDistanceThreshold;
    }

    /**
     * @brief This method sets the current mDistanceThreshold
     * @param DistanceThreshold The distance threshold considered
     */
    void SetDistanceThreshold(const double DistanceThreshold)
    {
        mDistanceThreshold = DistanceThreshold;
    }

    /**
     * @brief This method gets the current mEchoLevel
     * @return The echo level considered
     */
    SizeType& GetEchoLevel()
    {
        return mEchoLevel;
    }

    /**
     * @brief This method sets the current mEchoLevel
     * @param EchoLevel The echo level considered
     */
    void SetEchoLevel(const SizeType EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    /**
     * @brief This method gets the current mZeroToleranceFactor
     * @return The zero tolerance factor considered
     */
    double& GetZeroToleranceFactor()
    {
        return mZeroToleranceFactor;
    }

    /**
     * @brief This method sets the current mZeroToleranceFactor
     * @param ZeroToleranceFactor The zero tolerance factor considered
     */
    void SetZeroToleranceFactor(const double ZeroToleranceFactor)
    {
        mZeroToleranceFactor = ZeroToleranceFactor;
    }

    /**
     * @brief This method gets the current mConsiderDelaunator
     * @return If considering delaunator
     */
    bool& GetConsiderDelaunator()
    {
        return mConsiderDelaunator;
    }

    /**
     * @brief This method sets the current mConsiderDelaunator
     * @param ConsiderDelaunator If considering delaunator
     */
    void SetConsiderDelaunator(const bool ConsiderDelaunator)
    {
        mConsiderDelaunator = ConsiderDelaunator;
    }

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

    /**
     * @brief Get the integration method to consider
     */
    void GetIntegrationMethod();

    /**
     * @brief Get the integration method to consider
     */
    GeometryType::IntegrationPointsArrayType GetIntegrationTriangle();

    /**
     * @brief This method checks if the whole array is true
     * @param rAllInside The nodes that are inside or not the geometry
     * @return True if all the nodes are inside, false otherwise
     */
    template<SizeType TSizeCheck = TNumNodes>
    static inline bool CheckAllInside(const array_1d<bool, TSizeCheck>& rAllInside)
    {
        for (IndexType i_node = 0; i_node < TSizeCheck; ++i_node)
            if (!rAllInside[i_node])
                return false;

        return true;
    }

    /**
     * @brief This function intersects two lines in a 2D plane
     * @param rPointOrig1 The first point from the origin geometry
     * @param rPointOrig2 The second point from the origin geometry
     * @param rPointDest1 The first point in the destination geometry
     * @param rPointDest2 The second point in the destination geometry
     * @param rPointIntersection The intersection point if there is any
     * @return True if there is a intersection point, false otherwise
     */
    static inline bool Clipping2D(
        Point& rPointIntersection,
        const Point& rPointOrig1,
        const Point& rPointOrig2,
        const Point& rPointDest1,
        const Point& rPointDest2
        )
    {
        const array_1d<double, 3>& r_coord_point_orig1 = rPointOrig1.Coordinates();
        const array_1d<double, 3>& r_coord_point_orig2 = rPointOrig2.Coordinates();
        const array_1d<double, 3>& r_coord_point_dest1 = rPointDest1.Coordinates();
        const array_1d<double, 3>& r_coord_point_dest2 = rPointDest2.Coordinates();

        const double s_orig1_orig2_x = r_coord_point_orig2[0] - r_coord_point_orig1[0];
        const double s_orig1_orig2_y = r_coord_point_orig2[1] - r_coord_point_orig1[1];
        const double s_dest1_dest2_x = r_coord_point_dest2[0] - r_coord_point_dest1[0];
        const double s_dest1_dest2_y = r_coord_point_dest2[1] - r_coord_point_dest1[1];

        const double denom = s_orig1_orig2_x * s_dest1_dest2_y -
                             s_dest1_dest2_x * s_orig1_orig2_y;

        const double tolerance = 1.0e-15;
//         const double tolerance = std::numeric_limits<double>::epsilon();

        if (std::abs(denom) < tolerance) // NOTE: Collinear
            return false;

        const double s_orig1_dest1_x = r_coord_point_orig1[0] - r_coord_point_dest1[0];
        const double s_orig1_dest1_y = r_coord_point_orig1[1] - r_coord_point_dest1[1];

        const double s = (s_orig1_orig2_x * s_orig1_dest1_y - s_orig1_orig2_y * s_orig1_dest1_x)/denom;

        const double t = (s_dest1_dest2_x * s_orig1_dest1_y - s_dest1_dest2_y * s_orig1_dest1_x)/denom;

        if (s >= -tolerance && s <= (1.0 + tolerance) && t >= -tolerance && t <= (1.0 + tolerance)) {
            rPointIntersection.Coordinates()[0] = r_coord_point_orig1[0] + t * s_orig1_orig2_x;
            rPointIntersection.Coordinates()[1] = r_coord_point_orig1[1] + t * s_orig1_orig2_y;

            return true;
        } else
            return false;
    }

    /**
     * @brief This function calculates in 2D the normal vector to a given one
     * @param rVector The vector to compute the normal
     * @return normal The normal vector
     */
    static inline array_1d<double, 3> GetNormalVector2D(const array_1d<double, 3>& rVector)
    {
        array_1d<double, 3> normal;

        normal[0] = -rVector[1];
        normal[1] = rVector[0];
        normal[2] = 0.0;

        return normal;
    }

    /**
     * @brief This function calculates in 2D the angle between two points
     * @param rPointOrig1 The points from the origin geometry
     * @param rPointOrig2 The points in the destination geometry
     * @param rAxis1 The axis respect the angle is calculated
     * @param rAxis2 The normal to the previous axis
     * @return angle The angle formed
     */
    static inline double AnglePoints(
        const Point& rPointOrig1,
        const Point& rPointOrig2,
        const array_1d<double, 3>& rAxis1,
        const array_1d<double, 3>& rAxis2
        )
    {
        array_1d<double, 3> local_edge = rPointOrig2.Coordinates() - rPointOrig1.Coordinates();
        if (norm_2(local_edge) > 0.0)
            local_edge /= norm_2(local_edge);

        const double xi  = inner_prod(rAxis1, local_edge);
        const double eta = inner_prod(rAxis2, local_edge);

        return (std::atan2(eta, xi));
    }

    /**
     * @brief This function checks if two points are the same one
     * @param rPointOrig The points from the origin geometry
     * @param rPointDest The points in the destination geometry
     * @return check The check done
     */
    static inline bool CheckPoints(
        const Point& rPointOrig,
        const Point& rPointDest
        )
    {
//         const double tolerance = std::numeric_limits<double>::epsilon(); // NOTE: Giving some problems, too tight
        const double tolerance = 1.0e-15;
        return rPointDest.Distance(rPointOrig) < tolerance;
    }

    /**
     * @brief This functions calculates the determinant of a 2D triangle (using points) to check if invert the order
     * @param rPointOrig1 First point
     * @param rPointOrig2 Second point
     * @param rPointOrig3 Third point
     * @return The DetJ
     */
    static inline double FastTriagleCheck2D(
        const Point& rPointOrig1,
        const Point& rPointOrig2,
        const Point& rPointOrig3
        )
    {
        const double x10 = rPointOrig2.X() - rPointOrig1.X();
        const double y10 = rPointOrig2.Y() - rPointOrig1.Y();

        const double x20 = rPointOrig3.X() - rPointOrig1.X();
        const double y20 = rPointOrig3.Y() - rPointOrig1.Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|	            |=	|	          |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|

        return x10 * y20 - y10 * x20;
    }

    /**
     * @brief This function push backs the points that are inside
     * @param rPointList The intersection points
     * @param rAllInside The nodes that are already known as inside the other geometry
     * @param rThisGeometry The geometry considered
     */
    template<SizeType TSizeCheck = TNumNodes>
    inline void PushBackPoints(
        VectorPoints& rPointList,
        const array_1d<bool, TSizeCheck>& rAllInside,
        GeometryPointType& rThisGeometry
        )
    {
        for (IndexType i_node = 0; i_node < TSizeCheck; ++i_node) {
            if (rAllInside[i_node]) {
                // We check if the node already exists
                bool add_point = true;
                for (IndexType iter = 0; iter < rPointList.size(); ++iter)
                    if (CheckPoints(rThisGeometry[i_node], rPointList[iter])) add_point = false;

                if (add_point)
                    rPointList.push_back(rThisGeometry[i_node]);
            }
        }
    }

    /**
     * @brief This function push backs the points that are inside
     * @param rPointList The intersection points
     * @param rAllInside The nodes that are already known as inside the other geometry
     * @param rThisGeometry The geometry considered
     * @param rThisBelongs This determine where it belongs each intersection
     */
    template<SizeType TSizeCheck = TNumNodes>
    inline void PushBackPoints(
        VectorPointsBelong& rPointList,
        const array_1d<bool, TSizeCheck>& rAllInside,
        GeometryPointType& rThisGeometry,
        const PointBelongs& rThisBelongs
        )
    {
        for (IndexType i_node = 0; i_node < TSizeCheck; ++i_node) {
            if (rAllInside[i_node]) {
                // We check if the node already exists
                bool add_point = true;
                for (IndexType iter = 0; iter < rPointList.size(); ++iter) {
                    if (CheckPoints(rThisGeometry[i_node], rPointList[iter])) {
                        add_point = false;
                    }
                }

                if (add_point) {
                    const IndexType initial_index = (rThisBelongs == PointBelongs::Master) ? TNumNodes : 0;
                    rPointList.push_back(PointBelong<TNumNodes, TNumNodesMaster>(rThisGeometry[i_node].Coordinates(), static_cast<BelongType>(initial_index + i_node)));
                }
            }
        }
    }

    /**
     * @brief This function checks if the points of Geometry2 are inside Geometry1
     * @param rAllInside The nodes that are inside or not the geometry
     * @param rGeometry1 The geometry where the points are checked
     * @param rGeometry2 The geometry to check
     */
    template<SizeType TSizeCheck = TNumNodes>
    inline void CheckInside(
        array_1d<bool, TSizeCheck>& rAllInside,
        GeometryPointType& rGeometry1,
        GeometryPointType& rGeometry2,
        const double Tolerance
        )
    {
        for (IndexType i_node = 0; i_node < TSizeCheck; ++i_node) {
            GeometryType::CoordinatesArrayType projected_gp_local;
            rAllInside[i_node] = rGeometry1.IsInside( rGeometry2[i_node].Coordinates( ), projected_gp_local, Tolerance) ;
        }
    }

    /**
     * @brief This function computes the angles indexes
     * @param rPointList The intersection points
     * @param rNormal The normal vector
     */
    inline std::vector<IndexType> ComputeAnglesIndexes(
        PointListType& rPointList,
        const array_1d<double, 3>& rNormal
        ) const;

    /**
     * @brief This function computes the angles indexes
     * @param rPointList The intersection points
     * @param rSlaveGeometry The first (slave) geometry studied (projected)
     * @param rMasterGeometry The second (master) geometry studied (projected)
     * @param rRefCenter The reference point to rotate
     */
    inline void ComputeClippingIntersections(
        PointListType& rPointList,
        const GeometryPointType& rSlaveGeometry,
        const GeometryPointType& rMasterGeometry,
        const Point& rRefCenter
        );

    /**
     * @brief This function calculates the triangles intersections (this is a module, that can be used directly in the respective function)
     * @param rConditionsPointsSlave The final solution vector, containing all the nodes
     * @param rPointList The intersection points
     * @param rSlaveGeometry The first (slave) geometry studied (projected)
     * @param rMasterGeometry The second (master) geometry studied (projected)
     * @param rSlaveTangentXi The first vector used as base to rotate
     * @param rSlaveTangentEta The second vector used as base to rotate
     * @param rRefCenter The reference point to rotate
     * @param IsAllInside To simplify and consider the point_list directly
     * @return If there is intersection or not (true/false)
     */
    template <class TGeometryType = GeometryType>
    inline bool TriangleIntersections(
        ConditionArrayListType& rConditionsPointsSlave,
        PointListType& rPointList,
        const TGeometryType& rOriginalSlaveGeometry,
        const GeometryPointType& rSlaveGeometry,
        const GeometryPointType& rMasterGeometry,
        const array_1d<double, 3>& rSlaveTangentXi,
        const array_1d<double, 3>& rSlaveTangentEta,
        const Point& rRefCenter,
        const bool IsAllInside = false
        );

    /**
     * @brief This method checks if the center of the geometry is inside the slave geometry (to prevent convex geometries)
     * @param AuxiliarCenterLocalCoordinates These are the local coordinates corresponding to the center
     * @param NumNodes The number of nodes of the geometry
     * @return True if is inside false otherwise
     */
    static inline bool CheckCenterIsInside(
        const array_1d<double, 2>& rAuxiliarCenterLocalCoordinates,
        const SizeType NumNodes = TNumNodes
        );

    ///@}
private:
    ///@name Member Variables
    ///@{

    SizeType mIntegrationOrder;              /// The integration order to consider
    double mDistanceThreshold;               /// The distance where we directly  consider out of integration limits
    SizeType mEchoLevel;                     /// The echo level considered
    double mZeroToleranceFactor;             /// The zero tolerance factor considered
    IntegrationMethod mAuxIntegrationMethod; /// The auxiliar list of Gauss Points taken from the geometry
    bool mConsiderDelaunator;                /// If consider the DelaunatorUtilities in 3D in order to construct the triangles

    ///@}
};  // Class ExactMortarIntegrationUtility
}
