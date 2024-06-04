//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H
#define TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H

//// STL includes
#include <set>
#include <vector>
//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/triangle_mesh.hpp"

namespace queso
{
///@name QuESo Classes
///@{

class BRepOperator;
class TrimmedDomain;

/**
 * @class  TrimmedDomainOnPlane
 * @author Manuel Messmer
 * @brief Projects the trimmed domain onto a plane of the AABB and provides functions to triangluate
 *        the trimmed domain.
 * @details 1. Stores vertices and edges in different containers according to the
 *             orientation of the corresponding edge (see: CollectEdgesOnPlane(rTriangleMesh)).
 *             Edges are retrieved from EdgesOnPlane stored in rTriangleMesh).
 *          2. New edges are introduced, if neccessary, to close polyline (see. CloseContourEdges()).
 *          3. Domain on plane is triangulated (see: TriangulateDomain()).
 *          All three functions are called in pGetTriangulation().
 */
class TrimmedDomainOnPlane
{
public:
    ///@name Enum's
    ///@{

    enum Orientation {Positive, Negative, Vertical};

    ///@}
    ///@name Type Definitions
    ///@{

    typedef enum Orientation OrientationType;
    typedef std::array<double, 2> Point2DType;
    typedef Vector3d Point3DType;
    typedef Unique<TriangleMeshInterface> TriangleMeshPtrType;

    struct PointComparison {
        PointComparison(double Tolerance) : mTolerance(Tolerance){}
        bool operator() (const std::vector<Point2DType>::iterator& rLhs, const std::vector<Point2DType>::iterator& rRhs) const {
            if( std::abs( (*rLhs)[0] - (*rRhs)[0] ) < mTolerance ){ // If equal
                return (*rLhs)[1] < (*rRhs)[1] - mTolerance;
            }
            else {
                return (*rLhs)[0] <= (*rRhs)[0] - mTolerance;
            }
        }
    private:
        double mTolerance;
    };
    // Point set to enable fast search of dublicate vertices.
    typedef std::set<std::vector<Point2DType>::iterator, PointComparison> Point2DSetType;

    /**
     * @class  Edge2D
     * @author Manuel Messmer
     * @brief Simple edge with 2 vertices in 2D.
     */
    class Edge2D {
    public:
        ///@name Life Cycle
        ///@{
        /// Contructor
        Edge2D(IndexType V1, IndexType V2, const Point2DType& rNormal) : mV1(V1), mV2(V2), mNormal(rNormal) {}
        ///@}
        ///@name Public Operations
        ///@{
        IndexType V1() const {return mV1;}
        IndexType V2() const {return mV2;}
        void Set(IndexType V1, IndexType V2){mV1 = V1; mV2 = V2;}
        void SetVerticesOnUpperBoundary(bool V1, bool V2){ mVertexOnUpperBoundary = std::make_pair(V1, V2); }
        std::pair<bool, bool> IsVertexOnUpperBoundary() const { return mVertexOnUpperBoundary; }
        const std::vector<Point2DType> &GetSplitPoints() const {return mSplitPoints;}
        std::vector<Point2DType>& GetSplitPoints(){return mSplitPoints;}
        void ClearSplitPoints(){mSplitPoints.clear();}
        void AddSplitPoint(const Point2DType &rPoint){mSplitPoints.push_back(rPoint);}
        void Reserve(IndexType NewSize){mSplitPoints.reserve(NewSize);}
        Point2DType Normal() const {return mNormal; }
        ///@}
    private:
        ///@name Private Members
        ///@{
        IndexType mV1, mV2;
        Point2DType mNormal;
        std::vector<Point2DType> mSplitPoints{};
        std::pair<bool, bool> mVertexOnUpperBoundary{};
        ///@}
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TrimmedDomainOnPlane(IndexType PlaneIndex, bool UpperBoundary, const Point3DType &LowerBound, const Point3DType &UpperBound, const TrimmedDomain* pTrimmedDomain, bool SwitchOrientation = false)
        : mUpperBoundary(UpperBoundary), mLowerBound(LowerBound), mUpperBound(UpperBound), mpTrimmedDomain(pTrimmedDomain), mSwitchOrientation(SwitchOrientation)
    {
        // Orientation
        //
        //      _______                  DIRINDEX2
        //     /      /|                ´|`
        //    /_____ / |<-- lower plane  |-->DIRINDEX1
        //    |     |  /                /
        //    |     |</-- upper plane  DIRINDEX3
        //    |_____|/
        //
        mPlaneIndex = PlaneIndex;
        if (PlaneIndex == 2) {
            if( mSwitchOrientation ){
                DIRINDEX1 = 1;
                DIRINDEX2 = 0;
            } else {
                DIRINDEX1 = 0;
                DIRINDEX2 = 1;
            }
            DIRINDEX3 = 2;
        }
        else if (PlaneIndex == 1) {
            if( mSwitchOrientation ){
                DIRINDEX1 = 0;
                DIRINDEX2 = 2;
            } else {
                DIRINDEX1 = 2;
                DIRINDEX2 = 0;
            }
            DIRINDEX3 = 1;
        }
        else if (PlaneIndex == 0) {
            if( mSwitchOrientation ){
                DIRINDEX1 = 2;
                DIRINDEX2 = 1;
            } else {
                DIRINDEX1 = 1;
                DIRINDEX2 = 2;
            }
            DIRINDEX3 = 0;
        }
        else {
            QuESo_ERROR << "Wrong PlaneIndex.\n";
        }

        mSnapTolerance = RelativeSnapTolerance(mLowerBound, mUpperBound);

        // Instantiate vertices sets with mSnapTolerance.
        mVerticesSetPositive = MakeUnique<Point2DSetType>(PointComparison(mSnapTolerance));
        mVerticesSetNegative = MakeUnique<Point2DSetType>(PointComparison(mSnapTolerance));
        mVerticesSetVertical = MakeUnique<Point2DSetType>(PointComparison(mSnapTolerance));
    }

    ///@}
    ///@name Public Operations
    ///@{

    ///@brief Returns a triangulated mesh of trimmed domain.
    ///@param rTriangleMesh Clipped mesh inside trimmed domain. Must hold the edges on planes. No triangles on planes are allowed.
    ///@return TriangleMeshPtrType.
    ///@details see CollectEdgesOnPlane(), CloseContourEdges(), TriangulateDomain().
    //
    //     a_______b                 y
    //     /      /|                ´|`
    //  c /_____d/ |<-- plane lower  |-->x
    //    |     |  /                /
    //    |     |</-- plane upper  Z
    //    |_____|/
    //
    TriangleMeshPtrType pGetTriangulation(const TriangleMeshInterface &rTriangleMesh, const BRepOperator* pOperator) {
        CollectEdgesOnPlane(rTriangleMesh);
        CloseContourEdges(pOperator);
        return TriangulateDomain();
    }

    /// @brief Reserve all containers.
    /// @param NewSize
    void Reserve(IndexType NewSize) {
        mEdgesPositiveOriented.reserve(NewSize);
        mEdgesNegativeOriented.reserve(NewSize);
        mEdgesVertical.reserve(NewSize);
        mVerticesPositive.reserve(NewSize);
        mVerticesNegative.reserve(NewSize);
        mVerticesVertical.reserve(NewSize);
    }

private:

    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Coolects all edges and stores them in local containers. Sorts edges corresponding to their orientation.
    ///        Positive oriented, Negative oriented, and Vertically oriented.
    /// @details
    ///             ^                                 DIRINDEX2
    ///        -----|----- Positive oriented.             ^
    ///                                                   |---> DIRINDEX1
    ///        -----|----- Negative oriented.
    ///             V
    /// @param rTriangleMesh Clipped triangle mesh inside trimmed domain.
    void CollectEdgesOnPlane(const TriangleMeshInterface &rTriangleMesh);

    ///@brief Refine edges on trimmed domain on plane. See @details for more information.
    ///@details 1.) Checks if positive oriented edges span the entire AABB (we consider only the part that is inside the domain) in DIRINDEX1.
    ///             If not new edges are introduced at top of AABB (UpperBound[DIRINDEX2]) to fill gaps.
    ///        UB2                  UP2       _____
    ///         |  out  /   b |      |  out  /     |  new edges introduced since Point b is inside domain.
    ///         |------/      | ---> |------/      |
    ///         |  inside     |      |  inside     |
    ///        LB12          UB1    LB12           UP1
    ///
    ///         UB2
    ///         |  in   /   b |      |  in   /     |  no new edges introduced since Point b is outside.
    ///         |------/      | ---> |------/      |       LB1 - lower bound of AABB in DIRINDEX1
    ///         |  outsid     |      |  outside    |       UP1 - upper bound of AABB in DIRINDEX1
    ///        LB12          UB1    LB12           UP1
    ///
    ///         2.) Maps all negative oriented points onto positive oriented edges if overlap and vice versa.
    ///             Such that the vertices of postive and negative oriented edges are at same position in: DIRINDEX1.
    ///                 x--^-----x---^--x   positive oriented edges
    ///                 |  |     |   |                                     DIRINDEX2
    ///                 |  |     |   |                                         ^
    ///            x----v--x-----v---x      negative oriented edges            |---> DIRINDEX1
    void CloseContourEdges(const BRepOperator* pOperator);

    ///@brief Triangulates the trimmed domain on plane (see @details for more information).
    ///@details Assumes that the vertices of positive and negative oriented edges coincide along DIRINDEX1. See RefineEdges().
    ///
    ///          x___x   Positive Oriented      DIRINDEX2
    ///         /|   |                              ^
    ///   x___x/ |   |                              |----> DIRINDEX1
    ///   |   |  |   |
    ///   |   |  |   |
    ///   x___x__x___x    Negative Oriented      x-Vertices
    ///
    /// Each vertical stripe is individually triangulated. First we generate a polygon with 4 vertices and then triangulate it.
    /// If vertices of positve and negative edges coincide, only one triangle is created.
    ///
    ///          x___x   Positive Oriented
    ///         /|   |
    ///       x/_|___|   Negative Oriented       x-Vertices
    ///       Tri Poly
    ///
    ///@return TriangleMeshPtrType.
    TriangleMeshPtrType TriangulateDomain() const;

    /// @brief Returns intersecting edges with upper bound (mUpperBound[DIRINDEX2]). If at least one vertex is on mUpperBound[DIRINDEX2]
    ///        edge is added. Edge vertices are marked as On_Boundary. See: Edge2D::IsVertexOnUpperBoundary.
    ///        Also, edges are sorted from left to right along DIRINDEX1.
    /// @param [out] rEdges
    /// @param Orientation
    void FindIntersectingEdgesWithUpperBound(std::vector<Edge2D> &rEdges, OrientationType Orientation);

    ///@brief Return EdgeId of negative oriented edge that has same start and end point in DIRINDEX1.
    ///       If multiple edges are found, Id of closest partner is returned.
    ///       Return -1 if no edge is found.
    ///@param rV1 Left Point
    ///@param rV2 Right Point
    ///@param rNormal Normal of edge (rV1 and rV2).
    ///@param int EdgeId
    int FindNegativePartnerEdge(const Point2DType &rV1, const Point2DType &rV2, const Point2DType &rNormal) const;

    ///@brief Find intersecting point on edge with Orientation==OrientationDest with x=Point[DIRINDEX1] and mark as split point.
    ///@param rPoint Potential split point.
    ///@param OrientationDest Orientation of edges to split.
    void SetSplitPoint(const Point2DType &rPoint, OrientationType OrientationDest);

    ///@brief Splits all edges at their defined split point.
    ///@param Orientation Orientation of edges.
    void SplitEdgesAtSplitPoint(OrientationType Orientation);

    ////////////////////////
    /// Setter Functions ///
    ////////////////////////

    /// @brief Inserts edge into local containers.
    /// @param rV1 Point1 (3D-Point).
    /// @param rV2 Point2 (3D-Point).
    /// @param rNormal Corresponsing normal direction.
    void InsertEdge(const Point3DType& rV1, const Point3DType& rV2, const Point3DType &rNormal);

    /// @brief Inserts edge into local containers.
    /// @param rV1 Point1 (2D-Point).
    /// @param rV2 Point1 (2D-Point).
    /// @param rOrientation Current orientation.
    /// @return True if edge is inserted.
    bool InsertEdge(const Point2DType& rV1, const Point2DType& rV2, const Point2DType& rNormal, OrientationType rOrientation );

    ///@brief Simple and fast insertion of vertex to container. Check if point is already contained is omitted.
    ///@brief If uniqueness check of point is required, use: GetUniqueVertexIDs() + InsertVertex(Point2DType, IndexType, Orientation).
    ///@param rPoint NewPoint
    ///@param OrientationType Current orientation.
    ///@param Return Vertex index in container.
    IndexType InsertVertex(const Point2DType &rPoint, OrientationType Orientation);

    ///@brief Inserts point to vertex container by index. If index is known through GetUniqueVertexIDs, this function is used
    ///       to insert point at corresponding position to std::vector<Point2DType> containers.
    ///@param rPoint NewPoint
    ///@param NewIndex Index of new point.
    ///@param OrientationType Current orientation.
    void InsertVertex(const Point2DType &rPoint, IndexType NewIndex, OrientationType Orientation);

    /// @brief Adds vertices of pEdge that are marked as (IsVertexOnUpperBoundary) to container.
    ///        If left vertex is on upper boundary. Vertex is marked as right edge (std::pair.second = false).
    ///        If right vertex is on upper boundary. Vertex is marked as left edge (std::pair.second = true).
    /// @param pEdge Contains intersected vertices (see: IsVertexOnUpperBoundary).
    /// @param rVertices[out]
    void AddIntersectedVertexPositive(Edge2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);

    /// @brief Adds vertices of pEdge that are marked as (IsVertexOnUpperBoundary) to container.
    ///        If left vertex is on upper boundary. Vertex is marked as left edge (std::pair.second = true).
    ///        If right vertex is on upper boundary. Vertex is marked as right edge (std::pair.second = false).
    /// @param pEdge Contains intersected vertices (see: IsVertexOnUpperBoundary).
    /// @param rVertices[out]
    void AddIntersectedVertexNegative(Edge2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);

    /// @brief Adds vertices of pEdge that are marked as (IsVertexOnUpperBoundary) to container.
    ///        If normal[0] < ZEROTOL. Vertex is marked as left edge (std::pair.second = true).
    ///        If normal[0] > ZEROTOL. Vertex is marked as right edge (std::pair.second = false).
    /// @param pEdge Contains intersected vertices (see: IsVertexOnUpperBoundary).
    /// @param rVertices[out]
    void AddIntersectedVertexVertical(Edge2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);

    /// @brief Remove dublicate edges that are at the same position (DIRINDEX1) and have same normal.
    ///        If there are still two edges at the same position, remove both of them.
    /// @param rEdges
    /// @param rVertices
    void RemoveDublicateVerticalEdges(std::vector<Edge2D>& rEdges, std::vector<Point2DType>& rVertices);

    /// @brief Remove dublicate edges that are at the same position v1[0] and v2[0].
    /// @param rEdges
    /// @param rVertices
    void RemoveDublicateEdges(std::vector<Edge2D>& rEdges, std::vector<Point2DType>& rVertices);

    ////////////////////////
    /// Getter Functions ///
    ////////////////////////

    ///@brief Returns vertex IDs of point container, for two points.
    ///       Returns std::pair(0,0), if rV1==rV2 (according to PointComparison()).
    ///       If rV1 or rV2 are new points, a new index is created.
    ///@details std::vector<Point2DType> are the actual vertex containers. This function returns, the index of those vectors.
    ///         std::set<std::vector<Point2DType>::iterator> is used to allow a fast search of duplicate vertices.
    ///@param rV1 Vertex1.
    ///@param rV2 Vertex2.
    ///@param Orientation Current orientation.
    std::pair<IndexType, IndexType> GetUniqueVertexIDs(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation) const;

    /// @brief Returns vertex V1 of edge.
    /// @param EdgeId
    /// @param Orientation
    /// @return const Point2DType&
    const Point2DType& V1byEdgeId(IndexType EdgeId, OrientationType Orientation) const;

    /// @brief Returns vertex V2 of edge.
    /// @param EdgeId
    /// @param Orientation
    /// @return const Point2DType&
    const Point2DType& V2byEdgeId(IndexType EdgeId, OrientationType Orientation) const;

    ///@brief Returns Edges container. (const version)
    ///@param Orientation orientation of edges
    ///@return const std::vector<Edge2D>&
    const std::vector<Edge2D>& GetEdges(OrientationType Orientation) const;

    ///@brief Returns Edges container.  (non const version)
    ///@param Orientation orientation of edges
    ///@return std::vector<Edge2D>&
    std::vector<Edge2D>& GetEdges(OrientationType Orientation);

    ///@brief Returns vertices container. (const version)
    ///@param Orientation Orientation of vertices
    ///@return const std::vector<Point2DType>&
    const std::vector<Point2DType>& GetVertices(OrientationType Orientation) const;

    ///@brief Returns vertices container.  (non-const version)
    ///@param Orientation Orientation of vertices
    ///@return std::vector<Point2DType>&
    std::vector<Point2DType>& GetVertices(OrientationType Orientation);

    ///@brief Returns vertices set.  (non-const version)
    ///@param Orientation Orientation of vertices
    ///@return Point2DSetType&
    Point2DSetType& GetVerticesSet(OrientationType Orientation);

    ///@brief Returns vertices set.  (const version)
    ///@param Orientation Orientation of vertices
    ///@return Point2DSetType&
    const Point2DSetType& GetVerticesSet(OrientationType Orientation) const;

    ///@brief Returns number of edges.
    ///@param Positive Orientation of edges.
    ///@return IndexType
    IndexType GetNumberEdges(OrientationType Orientation) const;

    /// @brief Returns offset of plane to origin.
    /// @return double
    double GetPlanePosition() const;

    /// @brief Returns max value (DIRINDEX2) of all vertices in rPoints
    /// @param rPoints
    /// @return double
    double GetMaxDIRINDEX2(std::vector<Point2DType>& rPoints) const;

    ///@brief Return true if 2D point exists already in given PointSet.
    ///@param rPoint TestPoint
    ///@param rPointSet PointSet to search in.
    ///@return bool
    bool PointExists(const Point2DType& rPoint , const Point2DSetType& rPointSet) const;

    ///@}
    ///@name Private Members
    ///@{

    /// Edges container
    std::vector<Edge2D> mEdgesPositiveOriented{};
    std::vector<Edge2D> mEdgesNegativeOriented{};
    std::vector<Edge2D> mEdgesVertical{};

    /// Vertices container
    std::vector<Point2DType> mVerticesPositive{};
    std::vector<Point2DType> mVerticesNegative{};
    std::vector<Point2DType> mVerticesVertical{};

    /// Vertices sets. Used for fast search, if new vertex already exists.
    Unique<Point2DSetType> mVerticesSetPositive;
    Unique<Point2DSetType> mVerticesSetNegative;
    Unique<Point2DSetType> mVerticesSetVertical;

    /// Plane direction indices
    IndexType mPlaneIndex;
    IndexType DIRINDEX1; // In plane 1
    IndexType DIRINDEX2; // In plane 2
    IndexType DIRINDEX3; // Out of plane

    bool mUpperBoundary; // Is current plane upper bound?

    Point3DType mLowerBound; // Lower bound AABB
    Point3DType mUpperBound; // Upper bound AABB

    const TrimmedDomain* mpTrimmedDomain; // Reuiqred for IsInside()-Test.

    bool mSwitchOrientation; // Switch orientation of plane indices.

    double mSnapTolerance;
    ///@}
}; // End TrimmedDomainOnPlane

///@} // End QuESo Classes

} // End namespace queso
#endif // TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H