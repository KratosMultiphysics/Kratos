//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef CLIPPER_INCLUDE_H
#define CLIPPER_INCLUDE_H

//// STL includes
#include <cstddef>
#include <array>
#include <algorithm>
//// Project includes
#include "queso/includes/define.hpp"
#include "queso/embedding/polygon.h"
#include "queso/containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Clipper
 * @author Manuel Messmer
 * @brief  Provides functions to clip triangle by cube.
 * @todo Rename to Triangle clipper?
*/
class Clipper  {

public:
    ///@name Type Definitions
    ///@{

    typedef Polygon<9> PolygonType; // Intersections contains maximum 9 vertices.

    enum Side{IN_FRONT_OF_PLANE, BEHIND_PLANE, ON_PLANE};

    ///@}
    ///@name Operations
    ///@{

    ///@brief Clips triangle by AABB. Expects that input triangle does intersect one of the faces of the bounding box.
    ///       If triangle is fully contained within AABB, same triangle is returned. If triangle is fully outside of AABB nullptr is returned.
    ///@param rV1 Vertex 1 of Triangle.
    ///@param rV2 Vertex 2 of Triangle.
    ///@param rV3 Vertex 3 of Triangle.
    ///@param rNormal Normal vector of triangle. Is provided, such that recomputation can be omitted. Recomputation would also lose accuracy.
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@return Unique<Polygon> (Will contain maximal 9 vertices).
    static Unique<PolygonType> ClipTriangle(const PointType& rV1, const PointType& rV2, const PointType& rV3,
                 const PointType& rNormal, const PointType& rLowerBound, const PointType& rUpperBound);

    ///@}

private:

    /**
     * @class  Plane
     * @author Manuel Messmer
     * @brief  Plane used in clipper.
     **/
    struct Plane {

        /// Constructor
        ///@param Index
        /// 0 -> -x axis;
        /// 1 -> +x axis;
        /// 2 -> -y axis;
        /// 3 -> +y axis;
        /// 4 -> -z axis;
        /// 5 -> +z axis;
        ///@param Position Offset from origin.
        Plane(IndexType Index, double Position) : mPlaneIndex(Index), mPosition(Position){}

        ///@brief Returns index of normal direction.
        IndexType GetIndexOfNormalDirection() const{
            return mPlaneIndex/2;
        }

        ///@brief Return trus if normal points in negative direction.
        bool IsNegativeOriented() const {
            return (mPlaneIndex % 2) == 0;
        }

        /// Member variables
        IndexType mPlaneIndex;
        double mPosition;
    };

    ///@brief Clips the polygon by plane. Keeps part that is in front of plane.
    ///@details This is a specialization of the Sutherland-Hodgeman clipping algorithm
    ///         for axis-aligned planes. See Section 8.3 of Christer Ericson's "Real-Time Collision Detection".
    ///         In addition, the information if a vertex is located on a plane is stored on the vertex.
    ///@param rPrevPoly Input polygon. Contains vertices to be clipped.
    ///@param rCurrentPoly Output polygon.
    ///@param rPlane Plane.
    static void ClipPolygonByPlane(const PolygonType* rPrevPoly,
                         PolygonType* rCurrentPoly,
                         const Plane& rPlane);

    ///@brief Classifies Point with respect to axis aligned plane.
    ///@param rPoint Point to be classified.
    ///@param rPlane Plane.
    ///@return IndexType (IN_FRON_OF_PLANE, BEHIND_PLANE, ON_PLANE)
    static IndexType ClassifyPointToPlane(const PointType& rPoint,
                               const Plane& rPlane,
                               const double Eps = 10.0*ZEROTOL);

    ///@brief Returns intersection point between line and axis aligned plane.
    ///@param rA Point behind the plane.
    ///@param rB Point in front of the plane.
    ///@return rPlane Clipping plane.
    static PointType FindIntersectionPointOnPlane(const PointType& rA,
                                                  const PointType& rB,
                                                  const Plane& rPlane);


}; // End Clipper class
///@} // End QuESo classes

} // End namespace queso

#endif