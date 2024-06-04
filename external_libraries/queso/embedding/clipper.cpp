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

//// Project includes
#include "queso/embedding/clipper.h"
#include "queso/utilities/math_utilities.hpp"
#include "queso/io/io_utilities.h"
#include <cstdlib>

namespace queso {

typedef Clipper::PolygonType PolygonType;

Unique<PolygonType> Clipper::ClipTriangle(const PointType& rV1, const PointType& rV2, const PointType& rV3,
            const PointType& rNormal, const PointType& rLowerBound, const PointType& rUpperBound ){

    Unique<PolygonType> p_current_poly = MakeUnique<PolygonType>(rNormal);
    Unique<PolygonType> p_prev_poly = MakeUnique<PolygonType>(rNormal);

    // Return nullptr, if no point is on bounded side.
    // Note on_bounded_side means behind plane, since normal vectors point in outward direction.
    for( IndexType plane_index = 0; plane_index < 6; ++plane_index){
        // Mapping of plane_index: [-x, x, -y, y, -z, z]
        const double plane_position = (plane_index%2UL) ? rUpperBound[plane_index/2UL] : rLowerBound[plane_index/2UL];
        Plane plane(plane_index, plane_position);
        auto loc_v1 = ClassifyPointToPlane(rV1, plane, 10.0*ZEROTOL);
        auto loc_v2 = ClassifyPointToPlane(rV2, plane, 10.0*ZEROTOL);
        auto loc_v3 = ClassifyPointToPlane(rV3, plane, 10.0*ZEROTOL);
        IndexType count_on_plane = (loc_v1 == ON_PLANE) + (loc_v2 == ON_PLANE) + (loc_v3 == ON_PLANE);
        if( count_on_plane == 3 ){
            return nullptr; // Triangle is aligned with plane
        }
    }

    // Compute min/max values of bounding box around triangle.
    const PointType x_values{rV1[0], rV2[0], rV3[0]};
    const PointType y_values{rV1[1], rV2[1], rV3[1]};
    const PointType z_values{rV1[2], rV2[2], rV3[2]};

    auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
    auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
    auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());

    PointType min_tri{*x_min_max.first, *y_min_max.first, *z_min_max.first};
    PointType max_tri{*x_min_max.second, *y_min_max.second, *z_min_max.second};

    // Add triangle vertices to polygon.
    p_current_poly->AddVertex( rV1 );
    p_current_poly->AddVertex( rV2 );
    p_current_poly->AddVertex( rV3 );

    /// Loop over all three dimensions, and clip each plane (positive and negative direction.)
    for(IndexType dimension = 0; dimension < 3; ++dimension)
    {
        if(max_tri[dimension] > rLowerBound[dimension] )
        {
            Ptr::swap(p_prev_poly, p_current_poly);
            // Plane is oriented in negative direction.
            Plane plane(2 * dimension + 0, rLowerBound[dimension]);
            ClipPolygonByPlane(p_prev_poly.get(), p_current_poly.get(), plane);
        }
        if(min_tri[dimension] < rUpperBound[dimension])
        {
            Ptr::swap(p_prev_poly, p_current_poly);
            // Plane is positive in negative direction.
            Plane plane(2 * dimension + 1, rUpperBound[dimension]);
            ClipPolygonByPlane(p_prev_poly.get(), p_current_poly.get(), plane);
        }

    }
    return p_current_poly;

}

void Clipper::ClipPolygonByPlane(const PolygonType* pPrevPoly,
                         PolygonType* pCurrentPoly,
                         const Plane& rPlane) {

    pCurrentPoly->Clear();
    IndexType numVerts = pPrevPoly->NumVertices();

    if(numVerts == 0) {
        return;
    }

    // Initialize point a with the last vertex of the polygon
    const auto* a = &pPrevPoly->GetLastVertex();
    // Classify point
    IndexType aSide = ClassifyPointToPlane((*a).first, rPlane);

    /// Loop over all vertices.
    for(IndexType i = 0; i < numVerts; ++i) {
        // Classify next point.
        const auto& b = pPrevPoly->GetVertex(i);
        IndexType bSide = ClassifyPointToPlane(b.first, rPlane);

        // Switch regarding the relative location of a and b to current plane.
        switch(bSide)
        {
        case IN_FRONT_OF_PLANE: {
            if(aSide == BEHIND_PLANE) {
                // Set new point on plane, if both points a and b are on plane.
                // Only up to n=rPlane.mPlaneIndex are required, since others have not been checked yet.
                std::array<bool,6> p_on_plane = {false}; // Initalizes all elements to 'false'.
                std::generate_n(p_on_plane.begin(), rPlane.mPlaneIndex,
                    [&a, &b, j=-1]()mutable->bool {++j; return a->second[j] && b.second[j]; });

                IndexType index = pCurrentPoly->AddVertex(FindIntersectionPointOnPlane( (*a).first, b.first, rPlane), p_on_plane);
                pCurrentPoly->SetIndexOnPlane(index, rPlane.mPlaneIndex);
            }
            break;
        }
        case ON_PLANE: {
            if(aSide == BEHIND_PLANE) {
                IndexType index = pCurrentPoly->AddVertex(b);
                pCurrentPoly->SetIndexOnPlane(index, rPlane.mPlaneIndex);
            }

            break;
        }
        case BEHIND_PLANE: {
            switch(aSide)
            {
            case IN_FRONT_OF_PLANE: {
                // Set new point on plane, if both points a and b are on plane.
                // Only up to n=rPlane.mPlaneIndex are required, since others have not been checked yet.
                std::array<bool,6> p_on_plane = {false}; // Initalizes all elements to 'false'.
                std::generate_n(p_on_plane.begin(), rPlane.mPlaneIndex,
                    [&a, &b, j=-1]()mutable->bool {++j; return a->second[j] && b.second[j]; });

                IndexType index = pCurrentPoly->AddVertex(FindIntersectionPointOnPlane( (*a).first, b.first, rPlane), p_on_plane);
                pCurrentPoly->SetIndexOnPlane(index, rPlane.mPlaneIndex);
                pCurrentPoly->AddVertex(b); // This new point is not neccessarily on plane, but is classified in next iteration.
                break;
            }
            case ON_PLANE: {
                IndexType index = pCurrentPoly->AddVertex(*a);
                pCurrentPoly->SetIndexOnPlane(index, rPlane.mPlaneIndex);
                pCurrentPoly->AddVertex(b); // This new point is not neccessarily on plane, but is classified in next iteration.

                break;
            }
            case BEHIND_PLANE: {
                pCurrentPoly->AddVertex(b); // This new point is not neccessarily on plane, but is classified in next iteration.
                break;
            }
            }
            break;
        }
        }
        // swap a and b
        a = &b;
        aSide = bSide;
    }
}

IndexType Clipper::ClassifyPointToPlane(const PointType& rPoint,
                            const Plane& rPlane,
                            const double Eps) {

    // Note: Plane is axis aligned. Distance is directly given by coordinate.
    const bool is_negative_oriented = rPlane.IsNegativeOriented();
    const IndexType index = rPlane.GetIndexOfNormalDirection();

    const double dist = is_negative_oriented ? rPlane.mPosition - rPoint[index] : rPoint[index] - rPlane.mPosition;

    if(dist > Eps) {
        return IN_FRONT_OF_PLANE;
    }
    if(dist < -Eps) {
        return BEHIND_PLANE;
    }

    return ON_PLANE;
}

PointType Clipper::FindIntersectionPointOnPlane(const PointType& rA,
                                                const PointType& rB,
                                                const Plane& rPlane) {
    // Need to find a parameter t for the point pt, such that,
    // * 0 <= t <= 1
    // * pt = rA + t * (rB-rA)
    // * pt[index]  == val
    IndexType index = rPlane.GetIndexOfNormalDirection();
    double t = (rPlane.mPosition - rA[index]) / (rB[index] - rA[index]);

    if( !(0.0 <= t && t <= 1.0) ){
        QuESo_ERROR << "No valid intersection.\n";
    }

    PointType ret{ rA[0] + t*(rB[0] - rA[0]),
                   rA[1] + t*(rB[1] - rA[1]),
                   rA[2] + t*(rB[2] - rA[2]) };

    auto status = ClassifyPointToPlane(ret, rPlane, ZEROTOL);
    if( status != ON_PLANE){
        ret[index] = rPlane.mPosition;
    }
    //QuESo_ERROR_IF("Clipper::FindIntersectionPointOnPlane", status != ON_PLANE) << "Intersection Point is not on Boundary\n.";

    return ret;
}

} // End namespace queso
