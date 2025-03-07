#pragma once

#include <fcpw/core/bounding_volumes.h>

namespace fcpw {

// performs wide version of ray box intersection test
template<size_t WIDTH, size_t DIM>
inline MaskP<WIDTH> intersectWideBox(const VectorP<WIDTH, DIM>& bMin, const VectorP<WIDTH, DIM>& bMax,
                                     const enokiVector<DIM>& ro, const enokiVector<DIM>& rinvD, float rtMax,
                                     FloatP<WIDTH>& tMin, FloatP<WIDTH>& tMax)
{
    // vectorized slab test
    VectorP<WIDTH, DIM> t0 = (bMin - ro)*rinvD;
    VectorP<WIDTH, DIM> t1 = (bMax - ro)*rinvD;
    VectorP<WIDTH, DIM> tNear = enoki::min(t0, t1);
    VectorP<WIDTH, DIM> tFar = enoki::max(t0, t1);

    tMin = enoki::max(0.0f, enoki::hmax(tNear));
    tMax = enoki::min(rtMax, enoki::hmin(tFar));

    return tMin <= tMax;
}

// performs wide version of ray line segment intersection test
template<size_t WIDTH>
inline MaskP<WIDTH> intersectWideLineSegment(const Vector2P<WIDTH>& pa, const Vector2P<WIDTH>& pb,
                                             const enokiVector2& ro, const enokiVector2& rd, float rtMax,
                                             FloatP<WIDTH>& d, Vector2P<WIDTH>& pt, Vector2P<WIDTH>& n,
                                             FloatP<WIDTH>& t, bool checkForOcclusion)
{
    Vector2P<WIDTH> u = pa - ro;
    Vector2P<WIDTH> v = pb - pa;

    // track non-parallel line segments and rays
    FloatP<WIDTH> dv = rd[0]*v[1] - rd[1]*v[0];
    MaskP<WIDTH> active = enoki::abs(dv) > epsilon;
    FloatP<WIDTH> invDv = enoki::rcp(dv);

    // solve r.o + s*r.d = pa + t*(pb - pa) for s >= 0 && 0 <= t <= 1
    // t = (u x r.d)/(r.d x v)
    FloatP<WIDTH> ud = u[0]*rd[1] - u[1]*rd[0];
    t = ud*invDv;
    active &= t >= 0.0f && t <= 1.0f;

    // s = (u x v)/(r.d x v)
    FloatP<WIDTH> uv = u[0]*v[1] - u[1]*v[0];
    d = uv*invDv;
    active &= d >= 0.0f && d <= rtMax;
    if (checkForOcclusion) return active;
    pt = pa + t*v;
    n[0] = v[1];
    n[1] = -v[0];
    n = enoki::normalize(n);

    return active;
}

// performs wide version of ray triangle intersection test
template<size_t WIDTH>
inline MaskP<WIDTH> intersectWideTriangle(const Vector3P<WIDTH>& pa, const Vector3P<WIDTH>& pb, const Vector3P<WIDTH>& pc,
                                          const enokiVector3& ro, const enokiVector3& rd, float rtMax, FloatP<WIDTH>& d,
                                          Vector3P<WIDTH>& pt, Vector3P<WIDTH>& n, Vector2P<WIDTH>& t, bool checkForOcclusion)
{
    // vectorized Möller–Trumbore intersection algorithm
    Vector3P<WIDTH> v1 = pb - pa;
    Vector3P<WIDTH> v2 = pc - pa;
    Vector3P<WIDTH> p = enoki::cross(rd, v2);
    FloatP<WIDTH> det = enoki::dot(v1, p);

    MaskP<WIDTH> active = enoki::abs(det) > epsilon;
    FloatP<WIDTH> invDet = enoki::rcp(det);

    Vector3P<WIDTH> s = ro - pa;
    FloatP<WIDTH> v = enoki::dot(s, p)*invDet;
    active &= v >= 0.0f && v <= 1.0f;

    Vector3P<WIDTH> q = enoki::cross(s, v1);
    FloatP<WIDTH> w = enoki::dot(rd, q)*invDet;
    active &= w >= 0.0f && v + w <= 1.0f;

    d = enoki::dot(v2, q)*invDet;
    active &= d >= 0.0f && d <= rtMax;
    if (checkForOcclusion) return active;
    pt = pa + v1*v + v2*w;
    n = enoki::normalize(enoki::cross(v1, v2));
    t[0] = 1.0f - v - w;
    t[1] = v;

    return active;
}

// performs wide version of sphere box overlap test
template<size_t WIDTH, size_t DIM>
inline MaskP<WIDTH> overlapWideBox(const VectorP<WIDTH, DIM>& bMin, const VectorP<WIDTH, DIM>& bMax,
                                   const enokiVector<DIM>& sc, float sr2,
                                   FloatP<WIDTH>& d2Min, FloatP<WIDTH>& d2Max)
{
    VectorP<WIDTH, DIM> u = bMin - sc;
    VectorP<WIDTH, DIM> v = sc - bMax;
    d2Min = enoki::squared_norm(enoki::max(enoki::max(u, v), 0.0f));
    d2Max = enoki::squared_norm(enoki::min(u, v));

    return d2Min <= sr2;
}

// performs wide version of sphere box overlap test
template<size_t WIDTH, size_t DIM>
inline MaskP<WIDTH> overlapWideBox(const VectorP<WIDTH, DIM>& bMin, const VectorP<WIDTH, DIM>& bMax,
                                   const enokiVector<DIM>& sc, float sr2, FloatP<WIDTH>& d2Min)
{
    VectorP<WIDTH, DIM> u = bMin - sc;
    VectorP<WIDTH, DIM> v = sc - bMax;
    d2Min = enoki::squared_norm(enoki::max(enoki::max(u, v), 0.0f));

    return d2Min <= sr2;
}

// finds closest point on wide line segment to point
template<size_t WIDTH, size_t DIM>
inline FloatP<WIDTH> findClosestPointWideLineSegment(const VectorP<WIDTH, DIM>& pa, const VectorP<WIDTH, DIM>& pb,
                                                     const enokiVector<DIM>& x, VectorP<WIDTH, DIM>& pt, FloatP<WIDTH>& t)
{
    VectorP<WIDTH, DIM> u = pb - pa;
    VectorP<WIDTH, DIM> v = x - pa;

    // project x onto u
    FloatP<WIDTH> c1 = enoki::dot(u, v);
    FloatP<WIDTH> c2 = enoki::dot(u, u);
    MaskP<WIDTH> active1 = c1 <= 0.0f;
    MaskP<WIDTH> active2 = c2 <= c1;

    // compute closest point
    t = c1*enoki::rcp(c2);
    enoki::masked(t, active1) = 0.0f;
    enoki::masked(t, active2) = 1.0f;
    pt = pa + u*t;

    return enoki::norm(x - pt);
}

// finds closest point on wide triangle to point
template<size_t WIDTH>
inline FloatP<WIDTH> findClosestPointWideTriangle(const Vector3P<WIDTH>& pa, const Vector3P<WIDTH>& pb, const Vector3P<WIDTH>& pc,
                                                  const enokiVector3& x, Vector3P<WIDTH>& pt, Vector2P<WIDTH>& t)
{
    // check if x in vertex region outside pa
    Vector3P<WIDTH> ab = pb - pa;
    Vector3P<WIDTH> ac = pc - pa;
    Vector3P<WIDTH> ax = x - pa;
    FloatP<WIDTH> d1 = enoki::dot(ab, ax);
    FloatP<WIDTH> d2 = enoki::dot(ac, ax);
    MaskP<WIDTH> active1 = d1 <= 0.0f && d2 <= 0.0f;
    MaskP<WIDTH> active7 = active1;

    // barycentric coordinates (1, 0, 0)
    enoki::masked(pt, active1) = pa;
    enoki::masked(t[0], active1) = 1.0f;
    enoki::masked(t[1], active1) = 0.0f;
    if (enoki::all(active7)) return enoki::norm(x - pt);

    // check if x in vertex region outside pb
    Vector3P<WIDTH> bx = x - pb;
    FloatP<WIDTH> d3 = enoki::dot(ab, bx);
    FloatP<WIDTH> d4 = enoki::dot(ac, bx);
    MaskP<WIDTH> active2 = d3 >= 0.0f && d4 <= d3;
    active7 |= active2;

    // barycentric coordinates (0, 1, 0)
    enoki::masked(pt, active2) = pb;
    enoki::masked(t[0], active2) = 0.0f;
    enoki::masked(t[1], active2) = 1.0f;
    if (enoki::all(active7)) return enoki::norm(x - pt);

    // check if x in vertex region outside pc
    Vector3P<WIDTH> cx = x - pc;
    FloatP<WIDTH> d5 = enoki::dot(ab, cx);
    FloatP<WIDTH> d6 = enoki::dot(ac, cx);
    MaskP<WIDTH> active3 = d6 >= 0.0f && d5 <= d6;
    active7 |= active3;

    // barycentric coordinates (0, 0, 1)
    enoki::masked(pt, active3) = pc;
    enoki::masked(t[0], active3) = 0.0f;
    enoki::masked(t[1], active3) = 0.0f;
    if (enoki::all(active7)) return enoki::norm(x - pt);

    // check if x in edge region of ab, if so return projection of x onto ab
    FloatP<WIDTH> vc = d1*d4 - d3*d2;
    MaskP<WIDTH> active4 = vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f;
    active7 |= active4;

    // barycentric coordinates (1 - v, v, 0)
    FloatP<WIDTH> v = d1*enoki::rcp(d1 - d3);
    enoki::masked(pt, active4) = pa + ab*v;
    enoki::masked(t[0], active4) = 1.0f - v;
    enoki::masked(t[1], active4) = v;
    if (enoki::all(active7)) return enoki::norm(x - pt);

    // check if x in edge region of ac, if so return projection of x onto ac
    FloatP<WIDTH> vb = d5*d2 - d1*d6;
    MaskP<WIDTH> active5 = vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f;
    active7 |= active5;

    // barycentric coordinates (1 - w, 0, w)
    FloatP<WIDTH> w = d2*enoki::rcp(d2 - d6);
    enoki::masked(pt, active5) = pa + ac*w;
    enoki::masked(t[0], active5) = 1.0f - w;
    enoki::masked(t[1], active5) = 0.0f;
    if (enoki::all(active7)) return enoki::norm(x - pt);

    // check if x in edge region of bc, if so return projection of x onto bc
    FloatP<WIDTH> va = d3*d6 - d5*d4;
    MaskP<WIDTH> active6 = va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f;
    active7 |= active6;

    // barycentric coordinates (0, 1 - w, w)
    w = (d4 - d3)*enoki::rcp((d4 - d3) + (d5 - d6));
    enoki::masked(pt, active6) = pb + (pc - pb)*w;
    enoki::masked(t[0], active6) = 0.0f;
    enoki::masked(t[1], active6) = 1.0f - w;
    if (enoki::all(active7)) return enoki::norm(x - pt);

    // x inside face region. Compute pt through its barycentric coordinates (u, v, w)
    FloatP<WIDTH> denom = enoki::rcp(va + vb + vc);
    v = vb*denom;
    w = vc*denom;
    active7 = ~active7;

    enoki::masked(pt, active7) = pa + ab*v + ac*w; //= u*a + v*b + w*c, u = va*denom = 1.0f - v - w
    enoki::masked(t[0], active7) = 1.0f - v - w;
    enoki::masked(t[1], active7) = v;
    return enoki::norm(x - pt);
}

// computes orthonormal basis for direction in 3d
template<size_t WIDTH>
inline void computeWideOrthonormalBasis(const Vector3P<WIDTH>& n, Vector3P<WIDTH>& b1, Vector3P<WIDTH>& b2)
{
    // source: https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    FloatP<WIDTH> sign = enoki::copysign(1.0f, n[2]);
    const FloatP<WIDTH> a = -1.0f*enoki::rcp(sign + n[2]);
    const FloatP<WIDTH> b = n[0]*n[1]*a;

    b1[0] = 1.0f + sign*n[0]*n[0]*a;
    b1[1] = sign*b;
    b1[2] = -sign*n[0];
    b2[0] = b;
    b2[1] = sign + n[1]*n[1]*a;
    b2[2] = -n[1];
}

// projects vector onto plane defined by n
template<size_t WIDTH, size_t DIM>
inline FloatP<WIDTH> projectToWidePlane(const VectorP<WIDTH, DIM>& n, const VectorP<WIDTH, DIM>& e)
{
    std::cerr << "projectToWidePlane(): WIDTH: " << WIDTH << " DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return 0.0f;
}

// projects vector onto plane defined by n in 2d
template<size_t WIDTH>
inline FloatP<WIDTH> projectToWidePlane(const Vector2P<WIDTH>& n, const Vector2P<WIDTH>& e)
{
    // compute orthonormal basis
    Vector2P<WIDTH> b;
    b[0] = -n[1];
    b[1] = n[0];

    // compute maximal projection radius
    FloatP<WIDTH> r = enoki::dot(e, enoki::abs(b));
    return enoki::abs(r);
}

// projects vector onto plane defined by n in 3d
template<size_t WIDTH>
inline FloatP<WIDTH> projectToWidePlane(const Vector3P<WIDTH>& n, const Vector3P<WIDTH>& e)
{
    // compute orthonormal basis
    Vector3P<WIDTH> b1, b2;
    computeWideOrthonormalBasis<WIDTH>(n, b1, b2);

    // compute maximal projection radius
    FloatP<WIDTH> r1 = enoki::dot(e, enoki::abs(b1));
    FloatP<WIDTH> r2 = enoki::dot(e, enoki::abs(b2));
    return enoki::sqrt(r1*r1 + r2*r2);
}

// performs wide version of cone cone overlap test
template<size_t WIDTH, size_t DIM>
inline void overlapWideCone(const VectorP<WIDTH, DIM>& normalConeAxis, const FloatP<WIDTH>& normalConeHalfAngle,
                            const FloatP<WIDTH>& normalConeRadius, const enokiVector<DIM>& co,
                            const VectorP<WIDTH, DIM>& bMin, const VectorP<WIDTH, DIM>& bMax,
                            const FloatP<WIDTH>& distToBox, FloatP<WIDTH>& minAngleRange,
                            FloatP<WIDTH>& maxAngleRange, MaskP<WIDTH>& overlap)
{
    minAngleRange = 0.0f;
    maxAngleRange = M_PI_2;
    MaskP<WIDTH> activeQuery = overlap && normalConeHalfAngle < M_PI_2 && distToBox > epsilon;

    if (enoki::any(activeQuery)) {
        // compute the view cone axis
        VectorP<WIDTH, DIM> bc = 0.5f*(bMin + bMax);
        VectorP<WIDTH, DIM> viewConeAxis = bc - co;
        FloatP<WIDTH> l = enoki::norm(viewConeAxis);
        FloatP<WIDTH> lInv = enoki::rcp(l);
        viewConeAxis *= lInv;

        // compute the angle between the view and normal cone axes
        FloatP<WIDTH> dAxisAngle = enoki::acos(enoki::max(-1.0f, enoki::min(1.0f, enoki::dot(normalConeAxis, viewConeAxis)))); // [0, 180]

        // check if the view cone origin lies outside the normal cone's bounding sphere;
        // if it does, compute the view cone halfAngle and check for overlap
        MaskP<WIDTH> outsideSphere = l > normalConeRadius;
        MaskP<WIDTH> activeQueryOutsideSphere = activeQuery && outsideSphere;
        {
            FloatP<WIDTH> viewConeHalfAngle = enoki::asin(normalConeRadius*lInv);
            FloatP<WIDTH> halfAngleSum = normalConeHalfAngle + viewConeHalfAngle;
            minAngleRange = dAxisAngle - halfAngleSum;
            maxAngleRange = dAxisAngle + halfAngleSum;
            MaskP<WIDTH> angleInRange = M_PI_2 >= minAngleRange && M_PI_2 <= maxAngleRange;
            enoki::masked(overlap, activeQueryOutsideSphere) = halfAngleSum >= M_PI_2 || angleInRange;
        }
        if (enoki::all(activeQuery == activeQueryOutsideSphere)) return; // try to early out since the next two cases are less common

        // the view cone origin lies inside the box's bounding sphere, so check if
        // the plane defined by the view cone axis intersects the box; if it does, then
        // there's overlap since the view cone has a halfAngle greater than 90 degrees
        VectorP<WIDTH, DIM> be = bMax - bc;
        FloatP<WIDTH> d = enoki::dot(be, enoki::abs(viewConeAxis));
        FloatP<WIDTH> s = l - d;
        MaskP<WIDTH> intersectPlane = s <= 0.0f;
        MaskP<WIDTH> activeQueryInsideSphere = activeQuery && ~outsideSphere;

        // compute the view cone halfAngle by projecting the max extents of the box
        // onto the plane, and check for overlap
        d = projectToWidePlane(viewConeAxis, be);
        {
            FloatP<WIDTH> viewConeHalfAngle = enoki::atan2(d, s);
            FloatP<WIDTH> halfAngleSum = normalConeHalfAngle + viewConeHalfAngle;
            minAngleRange = dAxisAngle - halfAngleSum;
            maxAngleRange = dAxisAngle + halfAngleSum;
            MaskP<WIDTH> angleInRange = M_PI_2 >= minAngleRange && M_PI_2 <= maxAngleRange;
            enoki::masked(overlap, activeQueryInsideSphere && ~intersectPlane) = halfAngleSum >= M_PI_2 || angleInRange;
        }
    }
}

// checks whether vertex is a silhouette
template<size_t WIDTH>
inline MaskP<WIDTH> isWideSilhouetteVertex(const Vector2P<WIDTH>& n0, const Vector2P<WIDTH>& n1,
                                           const Vector2P<WIDTH>& viewDir, const FloatP<WIDTH>& d,
                                           bool flipNormalOrientation, float precision)
{
    float sign = flipNormalOrientation ? 1.0f : -1.0f;

    // vertex is a silhouette point if it concave and the query point lies on the vertex
    FloatP<WIDTH> det = n0[0]*n1[1] - n0[1]*n1[0];
    MaskP<WIDTH> isSilhouette = false;
    MaskP<WIDTH> active = d > precision;
    enoki::masked(isSilhouette, ~active) = sign*det > precision;

    // vertex is a silhouette point if the query point lies on the halfplane
    // defined by an adjacent line segment and the other segment is backfacing
    FloatP<WIDTH> invD = enoki::rcp(d);
    Vector2P<WIDTH> viewDirUnit = viewDir*invD;
    FloatP<WIDTH> dot0 = enoki::dot(viewDirUnit, n0);
    FloatP<WIDTH> dot1 = enoki::dot(viewDirUnit, n1);
    MaskP<WIDTH> isZeroDot0 = enoki::abs(dot0) <= precision;
    enoki::masked(isSilhouette, active && isZeroDot0) = sign*dot1 > precision;

    active &= ~isZeroDot0;
    MaskP<WIDTH> isZeroDot1 = enoki::abs(dot1) <= precision;
    enoki::masked(isSilhouette, active && isZeroDot1) = sign*dot0 > precision;

    // vertex is a silhouette point if an adjacent line segment is frontfacing
    // w.r.t. the query point and the other segment is backfacing
    active &= ~isZeroDot1;
    enoki::masked(isSilhouette, active) = dot0*dot1 < 0.0f;

    return isSilhouette;
}

// checks whether edge is a silhouette
template<size_t WIDTH>
inline MaskP<WIDTH> isWideSilhouetteEdge(const Vector3P<WIDTH>& pb, const Vector3P<WIDTH>& pc,
                                         const Vector3P<WIDTH>& n0, const Vector3P<WIDTH>& n1,
                                         const Vector3P<WIDTH>& viewDir, const FloatP<WIDTH>& d,
                                         bool flipNormalOrientation, float precision)
{
    float sign = flipNormalOrientation ? 1.0f : -1.0f;

    // edge is a silhouette if it concave and the query point lies on the edge
    Vector3P<WIDTH> edgeDir = enoki::normalize(pc - pb);
    FloatP<WIDTH> signedDihedralAngle = enoki::atan2(enoki::dot(edgeDir, enoki::cross(n0, n1)), enoki::dot(n0, n1));
    MaskP<WIDTH> isSilhouette = false;
    MaskP<WIDTH> active = d > precision;
    enoki::masked(isSilhouette, ~active) = sign*signedDihedralAngle > precision;

    // edge is a silhouette if the query point lies on the halfplane defined
    // by an adjacent triangle and the other triangle is backfacing
    FloatP<WIDTH> invD = enoki::rcp(d);
    Vector3P<WIDTH> viewDirUnit = viewDir*invD;
    FloatP<WIDTH> dot0 = enoki::dot(viewDirUnit, n0);
    FloatP<WIDTH> dot1 = enoki::dot(viewDirUnit, n1);
    MaskP<WIDTH> isZeroDot0 = enoki::abs(dot0) <= precision;
    enoki::masked(isSilhouette, active && isZeroDot0) = sign*dot1 > precision;

    active &= ~isZeroDot0;
    MaskP<WIDTH> isZeroDot1 = enoki::abs(dot1) <= precision;
    enoki::masked(isSilhouette, active && isZeroDot1) = sign*dot0 > precision;

    // edge is a silhouette if an adjacent triangle is frontfacing w.r.t. the
    // query point and the other triangle is backfacing
    active &= ~isZeroDot1;
    enoki::masked(isSilhouette, active) = dot0*dot1 < 0.0f;

    return isSilhouette;
}

} // namespace fcpw