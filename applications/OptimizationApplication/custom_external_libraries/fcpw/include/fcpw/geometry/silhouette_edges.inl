namespace fcpw {

inline SilhouetteEdge::SilhouetteEdge()
{
    indices[0] = -1;
    indices[1] = -1;
    indices[2] = -1;
    indices[3] = -1;
    soup = nullptr;
    pIndex = -1;
}

inline BoundingBox<3> SilhouetteEdge::boundingBox() const
{
    const Vector3& pa = soup->positions[indices[1]];
    const Vector3& pb = soup->positions[indices[2]];

    BoundingBox<3> box(pa);
    box.expandToInclude(pb);

    return box;
}

inline Vector3 SilhouetteEdge::centroid() const
{
    const Vector3& pa = soup->positions[indices[1]];
    const Vector3& pb = soup->positions[indices[2]];

    return (pa + pb)*0.5f;
}

inline float SilhouetteEdge::surfaceArea() const
{
    const Vector3& pa = soup->positions[indices[1]];
    const Vector3& pb = soup->positions[indices[2]];

    return (pb - pa).norm();
}

inline bool SilhouetteEdge::hasFace(int fIndex) const
{
    return fIndex == 0 ? indices[3] != -1 : indices[0] != -1;
}

inline Vector3 SilhouetteEdge::normal(int fIndex, bool normalize) const
{
    int i, j, k;
    if (fIndex == 0) {
        i = 3;
        j = 1;
        k = 2;

    } else {
        i = 0;
        j = 2;
        k = 1;
    }

    const Vector3& pa = soup->positions[indices[j]];
    const Vector3& pb = soup->positions[indices[k]];
    const Vector3& pc = soup->positions[indices[i]];

    Vector3 v1 = pb - pa;
    Vector3 v2 = pc - pa;

    Vector3 n = v1.cross(v2);
    return normalize ? n.normalized() : n;
}

inline Vector3 SilhouetteEdge::normal() const
{
    if (soup->eNormals.size() > 0) {
        return soup->eNormals[pIndex];
    }

    Vector3 n = Vector3::Zero();
    if (hasFace(0)) n += normal(0, false);
    if (hasFace(1)) n += normal(1, false);

    return n.normalized();
}

inline bool isSilhouetteEdge(const Vector3& pa, const Vector3& pb,
                             const Vector3& n0, const Vector3& n1, const Vector3& viewDir,
                             float d, bool flipNormalOrientation, float precision)
{
    float sign = flipNormalOrientation ? 1.0f : -1.0f;

    // edge is a silhouette if it concave and the query point lies on the edge
    if (d <= precision) {
        Vector3 edgeDir = (pb - pa).normalized();
        float signedDihedralAngle = std::atan2(edgeDir.dot(n0.cross(n1)), n0.dot(n1));
        return sign*signedDihedralAngle > precision;
    }

    // edge is a silhouette if the query point lies on the halfplane defined
    // by an adjacent triangle and the other triangle is backfacing
    Vector3 viewDirUnit = viewDir/d;
    float dot0 = viewDirUnit.dot(n0);
    float dot1 = viewDirUnit.dot(n1);

    bool isZeroDot0 = std::fabs(dot0) <= precision;
    if (isZeroDot0) return sign*dot1 > precision;

    bool isZeroDot1 = std::fabs(dot1) <= precision;
    if (isZeroDot1) return sign*dot0 > precision;

    // edge is a silhouette if an adjacent triangle is frontfacing w.r.t. the
    // query point and the other triangle is backfacing
    return dot0*dot1 < 0.0f;
}

inline bool SilhouetteEdge::findClosestSilhouettePoint(const BoundingSphere<3>& s, Interaction<3>& i,
                                                       bool flipNormalOrientation, float squaredMinRadius,
                                                       float precision) const
{
    if (squaredMinRadius >= s.r2) return false;

    // compute closest edge point from query point
    const Vector3& pa = soup->positions[indices[1]];
    const Vector3& pb = soup->positions[indices[2]];
    float d = findClosestPointLineSegment<3>(pa, pb, s.c, i.p, i.uv[0]);
    if (d*d > s.r2) return false;

    // check if edge is a silhouette from view direction
    bool isSilhouette = !hasFace(0) || !hasFace(1);
    if (!isSilhouette) {
        Vector3 n0 = normal(0);
        Vector3 n1 = normal(1);
        Vector3 viewDir = s.c - i.p;
        isSilhouette = isSilhouetteEdge(pa, pb, n0, n1, viewDir, d, flipNormalOrientation, precision);
    }

    if (isSilhouette && d*d <= s.r2) {
        i.d = d;
        i.primitiveIndex = pIndex;
        i.uv[1] = -1;

        return true;
    }

    return false;
}

} // namespace fcpw