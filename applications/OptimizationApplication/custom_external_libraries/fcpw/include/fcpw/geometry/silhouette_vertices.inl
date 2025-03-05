namespace fcpw {

inline SilhouetteVertex::SilhouetteVertex()
{
    indices[0] = -1;
    indices[1] = -1;
    indices[2] = -1;
    soup = nullptr;
    pIndex = -1;
}

inline BoundingBox<2> SilhouetteVertex::boundingBox() const
{
    const Vector2& p = soup->positions[indices[1]];
    return BoundingBox<2>(p);
}

inline Vector2 SilhouetteVertex::centroid() const
{
    const Vector2& p = soup->positions[indices[1]];
    return p;
}

inline float SilhouetteVertex::surfaceArea() const
{
    return 0.0f;
}

inline bool SilhouetteVertex::hasFace(int fIndex) const
{
    return fIndex == 0 ? indices[2] != -1 : indices[0] != -1;
}

inline Vector2 SilhouetteVertex::normal(int fIndex, bool normalize) const
{
    int i = fIndex == 0 ? 1 : 0;
    const Vector2& pa = soup->positions[indices[i + 0]];
    const Vector2& pb = soup->positions[indices[i + 1]];

    Vector2 s = pb - pa;
    Vector2 n(s[1], -s[0]);

    return normalize ? n.normalized() : n;
}

inline Vector2 SilhouetteVertex::normal() const
{
    if (soup->vNormals.size() > 0) {
        return soup->vNormals[indices[1]];
    }

    Vector2 n = Vector2::Zero();
    if (hasFace(0)) n += normal(0, false);
    if (hasFace(1)) n += normal(1, false);

    return n.normalized();
}

inline bool isSilhouetteVertex(const Vector2& n0, const Vector2& n1, const Vector2& viewDir,
                               float d, bool flipNormalOrientation, float precision)
{
    float sign = flipNormalOrientation ? 1.0f : -1.0f;

    // vertex is a silhouette point if it concave and the query point lies on the vertex
    if (d <= precision) {
        float det = n0[0]*n1[1] - n0[1]*n1[0];
        return sign*det > precision;
    }

    // vertex is a silhouette point if the query point lies on the halfplane
    // defined by an adjacent line segment and the other segment is backfacing
    Vector2 viewDirUnit = viewDir/d;
    float dot0 = viewDirUnit.dot(n0);
    float dot1 = viewDirUnit.dot(n1);

    bool isZeroDot0 = std::fabs(dot0) <= precision;
    if (isZeroDot0) return sign*dot1 > precision;

    bool isZeroDot1 = std::fabs(dot1) <= precision;
    if (isZeroDot1) return sign*dot0 > precision;

    // vertex is a silhouette point if an adjacent line segment is frontfacing
    // w.r.t. the query point and the other segment is backfacing
    return dot0*dot1 < 0.0f;
}

inline bool SilhouetteVertex::findClosestSilhouettePoint(const BoundingSphere<2>& s, Interaction<2>& i,
                                                         bool flipNormalOrientation, float squaredMinRadius,
                                                         float precision) const
{
    if (squaredMinRadius >= s.r2) return false;

    // compute view direction
    const Vector2& p = soup->positions[indices[1]];
    Vector2 viewDir = s.c - p;
    float d = viewDir.norm();
    if (d*d > s.r2) return false;

    // check if vertex is a silhouette point from view direction
    bool isSilhouette = !hasFace(0) || !hasFace(1);
    if (!isSilhouette) {
        Vector2 n0 = normal(0);
        Vector2 n1 = normal(1);
        isSilhouette = isSilhouetteVertex(n0, n1, viewDir, d, flipNormalOrientation, precision);
    }

    if (isSilhouette && d*d <= s.r2) {
        i.d = d;
        i.p = p;
        i.primitiveIndex = pIndex;
        i.uv[0] = -1;

        return true;
    }

    return false;
}

} // namespace fcpw