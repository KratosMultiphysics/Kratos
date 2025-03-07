namespace fcpw {

inline LineSegment::LineSegment()
{
    indices[0] = -1;
    indices[1] = -1;
    soup = nullptr;
    pIndex = -1;
}

inline BoundingBox<2> LineSegment::boundingBox() const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    BoundingBox<2> box(pa);
    box.expandToInclude(pb);

    return box;
}

inline Vector2 LineSegment::centroid() const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    return (pa + pb)*0.5f;
}

inline float LineSegment::surfaceArea() const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    return (pb - pa).norm();
}

inline float LineSegment::signedVolume() const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    return 0.5f*(pa[0]*pb[1] - pa[1]*pb[0]);
}

inline Vector2 LineSegment::normal(bool normalize) const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    Vector2 s = pb - pa;
    Vector2 n(s[1], -s[0]);

    return normalize ? n.normalized() : n;
}

inline Vector2 LineSegment::normal(int vIndex) const
{
    if (soup->vNormals.size() > 0 && vIndex >= 0) {
        return soup->vNormals[indices[vIndex]];
    }

    return normal(true);
}

inline Vector2 LineSegment::normal(const Vector1& uv) const
{
    int vIndex = -1;
    if (uv[0] <= epsilon) vIndex = 0;
    else if (uv[0] >= oneMinusEpsilon) vIndex = 1;

    return normal(vIndex);
}

inline Vector1 LineSegment::barycentricCoordinates(const Vector2& p) const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    return Vector1((p - pa).norm()/(pb - pa).norm());
}

inline float LineSegment::samplePoint(const Vector2& randNums, Vector1& uv, Vector2& p, Vector2& n) const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    Vector2 s = pb - pa;
    float area = s.norm();
    float u = randNums[1];
    uv[0] = u;
    p = pa + u*s;
    n[0] = s[1];
    n[1] = -s[0];
    n /= area;

    return 1.0f/area;
}

inline void LineSegment::split(int dim, float splitCoord, BoundingBox<2>& boxLeft,
                               BoundingBox<2>& boxRight) const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    if (pa[dim] <= splitCoord) {
        if (pb[dim] <= splitCoord) {
            boxLeft = BoundingBox<2>(pa);
            boxLeft.expandToInclude(pb);
            boxRight = BoundingBox<2>();

        } else {
            Vector2 u = pb - pa;
            float t = std::clamp((splitCoord - pa[dim])/u[dim], 0.0f, 1.0f);

            boxLeft = BoundingBox<2>(pa + u*t);
            boxRight = boxLeft;
            boxLeft.expandToInclude(pa);
            boxRight.expandToInclude(pb);
        }

    } else {
        if (pb[dim] >= splitCoord) {
            boxRight = BoundingBox<2>(pa);
            boxRight.expandToInclude(pb);
            boxLeft = BoundingBox<2>();

        } else {
            Vector2 u = pb - pa;
            float t = std::clamp((splitCoord - pa[dim])/u[dim], 0.0f, 1.0f);

            boxRight = BoundingBox<2>(pa + u*t);
            boxLeft = boxRight;
            boxRight.expandToInclude(pa);
            boxLeft.expandToInclude(pb);
        }
    }
}

inline bool LineSegment::intersect(const Ray<2>& r, Interaction<2>& i, bool checkForOcclusion) const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    Vector2 u = pa - r.o;
    Vector2 v = pb - pa;

    // return if line segment and ray are parallel
    float dv = r.d[0]*v[1] - r.d[1]*v[0];
    if (std::fabs(dv) <= epsilon) return 0;

    // solve r.o + t*r.d = pa + s*(pb - pa) for t >= 0 && 0 <= s <= 1
    // s = (u x r.d)/(r.d x v)
    float ud = u[0]*r.d[1] - u[1]*r.d[0];
    float s = ud/dv;

    if (s >= 0.0f && s <= 1.0f) {
        // t = (u x v)/(r.d x v)
        float uv = u[0]*v[1] - u[1]*v[0];
        float t = uv/dv;

        if (t >= 0.0f && t <= r.tMax) {
            if (checkForOcclusion) return true;
            i.d = t;
            i.p = pa + s*v;
            i.n = Vector2(v[1], -v[0]).normalized();
            i.uv[0] = s;
            i.primitiveIndex = pIndex;

            return true;
        }
    }

    return false;   
}

inline int LineSegment::intersect(const Ray<2>& r, std::vector<Interaction<2>>& is,
                                  bool checkForOcclusion, bool recordAllHits) const
{
    is.clear();
    Interaction<2> i;
    bool hit = intersect(r, i, checkForOcclusion);

    if (hit) {
        if (checkForOcclusion) return 1;

        is.emplace_back(i);
        return 1;
    }

    return 0;
}

template<size_t DIM>
inline float findClosestPointLineSegment(const Vector<DIM>& pa, const Vector<DIM>& pb,
                                         const Vector<DIM>& x, Vector<DIM>& pt, float& t)
{
    Vector<DIM> u = pb - pa;
    Vector<DIM> v = x - pa;

    float c1 = u.dot(v);
    if (c1 <= 0.0f) {
        pt = pa;
        t = 0.0f;

        return (x - pt).norm();
    }

    float c2 = u.dot(u);
    if (c2 <= c1) {
        pt = pb;
        t = 1.0f;

        return (x - pt).norm();
    }

    t = c1/c2;
    pt = pa + u*t;

    return (x - pt).norm();
}

inline bool LineSegment::intersect(const BoundingSphere<2>& s, Interaction<2>& i,
                                   bool recordSurfaceArea) const
{
    bool found = findClosestPoint(s, i);
    if (found) {
        i.d = recordSurfaceArea ? surfaceArea() : 1.0f;
        return true;
    }

    return false;
}

inline bool LineSegment::findClosestPoint(const BoundingSphere<2>& s, Interaction<2>& i) const
{
    const Vector2& pa = soup->positions[indices[0]];
    const Vector2& pb = soup->positions[indices[1]];

    float d = findClosestPointLineSegment<2>(pa, pb, s.c, i.p, i.uv[0]);
    if (d*d <= s.r2) {
        i.d = d;
        i.primitiveIndex = pIndex;

        return true;
    }

    return false;
}

} // namespace fcpw