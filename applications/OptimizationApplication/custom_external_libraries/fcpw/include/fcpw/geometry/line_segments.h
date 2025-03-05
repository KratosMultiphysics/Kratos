#pragma once

#include <fcpw/geometry/polygon_soup.h>

namespace fcpw {

class LineSegment: public GeometricPrimitive<2> {
public:
    // constructor
    LineSegment();

    // returns bounding box
    BoundingBox<2> boundingBox() const;

    // returns centroid
    Vector2 centroid() const;

    // returns surface area
    float surfaceArea() const;

    // returns signed volume
    float signedVolume() const;

    // returns normal
    Vector2 normal(bool normalize=false) const;

    // returns the normalized normal based on the local parameterization
    Vector2 normal(const Vector1& uv) const;

    // returns barycentric coordinates
    Vector1 barycentricCoordinates(const Vector2& p) const;

    // samples a random point on the geometric primitive and returns sampling pdf
    float samplePoint(const Vector2& randNums, Vector1& uv, Vector2& p, Vector2& n) const;

    // splits the line segment along the provided coordinate and axis
    void split(int dim, float splitCoord, BoundingBox<2>& boxLeft, BoundingBox<2>& boxRight) const;

    // intersects with ray
    bool intersect(const Ray<2>& r, Interaction<2>& i, bool checkForOcclusion=false) const;

    // intersects with ray
    int intersect(const Ray<2>& r, std::vector<Interaction<2>>& is,
                  bool checkForOcclusion=false, bool recordAllHits=false) const;

    // intersects with sphere
    bool intersect(const BoundingSphere<2>& s, Interaction<2>& i,
                   bool recordSurfaceArea=false) const;

    // finds closest point to sphere center
    bool findClosestPoint(const BoundingSphere<2>& s, Interaction<2>& i) const;

    // get and set index
    int getIndex() const { return pIndex; }
    void setIndex(int index) { pIndex = index; }

    // members
    int indices[2];
    int pIndex;
    const PolygonSoup<2> *soup;

protected:
    // returns normalized vertex normal if available;
    // otherwise computes normalized segment normal
    Vector2 normal(int vIndex) const;
};

} // namespace fcpw

#include "line_segments.inl"