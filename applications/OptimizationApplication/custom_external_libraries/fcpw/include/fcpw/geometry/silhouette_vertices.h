#pragma once

#include <fcpw/geometry/polygon_soup.h>

namespace fcpw {

class SilhouetteVertex: public SilhouettePrimitive<2> {
public:
    // constructor
    SilhouetteVertex();

    // returns bounding box
    BoundingBox<2> boundingBox() const;

    // returns centroid
    Vector2 centroid() const;

    // returns surface area
    float surfaceArea() const;

    // checks whether silhouette has adjacent face
    bool hasFace(int fIndex) const;

    // returns normal of adjacent face
    Vector2 normal(int fIndex, bool normalize=true) const;

    // returns normalized silhouette normal
    Vector2 normal() const;

    // finds closest silhouette point to sphere center
    bool findClosestSilhouettePoint(const BoundingSphere<2>& s, Interaction<2>& i,
                                    bool flipNormalOrientation=false,
                                    float squaredMinRadius=0.0f,
                                    float precision=1e-3f) const;

    // get and set index
    int getIndex() const { return pIndex; }
    void setIndex(int index) { pIndex = index; }

    // members
    int indices[3];
    int pIndex;
    const PolygonSoup<2> *soup;
};

} // namespace fcpw

#include "silhouette_vertices.inl"