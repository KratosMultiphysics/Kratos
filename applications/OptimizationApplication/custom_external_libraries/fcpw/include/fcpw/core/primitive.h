#pragma once

#include <fcpw/core/ray.h>
#include <fcpw/core/bounding_volumes.h>
#include <fcpw/core/interaction.h>

namespace fcpw {

template<size_t DIM>
class Primitive {
public:
    // constructor
    Primitive() {}

    // destructor
    virtual ~Primitive() {}

    // returns bounding box
    virtual BoundingBox<DIM> boundingBox() const = 0;

    // returns centroid
    virtual Vector<DIM> centroid() const = 0;

    // returns surface area
    virtual float surfaceArea() const = 0;

    // returns signed volume
    virtual float signedVolume() const = 0;

    // splits the primitive along the provided coordinate and axis
    virtual void split(int dim, float splitCoord, BoundingBox<DIM>& boxLeft,
                       BoundingBox<DIM>& boxRight) const = 0;

    // get and set index
    virtual int getIndex() const = 0;
    virtual void setIndex(int index) = 0;
};

template<size_t DIM>
class GeometricPrimitive: public Primitive<DIM> {
public:
    // returns normal
    virtual Vector<DIM> normal(bool normalize=false) const = 0;

    // returns the normalized normal based on the local parameterization
    virtual Vector<DIM> normal(const Vector<DIM - 1>& uv) const = 0;

    // returns barycentric coordinates
    virtual Vector<DIM - 1> barycentricCoordinates(const Vector<DIM>& p) const = 0;

    // samples a random point on the geometric primitive and returns sampling pdf
    virtual float samplePoint(const Vector<DIM>& randNums, Vector<DIM - 1>& uv,
                              Vector<DIM>& p, Vector<DIM>& n) const = 0;

    // intersects with ray
    virtual bool intersect(const Ray<DIM>& r, Interaction<DIM>& i, bool checkForOcclusion=false) const = 0;

    // intersects with ray
    virtual int intersect(const Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                          bool checkForOcclusion=false, bool recordAllHits=false) const = 0;

    // intersects with sphere
    virtual bool intersect(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                           bool recordSurfaceArea=false) const = 0;

    // finds closest point to sphere center
    virtual bool findClosestPoint(const BoundingSphere<DIM>& s, Interaction<DIM>& i) const = 0;
};

template<size_t DIM>
class SilhouettePrimitive: public Primitive<DIM> {
public:
    // checks whether silhouette has adjacent face
    virtual bool hasFace(int fIndex) const = 0;

    // returns normal of adjacent face
    virtual Vector<DIM> normal(int fIndex, bool normalize=true) const = 0;

    // returns normalized silhouette normal
    virtual Vector<DIM> normal() const = 0;

    // returns signed volume
    float signedVolume() const {
        return 0.0f;
    }

    // splits the primitive along the provided coordinate and axis
    void split(int dim, float splitCoord, BoundingBox<DIM>& boxLeft,
               BoundingBox<DIM>& boxRight) const {
        // do nothing
    }

    // finds closest silhouette point to sphere center
    virtual bool findClosestSilhouettePoint(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                            bool flipNormalOrientation=false, float squaredMinRadius=0.0f,
                                            float precision=1e-3f) const = 0;
};

template<size_t DIM>
class Aggregate: public Primitive<DIM> {
public:
    // splits the primitive along the provided coordinate and axis
    void split(int dim, float splitCoord, BoundingBox<DIM>& boxLeft,
               BoundingBox<DIM>& boxRight) const {
        BoundingBox<DIM> box = this->boundingBox();

        if (box.pMin[dim] <= splitCoord) {
            boxLeft = box;
            boxLeft.pMax[dim] = splitCoord;
        }

        if (box.pMax[dim] >= splitCoord) {
            boxRight = box;
            boxRight.pMin[dim] = splitCoord;
        }
    }

    // intersects with ray
    bool intersect(Ray<DIM>& r, Interaction<DIM>& i, bool checkForOcclusion=false) const {
        int nodesVisited = 0;
        return this->intersectFromNode(r, i, 0, this->pIndex, nodesVisited, checkForOcclusion);
    }

    // intersects with ray
    int intersect(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                  bool checkForOcclusion=false, bool recordAllHits=false) const {
        int nodesVisited = 0;
        return this->intersectFromNode(r, is, 0, this->pIndex, nodesVisited,
                                       checkForOcclusion, recordAllHits);
    }

    // intersects with sphere
    int intersect(const BoundingSphere<DIM>& s, std::vector<Interaction<DIM>>& is,
                  bool recordOneHit=false) const {
        int nodesVisited = 0;
        return this->intersectFromNode(s, is, 0, this->pIndex, nodesVisited, recordOneHit);
    }

    // intersects with sphere
    int intersect(const BoundingSphere<DIM>& s, Interaction<DIM>& i, const Vector<DIM>& randNums,
                  const std::function<float(float)>& branchTraversalWeight={}) const {
        int nodesVisited = 0;
        return this->intersectFromNode(s, i, randNums, 0, this->pIndex,
                                       nodesVisited, branchTraversalWeight);
    }

    // finds closest point to sphere center
    bool findClosestPoint(BoundingSphere<DIM>& s, Interaction<DIM>& i, bool recordNormal=false) const {
        int nodesVisited = 0;
        return this->findClosestPointFromNode(s, i, 0, this->pIndex, nodesVisited, recordNormal);
    }

    // finds closest silhouette point to sphere center
    bool findClosestSilhouettePoint(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                    bool flipNormalOrientation=false, float squaredMinRadius=0.0f,
                                    float precision=1e-3f, bool recordNormal=false) const {
        int nodesVisited = 0;
        return this->findClosestSilhouettePointFromNode(s, i, 0, this->pIndex, nodesVisited,
                                                        flipNormalOrientation, squaredMinRadius,
                                                        precision, recordNormal);
    }

    // performs inside outside test for x
    // NOTE: assumes aggregate bounds watertight shape
    bool contains(const Vector<DIM>& x, bool useRayIntersection=true) const {
        if (useRayIntersection) {
            // do two intersection tests for robustness
            Vector<DIM> direction1 = Vector<DIM>::Zero();
            Vector<DIM> direction2 = Vector<DIM>::Zero();
            direction1[0] = 1;
            direction2[1] = 1;

            std::vector<Interaction<DIM>> is1;
            Ray<DIM> r1(x, direction1);
            int hits1 = this->intersect(r1, is1, false, true);

            std::vector<Interaction<DIM>> is2;
            Ray<DIM> r2(x, direction2);
            int hits2 = this->intersect(r2, is2, false, true);

            return hits1%2 == 1 && hits2%2 == 1;
        }

        Interaction<DIM> i;
        BoundingSphere<DIM> s(x, maxFloat);
        bool found = this->findClosestPoint(s, i, true);

        return i.signedDistance(x) < 0;
    }

    // checks whether there is a line of sight between xi and xj
    bool hasLineOfSight(const Vector<DIM>& xi, const Vector<DIM>& xj) const {
        Vector<DIM> direction = xj - xi;
        float dNorm = direction.norm();
        direction /= dNorm;

        Interaction<DIM> i;
        Ray<DIM> r(xi, direction, dNorm);
        int hits = this->intersect(r, i, true);

        return hits == 0;
    }

    // clamps x to the closest primitive this aggregate bounds
    void clampToBoundary(Vector<DIM>& x, float distanceUpperBound) const {
        Interaction<DIM> i;
        BoundingSphere<DIM> s(x, distanceUpperBound*distanceUpperBound);
        bool found = this->findClosestPoint(s, i);

        if (!found) {
            std::cerr << "Aggregate::clampToBoundary(): Cannot clamp to boundary since no "
                      << "closest point was found inside distance bound: " << distanceUpperBound
                      << std::endl;
        }

        if (i.distanceInfo == DistanceInfo::Bounded) {
            std::cerr << "Aggregate::clampToBoundary(): Cannot clamp to boundary since exact "
                      << "distance isn't available" << std::endl;
        }

        x = i.p;
    }

    // refits the aggregate
    virtual void refit() = 0;

    // intersects with ray, starting the traversal at the specified node in an aggregate
    // NOTE: interaction is invalid when checkForOcclusion is enabled
    virtual bool intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i,
                                   int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                                   bool checkForOcclusion=false) const = 0;

    // intersects with ray, starting the traversal at the specified node in an aggregate
    // NOTE: interactions are invalid when checkForOcclusion is enabled
    virtual int intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                  int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                                  bool checkForOcclusion=false, bool recordAllHits=false) const = 0;

    // intersects with sphere, starting the traversal at the specified node in an aggregate
    // NOTE: interactions contain primitive index
    virtual int intersectFromNode(const BoundingSphere<DIM>& s,
                                  std::vector<Interaction<DIM>>& is,
                                  int nodeStartIndex, int aggregateIndex,
                                  int& nodesVisited, bool recordOneHit=false) const = 0;

    // intersects with sphere, starting the traversal at the specified node in an aggregate
    // NOTE: interactions contain primitive index
    virtual int intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                  const Vector<DIM>& randNums, int nodeStartIndex,
                                  int aggregateIndex, int& nodesVisited,
                                  const std::function<float(float)>& branchTraversalWeight={}) const = 0;

    // finds closest point to sphere center, starting the traversal at the specified node in an aggregate
    virtual bool findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                          int nodeStartIndex, int aggregateIndex,
                                          int& nodesVisited, bool recordNormal=false) const = 0;

    // finds closest silhouette point to sphere center, starting the traversal at the specified node in an aggregate
    virtual bool findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                    int nodeStartIndex, int aggregateIndex,
                                                    int& nodesVisited, bool flipNormalOrientation=false,
                                                    float squaredMinRadius=0.0f, float precision=1e-3f,
                                                    bool recordNormal=false) const = 0;

    // get and set index
    int getIndex() const { return pIndex; }
    void setIndex(int index) { pIndex = index; }

    // member
    int pIndex;
};

template<size_t DIM>
class TransformedAggregate: public Aggregate<DIM> {
public:
    // constructor
    TransformedAggregate(const std::shared_ptr<Aggregate<DIM>>& aggregate_,
                         const Transform<DIM>& transform_):
                         aggregate(aggregate_), t(transform_), tInv(t.inverse()),
                         det(t.matrix().determinant()), sqrtDet(std::sqrt(det)) {}

    // returns bounding box
    BoundingBox<DIM> boundingBox() const {
        return aggregate->boundingBox().transform(t);
    }

    // returns centroid
    Vector<DIM> centroid() const {
        return t*aggregate->centroid();
    }

    // returns surface area
    float surfaceArea() const {
        // NOTE: this is an approximate estimate
        return sqrtDet*aggregate->surfaceArea();
    }

    // returns signed volume
    float signedVolume() const {
        // NOTE: this is an approximate estimate
        return det*aggregate->signedVolume();
    }

    // refits the transformed aggregate
    // NOTE: refitting of transformed aggregates is currently quite inefficient, since the
    // shared aggregate is refit every time this function is called
    void refit() {
        aggregate->refit();
    }

    // intersects with ray, starting the traversal at the specified node in an aggregate
    bool intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i,
                           int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                           bool checkForOcclusion=false) const {
        // apply inverse transform to ray
        Ray<DIM> rInv = r.transform(tInv);

        // intersect
        bool hit = aggregate->intersectFromNode(rInv, i, nodeStartIndex, aggregateIndex,
                                                nodesVisited, checkForOcclusion);

        // apply transform to ray and interactions
        r.tMax = rInv.transform(t).tMax;
        if (hit) i.applyTransform(t, tInv, r.o);

        nodesVisited++;
        return hit;
    }

    // intersects with ray, starting the traversal at the specified node in an aggregate
    int intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                          int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                          bool checkForOcclusion=false, bool recordAllHits=false) const {
        // apply inverse transform to ray
        Ray<DIM> rInv = r.transform(tInv);

        // intersect
        int hits = aggregate->intersectFromNode(rInv, is, nodeStartIndex, aggregateIndex,
                                                nodesVisited, checkForOcclusion, recordAllHits);

        // apply transform to ray and interactions
        r.tMax = rInv.transform(t).tMax;
        if (hits > 0) {
            for (int i = 0; i < (int)is.size(); i++) {
                is[i].applyTransform(t, tInv, r.o);
            }
        }

        nodesVisited++;
        return hits;
    }

    // intersects with sphere, starting the traversal at the specified node in an aggregate
    int intersectFromNode(const BoundingSphere<DIM>& s,
                          std::vector<Interaction<DIM>>& is,
                          int nodeStartIndex, int aggregateIndex,
                          int& nodesVisited, bool recordOneHit=false) const {
        // apply inverse transform to sphere
        BoundingSphere<DIM> sInv = s.transform(tInv);

        // intersect
        int hits = aggregate->intersectFromNode(sInv, is, nodeStartIndex, aggregateIndex,
                                                nodesVisited, recordOneHit);

        nodesVisited++;
        return hits;
    }

    // intersects with sphere, starting the traversal at the specified node in an aggregate
    int intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                          const Vector<DIM>& randNums, int nodeStartIndex,
                          int aggregateIndex, int& nodesVisited,
                          const std::function<float(float)>& branchTraversalWeight={}) const {
        // apply inverse transform to sphere
        BoundingSphere<DIM> sInv = s.transform(tInv);

        // intersect
        int hits = aggregate->intersectFromNode(sInv, i, randNums, nodeStartIndex, aggregateIndex,
                                                nodesVisited, branchTraversalWeight);

        // apply transform to interaction
        if (hits > 0) {
            i.applyTransform(t, tInv, s.c, false);
        }

        nodesVisited++;
        return hits;
    }

    // finds closest point to sphere center, starting the traversal at the specified node in an aggregate
    bool findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                  int nodeStartIndex, int aggregateIndex,
                                  int& nodesVisited, bool recordNormal=false) const {
        // apply inverse transform to sphere
        BoundingSphere<DIM> sInv = s.transform(tInv);

        // find closest point
        bool found = aggregate->findClosestPointFromNode(sInv, i, nodeStartIndex, aggregateIndex,
                                                         nodesVisited, recordNormal);

        // apply transform to sphere and interaction
        s.r2 = sInv.transform(t).r2;
        if (found) i.applyTransform(t, tInv, s.c);

        nodesVisited++;
        return found;
    }

    // finds closest silhouette point to sphere center, starting the traversal at the specified node in an aggregate
    bool findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                            int nodeStartIndex, int aggregateIndex,
                                            int& nodesVisited, bool flipNormalOrientation=false,
                                            float squaredMinRadius=0.0f, float precision=1e-3f,
                                            bool recordNormal=false) const {
        // apply inverse transform to sphere
        BoundingSphere<DIM> sInv = s.transform(tInv);
        BoundingSphere<DIM> sMin(s.c, squaredMinRadius);
        BoundingSphere<DIM> sMinInv = sMin.transform(tInv);

        // find closest silhouette point
        bool found = aggregate->findClosestSilhouettePointFromNode(sInv, i, nodeStartIndex, aggregateIndex,
                                                                   nodesVisited, flipNormalOrientation,
                                                                   sMinInv.r2, precision, recordNormal);

        // apply transform to sphere and interaction
        s.r2 = sInv.transform(t).r2;
        if (found) i.applyTransform(t, tInv, s.c);

        nodesVisited++;
        return found;
    }

    // performs inside outside test for x
    bool contains(const Vector<DIM>& x, bool useRayIntersection=true) const {
        return aggregate->contains(tInv*x, useRayIntersection);
    }

    // checks whether there is a line of sight between xi and xj
    bool hasLineOfSight(const Vector<DIM>& xi, const Vector<DIM>& xj) const {
        return aggregate->hasLineOfSight(tInv*xi, tInv*xj);
    }

    // clamps x to the closest primitive this aggregate bounds
    void clampToBoundary(Vector<DIM>& x, float distanceUpperBound) const {
        // apply inverse transform to x and distance bound
        Vector<DIM> xInv = tInv*x;
        if (distanceUpperBound < maxFloat) {
            Vector<DIM> direction = Vector<DIM>::Zero();
            direction[0] = 1;
            distanceUpperBound = (tInv*(x + direction*distanceUpperBound) - xInv).norm();
        }

        // clamp in object space and apply transform to x
        aggregate->clampToBoundary(xInv, distanceUpperBound);
        x = t*xInv;
    }

private:
    // members
    std::shared_ptr<Aggregate<DIM>> aggregate;
    Transform<DIM> t, tInv;
    float det, sqrtDet;
};

} // namespace fcpw