#pragma once

#include <fcpw/core/primitive.h>

namespace fcpw {

enum class BooleanOperation {
    Union,
    Intersection,
    Difference,
    None
};

template<size_t DIM,
         typename PrimitiveTypeLeft=Primitive<DIM>,
         typename PrimitiveTypeRight=Primitive<DIM>>
class CsgNode: public Aggregate<DIM> {
public:
    // constructor
    CsgNode(std::unique_ptr<PrimitiveTypeLeft> left_,
            std::unique_ptr<PrimitiveTypeRight> right_,
            const BooleanOperation& operation_);

    // refits the aggregate
    void refit();

    // returns bounding box
    BoundingBox<DIM> boundingBox() const;

    // returns centroid
    Vector<DIM> centroid() const;

    // returns surface area
    float surfaceArea() const;

    // returns signed volume
    float signedVolume() const;

    // intersects with ray, starting the traversal at the specified node in an aggregate
    // NOTE: interaction is invalid when checkForOcclusion is enabled
    bool intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                           int aggregateIndex, int& nodesVisited,
                           bool checkForOcclusion=false) const;

    // intersects with ray, starting the traversal at the specified node in an aggregate
    // NOTE: interactions are invalid when checkForOcclusion is enabled
    int intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                          int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                          bool checkForOcclusion=false, bool recordAllHits=false) const;

    // intersects with sphere, starting the traversal at the specified node in an aggregate
    // NOTE: interactions contain primitive index
    int intersectFromNode(const BoundingSphere<DIM>& s,
                          std::vector<Interaction<DIM>>& is,
                          int nodeStartIndex, int aggregateIndex,
                          int& nodesVisited, bool recordOneHit=false) const;

    // intersects with sphere, starting the traversal at the specified node in an aggregate
    // NOTE: interactions contain primitive index
    int intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                          const Vector<DIM>& randNums, int nodeStartIndex,
                          int aggregateIndex, int& nodesVisited,
                          const std::function<float(float)>& branchTraversalWeight={}) const;

    // finds closest point to sphere center, starting the traversal at the specified node in an aggregate
    bool findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                  int nodeStartIndex, int aggregateIndex,
                                  int& nodesVisited, bool recordNormal=false) const;

    // finds closest silhouette point to sphere center, starting the traversal at the specified node in an aggregate
    bool findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                            int nodeStartIndex, int aggregateIndex,
                                            int& nodesVisited, bool flipNormalOrientation=false,
                                            float squaredMinRadius=0.0f, float precision=1e-3f,
                                            bool recordNormal=false) const;

private:
    // computes bounding box in world sparce
    void computeBoundingBox();

    // computes interactions for ray intersection
    void computeInteractions(const std::vector<Interaction<DIM>>& isLeft,
                             const std::vector<Interaction<DIM>>& isRight,
                             std::vector<Interaction<DIM>>& is) const;

    // members
    std::unique_ptr<PrimitiveTypeLeft> left;
    std::unique_ptr<PrimitiveTypeRight> right;
    BooleanOperation operation;
    BoundingBox<DIM> box;
    bool leftPrimitiveTypeIsAggregate, rightPrimitiveTypeIsAggregate;
};

} // namespace fcpw

#include "csg_node.inl"