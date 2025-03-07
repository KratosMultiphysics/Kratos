#pragma once

#include <fcpw/core/primitive.h>
#include <tuple>
#define FCPW_BVH_MAX_DEPTH 64

namespace fcpw {
// modified version of https://github.com/brandonpelfrey/Fast-BVH and
// https://github.com/straaljager/GPU-path-tracing-with-CUDA-tutorial-4

enum class CostHeuristic {
    LongestAxisCenter,
    SurfaceArea,
    OverlapSurfaceArea,
    Volume,
    OverlapVolume
};

template<size_t DIM>
struct BvhNode {
    // constructor
    BvhNode(): nReferences(0) {}

    // members
    BoundingBox<DIM> box;
    union {
        int referenceOffset;
        int secondChildOffset;
    };
    int nReferences;
};

template<size_t DIM>
struct SnchNode {
    // constructor
    SnchNode(): nReferences(0), nSilhouetteReferences(0) {}

    // members
    BoundingBox<DIM> box;
    BoundingCone<DIM> cone;
    union {
        int referenceOffset;
        int secondChildOffset;
    };
    int silhouetteReferenceOffset;
    int nReferences;
    int nSilhouetteReferences;
};

struct TraversalStack {
    // constructors
    TraversalStack(): node(-1), distance(0.0f) {}
    TraversalStack(int node_, float distance_): node(node_), distance(distance_) {}

    // members
    int node; // node index
    float distance; // minimum distance (parametric, squared, ...) to this node
};

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
using SortPositionsFunc = std::function<void(const std::vector<NodeType>&,
                                             std::vector<PrimitiveType *>&,
                                             std::vector<SilhouetteType *>&)>;

template<size_t DIM,
         typename NodeType=BvhNode<DIM>,
         typename PrimitiveType=Primitive<DIM>,
         typename SilhouetteType=SilhouettePrimitive<DIM>>
class Bvh: public Aggregate<DIM> {
public:
    // constructor
    Bvh(const CostHeuristic& costHeuristic_,
        std::vector<PrimitiveType *>& primitives_,
        std::vector<SilhouetteType *>& silhouettes_,
        SortPositionsFunc<DIM, NodeType, PrimitiveType, SilhouetteType> sortPositions_={},
        const std::function<bool(float, int)>& ignoreSilhouette_={},
        bool packLeaves_=false, int leafSize_=4, int nBuckets_=8);

    // refits the aggregate
    void refit();

    // prints statistics
    void printStats() const;

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

protected:
    // computes split cost based on heuristic
    float computeSplitCost(const BoundingBox<DIM>& boxLeft,
                           const BoundingBox<DIM>& boxRight,
                           int nReferencesLeft, int nReferencesRight,
                           int depth) const;

    // computes object split
    float computeObjectSplit(const BoundingBox<DIM>& nodeBoundingBox,
                             const BoundingBox<DIM>& nodeCentroidBox,
                             const std::vector<BoundingBox<DIM>>& referenceBoxes,
                             const std::vector<Vector<DIM>>& referenceCentroids,
                             int depth, int nodeStart, int nodeEnd,
                             int& splitDim, float& splitCoord);

    // performs object split
    int performObjectSplit(int nodeStart, int nodeEnd, int splitDim, float splitCoord,
                           std::vector<BoundingBox<DIM>>& referenceBoxes,
                           std::vector<Vector<DIM>>& referenceCentroids);

    // helper function to build binary tree
    void buildRecursive(std::vector<BoundingBox<DIM>>& referenceBoxes,
                        std::vector<Vector<DIM>>& referenceCentroids,
                        std::vector<NodeType>& buildNodes,
                        int parent, int start, int end, int depth);

    // builds binary tree
    void build();

    // assigns geometric data (e.g. cones) to nodes
    void assignGeometricDataToNodes(const std::function<bool(float, int)>& ignoreSilhouette);

    // processes subtree for intersection
    bool processSubtreeForIntersection(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                                       int aggregateIndex, bool checkForOcclusion,
                                       TraversalStack *subtree, float *boxHits,
                                       bool& didHit, int& nodesVisited) const;

    // processes subtree for intersection
    bool processSubtreeForIntersection(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                       int nodeStartIndex, int aggregateIndex, bool checkForOcclusion,
                                       bool recordAllHits, TraversalStack *subtree,
                                       float *boxHits, int& hits, int& nodesVisited) const;

    // processes subtree for intersection
    float processSubtreeForIntersection(const BoundingSphere<DIM>& s, std::vector<Interaction<DIM>>& is,
                                        int nodeStartIndex, int aggregateIndex, bool recordOneHit,
                                        TraversalStack *subtree, float *boxHits, int& hits, int& nodesVisited) const;

    // processes subtree for intersection
    void processSubtreeForIntersection(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                       const Vector<DIM>& randNums, int nodeStartIndex, int aggregateIndex,
                                       const std::function<float(float)>& branchTraversalWeight,
                                       int nodeIndex, float traversalPdf, float *boxHits,
                                       int& hits, int& nodesVisited) const;

    // processes subtree for closest point
    void processSubtreeForClosestPoint(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                       int nodeStartIndex, int aggregateIndex,
                                       bool recordNormal, TraversalStack *subtree,
                                       float *boxHits, bool& notFound, int& nodesVisited) const;

    // members
    CostHeuristic costHeuristic;
    int nNodes, nLeafs, leafSize, nBuckets, maxDepth, depthGuess;
    std::vector<std::pair<BoundingBox<DIM>, int>> buckets, rightBuckets;
    std::vector<PrimitiveType *>& primitives;
    std::vector<SilhouetteType *>& silhouettes;
    std::vector<SilhouetteType *> silhouetteRefs;
    std::vector<NodeType> flatTree;
    bool packLeaves, primitiveTypeIsAggregate;

    template<size_t A, size_t B, typename C, typename D, typename E, typename F, typename G>
    friend class Mbvh;

    template<size_t A, typename B, typename C, typename D, typename E, typename F, typename G>
    friend class CPUBvhDataExtractor;
    template<size_t A, typename B, typename C, typename D>
    friend class CPUBvhRefitDataExtractor;
};

} // namespace fcpw

#include "bvh.inl"