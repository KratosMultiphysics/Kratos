#pragma once

#include <fcpw/aggregates/bvh.h>
#ifdef FCPW_USE_EIGHT_WIDE_BRANCHING
    #define FCPW_MBVH_BRANCHING_FACTOR 8
    #define FCPW_MBVH_MAX_DEPTH 154
#else
    #define FCPW_MBVH_BRANCHING_FACTOR 4
    #define FCPW_MBVH_MAX_DEPTH 96
#endif

namespace fcpw {

template<size_t DIM>
struct MbvhNode {
    // constructor
    MbvhNode(): boxMin(FloatP<FCPW_MBVH_BRANCHING_FACTOR>(maxFloat)),
                boxMax(FloatP<FCPW_MBVH_BRANCHING_FACTOR>(minFloat)),
                child(maxInt) {}

    // members
    VectorP<FCPW_MBVH_BRANCHING_FACTOR, DIM> boxMin, boxMax;
    IntP<FCPW_MBVH_BRANCHING_FACTOR> child; // use sign to differentiate between inner and leaf nodes
};

template<size_t DIM>
struct MsnchNode {
    // constructor
    MsnchNode(): boxMin(FloatP<FCPW_MBVH_BRANCHING_FACTOR>(maxFloat)),
                 boxMax(FloatP<FCPW_MBVH_BRANCHING_FACTOR>(minFloat)),
                 coneAxis(FloatP<FCPW_MBVH_BRANCHING_FACTOR>(0.0f)),
                 coneHalfAngle(M_PI), coneRadius(0.0f),
                 child(maxInt), silhouetteChild(maxInt) {}

    // members
    VectorP<FCPW_MBVH_BRANCHING_FACTOR, DIM> boxMin, boxMax;
    VectorP<FCPW_MBVH_BRANCHING_FACTOR, DIM> coneAxis;
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> coneHalfAngle;
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> coneRadius;
    IntP<FCPW_MBVH_BRANCHING_FACTOR> child; // use sign to differentiate between inner and leaf nodes
    IntP<FCPW_MBVH_BRANCHING_FACTOR> silhouetteChild; // use sign to differentiate between inner and silhouette leaf nodes
};

template<size_t WIDTH, size_t DIM>
struct MbvhLeafNode {
    // members
    VectorP<WIDTH, DIM> positions[DIM];
    IntP<WIDTH> primitiveIndex;
};

template<size_t WIDTH, size_t DIM>
struct MbvhSilhouetteLeafNode {
    // members
    VectorP<WIDTH, DIM> positions[DIM + 1];
    IntP<WIDTH> primitiveIndex;
    MaskP<WIDTH> missingFace;
};

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType=Primitive<DIM>,
         typename SilhouetteType=SilhouettePrimitive<DIM>,
         typename NodeType=MbvhNode<DIM>,
         typename LeafNodeType=MbvhLeafNode<WIDTH, DIM>,
         typename SilhouetteLeafNodeType=MbvhSilhouetteLeafNode<WIDTH, DIM>>
class Mbvh: public Aggregate<DIM> {
public:
    // constructor
    Mbvh(std::vector<PrimitiveType *>& primitives_,
         std::vector<SilhouetteType *>& silhouettes_);

    // initialize Mbvh (clears any previously stored data)
    template<typename BvhNodeType>
    void initialize(const Bvh<DIM, BvhNodeType, PrimitiveType, SilhouetteType> *bvh);

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
    // collapses bvh into a mbvh
    template<typename BvhNodeType>
    int collapseBvh(const Bvh<DIM, BvhNodeType, PrimitiveType, SilhouetteType> *bvh,
                    int bvhNodeIndex, int parent, int depth);

    // determines whether mbvh node is a leaf node
    bool isLeafNode(const NodeType& node) const;

    // populates leaf nodes
    void populateLeafNodes();

    // populates leaf nodes
    void populateSilhouetteLeafNodes();

    // members
    int nNodes, nLeafs, nSilhouetteLeafs, maxDepth;
    std::vector<PrimitiveType *>& primitives;
    std::vector<SilhouetteType *>& silhouettes;
    std::vector<SilhouetteType *> silhouetteRefs;
    std::vector<NodeType> flatTree;
    std::vector<LeafNodeType> leafNodes;
    std::vector<SilhouetteLeafNodeType> silhouetteLeafNodes;
    bool primitiveTypeIsAggregate;
    bool primitiveTypeSupportsVectorizedQueries;
    bool silhouetteTypeSupportsVectorizedQueries;
    enoki::Array<int, DIM> range;
};

} // namespace fcpw

#include "mbvh.inl"