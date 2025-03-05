#include <fcpw/core/wide_query_operations.h>

namespace fcpw {

template<typename BvhNodeType, typename MbvhNodeType>
inline void assignGeometricDataToNode(const BvhNodeType& bvhNode, MbvhNodeType& mbvhNode, int index)
{
    // do nothing
}

template<size_t DIM>
inline void assignGeometricDataToNode(const SnchNode<DIM>& bvhNode, MsnchNode<DIM>& mbvhNode, int index)
{
    // assign bvh node's bounding cone to mbvh node
    for (size_t j = 0; j < DIM; j++) {
        mbvhNode.coneAxis[j][index] = bvhNode.cone.axis[j];
    }

    mbvhNode.coneHalfAngle[index] = bvhNode.cone.halfAngle;
    mbvhNode.coneRadius[index] = bvhNode.cone.radius;
}

template<typename BvhNodeType, typename MbvhNodeType>
inline void assignSilhouetteLeafRangeToNode(const BvhNodeType& bvhNode, MbvhNodeType& mbvhNode,
                                            size_t WIDTH, int& nSilhouetteLeafs)
{
    // do nothing
}

template<size_t DIM>
inline void assignSilhouetteLeafRangeToNode(const SnchNode<DIM>& bvhNode, MsnchNode<DIM>& mbvhNode,
                                            size_t WIDTH, int& nSilhouetteLeafs)
{
    if (bvhNode.nSilhouetteReferences > 0) {
        mbvhNode.silhouetteChild[0] = -(nSilhouetteLeafs + 1); // negative value indicates that node is a leaf
        mbvhNode.silhouetteChild[1] = bvhNode.nSilhouetteReferences/WIDTH;
        if (bvhNode.nSilhouetteReferences%WIDTH != 0) mbvhNode.silhouetteChild[1] += 1;
        mbvhNode.silhouetteChild[2] = bvhNode.silhouetteReferenceOffset;
        mbvhNode.silhouetteChild[3] = bvhNode.nSilhouetteReferences;
        nSilhouetteLeafs += mbvhNode.silhouetteChild[1];

    } else {
        mbvhNode.silhouetteChild[0] = 0;
        mbvhNode.silhouetteChild[1] = 0;
        mbvhNode.silhouetteChild[2] = 0;
        mbvhNode.silhouetteChild[3] = 0;
    }
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
template<typename BvhNodeType>
inline int Mbvh<WIDTH, DIM,
                PrimitiveType,
                SilhouetteType,
                NodeType,
                LeafNodeType,
                SilhouetteLeafNodeType>::collapseBvh(const Bvh<DIM, BvhNodeType, PrimitiveType, SilhouetteType> *bvh,
                                                     int bvhNodeIndex, int parent, int depth)
{
    const BvhNodeType& bvhNode = bvh->flatTree[bvhNodeIndex];
    maxDepth = std::max(depth, maxDepth);

    // create mbvh node
    NodeType mbvhNode;
    int mbvhNodeIndex = nNodes++;
    flatTree.emplace_back(mbvhNode);

    if (bvhNode.nReferences > 0) {
        // bvh node is a leaf node; assign mbvh node its reference indices
        NodeType& mbvhNode = flatTree[mbvhNodeIndex];
        mbvhNode.child[0] = -(nLeafs + 1); // negative value indicates that node is a leaf
        mbvhNode.child[1] = bvhNode.nReferences/WIDTH;
        if (bvhNode.nReferences%WIDTH != 0) mbvhNode.child[1] += 1;
        mbvhNode.child[2] = bvhNode.referenceOffset;
        mbvhNode.child[3] = bvhNode.nReferences;
        nLeafs += mbvhNode.child[1];
        assignSilhouetteLeafRangeToNode(bvhNode, mbvhNode, WIDTH, nSilhouetteLeafs);

    } else {
        // bvh node is an inner node, flatten it
        int nNodesToCollapse = 2;
        int nodesToCollapse[FCPW_MBVH_BRANCHING_FACTOR];
        nodesToCollapse[0] = bvhNodeIndex + bvhNode.secondChildOffset;
        nodesToCollapse[1] = bvhNodeIndex + 1;
        bool noMoreNodesToCollapse = false;

        while (nNodesToCollapse < FCPW_MBVH_BRANCHING_FACTOR && !noMoreNodesToCollapse) {
            // find the (non-leaf) node entry with the largest surface area
            float maxSurfaceArea = minFloat;
            int maxIndex = -1;

            for (int i = 0; i < nNodesToCollapse; i++) {
                int bvhNodeIndex = nodesToCollapse[i];
                const BvhNodeType& bvhNode = bvh->flatTree[bvhNodeIndex];

                if (bvhNode.nReferences == 0) {
                    float surfaceArea = bvhNode.box.surfaceArea();

                    if (maxSurfaceArea < surfaceArea) {
                        maxSurfaceArea = surfaceArea;
                        maxIndex = i;
                    }
                }
            }

            if (maxIndex == -1) {
                // no more nodes to collapse
                noMoreNodesToCollapse = true;

            } else {
                // remove the selected node from the list, and add its two children
                int bvhNodeIndex = nodesToCollapse[maxIndex];
                const BvhNodeType& bvhNode = bvh->flatTree[bvhNodeIndex];

                nodesToCollapse[maxIndex] = bvhNodeIndex + bvhNode.secondChildOffset;
                nodesToCollapse[nNodesToCollapse] = bvhNodeIndex + 1;
                nNodesToCollapse++;
            }
        }

        // collapse the nodes
        std::sort(nodesToCollapse, nodesToCollapse + nNodesToCollapse);
        for (int i = 0; i < nNodesToCollapse; i++) {
            int bvhNodeIndex = nodesToCollapse[i];
            const BvhNodeType& bvhNode = bvh->flatTree[bvhNodeIndex];

            // assign bvh node's bounding box to mbvh node
            for (size_t j = 0; j < DIM; j++) {
                flatTree[mbvhNodeIndex].boxMin[j][i] = bvhNode.box.pMin[j];
                flatTree[mbvhNodeIndex].boxMax[j][i] = bvhNode.box.pMax[j];
            }

            // assign bvh node's geometric data (e.g. cones) to mbvh node
            assignGeometricDataToNode(bvhNode, flatTree[mbvhNodeIndex], i);

            // collapse bvh node
            flatTree[mbvhNodeIndex].child[i] = collapseBvh(bvh, bvhNodeIndex, mbvhNodeIndex, depth + 1);
        }
    }

    return mbvhNodeIndex;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline bool Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::isLeafNode(const NodeType& node) const
{
    return node.child[0] < 0;
}

template<typename PrimitiveType,
         typename NodeType,
         typename LeafNodeType>
inline void populateLeafNode(const NodeType& node,
                             const std::vector<PrimitiveType *>& primitives,
                             std::vector<LeafNodeType>& leafNodes, size_t WIDTH)
{
    std::cerr << "populateLeafNode(): WIDTH: " << WIDTH << " not supported" << std::endl;
    exit(EXIT_FAILURE);
}

template<typename NodeType, typename LeafNodeType>
inline void populateLeafNode(const NodeType& node,
                             const std::vector<LineSegment *>& primitives,
                             std::vector<LeafNodeType>& leafNodes, size_t WIDTH)
{
    int leafOffset = -node.child[0] - 1;
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];

    // populate leaf node with line segments
    for (int p = 0; p < nReferences; p++) {
        int referenceIndex = referenceOffset + p;
        int leafIndex = leafOffset + p/WIDTH;
        int w = p%WIDTH;

        const LineSegment *lineSegment = primitives[referenceIndex];
        const Vector2& pa = lineSegment->soup->positions[lineSegment->indices[0]];
        const Vector2& pb = lineSegment->soup->positions[lineSegment->indices[1]];

        leafNodes[leafIndex].primitiveIndex[w] = lineSegment->getIndex();
        for (int i = 0; i < 2; i++) {
            leafNodes[leafIndex].positions[0][i][w] = pa[i];
            leafNodes[leafIndex].positions[1][i][w] = pb[i];
        }
    }
}

template<typename NodeType, typename LeafNodeType>
inline void populateLeafNode(const NodeType& node,
                             const std::vector<Triangle *>& primitives,
                             std::vector<LeafNodeType>& leafNodes, size_t WIDTH)
{
    int leafOffset = -node.child[0] - 1;
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];

    // populate leaf node with triangles
    for (int p = 0; p < nReferences; p++) {
        int referenceIndex = referenceOffset + p;
        int leafIndex = leafOffset + p/WIDTH;
        int w = p%WIDTH;

        const Triangle *triangle = primitives[referenceIndex];
        const Vector3& pa = triangle->soup->positions[triangle->indices[0]];
        const Vector3& pb = triangle->soup->positions[triangle->indices[1]];
        const Vector3& pc = triangle->soup->positions[triangle->indices[2]];

        leafNodes[leafIndex].primitiveIndex[w] = triangle->getIndex();
        for (int i = 0; i < 3; i++) {
            leafNodes[leafIndex].positions[0][i][w] = pa[i];
            leafNodes[leafIndex].positions[1][i][w] = pb[i];
            leafNodes[leafIndex].positions[2][i][w] = pc[i];
        }
    }
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline void Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::populateLeafNodes()
{
    if (primitiveTypeSupportsVectorizedQueries) {
        leafNodes.resize(nLeafs);

        for (int i = 0; i < nNodes; i++) {
            const NodeType& node = flatTree[i];
            if (isLeafNode(node)) populateLeafNode(node, primitives, leafNodes, WIDTH);
        }
    }
}

template<typename SilhouetteType,
         typename NodeType,
         typename SilhouetteLeafNodeType>
inline void populateSilhouetteLeafNode(const NodeType& node,
                                       const std::vector<SilhouetteType *>& silhouettes,
                                       std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                       size_t WIDTH)
{
    std::cerr << "populateSilhouetteLeafNode(): WIDTH: " << WIDTH << " not supported" << std::endl;
    exit(EXIT_FAILURE);
}

template<typename SilhouetteLeafNodeType>
inline void populateSilhouetteLeafNode(const MsnchNode<2>& node,
                                       const std::vector<SilhouetteVertex *>& silhouettes,
                                       std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                       size_t WIDTH)
{
    int silhouetteLeafOffset = -node.silhouetteChild[0] - 1;
    int silhouetteReferenceOffset = node.silhouetteChild[2];
    int nSilhouetteReferences = node.silhouetteChild[3];

    // populate silhouette leaf node with silhouette vertices
    for (int p = 0; p < nSilhouetteReferences; p++) {
        int referenceIndex = silhouetteReferenceOffset + p;
        int leafIndex = silhouetteLeafOffset + p/WIDTH;
        int w = p%WIDTH;

        const SilhouetteVertex *silhouetteVertex = silhouettes[referenceIndex];
        SilhouetteLeafNodeType& silhouetteLeafNode = silhouetteLeafNodes[leafIndex];
        silhouetteLeafNode.primitiveIndex[w] = silhouetteVertex->getIndex();
        silhouetteLeafNode.missingFace[w] = !silhouetteVertex->hasFace(0) || !silhouetteVertex->hasFace(1);

        const Vector2& pb = silhouetteVertex->soup->positions[silhouetteVertex->indices[1]];
        for (int i = 0; i < 2; i++) {
            silhouetteLeafNode.positions[1][i][w] = pb[i];
        }

        if (silhouetteVertex->hasFace(0)) {
            Vector2 n0 = silhouetteVertex->normal(0);
            for (int i = 0; i < 2; i++) {
                silhouetteLeafNode.positions[0][i][w] = n0[i];
            }
        }

        if (silhouetteVertex->hasFace(1)) {
            Vector2 n1 = silhouetteVertex->normal(1);
            for (int i = 0; i < 2; i++) {
                silhouetteLeafNode.positions[2][i][w] = n1[i];
            }
        }
    }
}

template<typename SilhouetteLeafNodeType>
inline void populateSilhouetteLeafNode(const MsnchNode<3>& node,
                                       const std::vector<SilhouetteEdge *>& silhouettes,
                                       std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                       size_t WIDTH)
{
    int silhouetteLeafOffset = -node.silhouetteChild[0] - 1;
    int silhouetteReferenceOffset = node.silhouetteChild[2];
    int nSilhouetteReferences = node.silhouetteChild[3];

    // populate silhouette leaf node with silhouette edges
    for (int p = 0; p < nSilhouetteReferences; p++) {
        int referenceIndex = silhouetteReferenceOffset + p;
        int leafIndex = silhouetteLeafOffset + p/WIDTH;
        int w = p%WIDTH;

        const SilhouetteEdge *silhouetteEdge = silhouettes[referenceIndex];
        SilhouetteLeafNodeType& silhouetteLeafNode = silhouetteLeafNodes[leafIndex];
        silhouetteLeafNode.primitiveIndex[w] = silhouetteEdge->getIndex();
        silhouetteLeafNode.missingFace[w] = !silhouetteEdge->hasFace(0) || !silhouetteEdge->hasFace(1);

        const Vector3& pb = silhouetteEdge->soup->positions[silhouetteEdge->indices[1]];
        const Vector3& pc = silhouetteEdge->soup->positions[silhouetteEdge->indices[2]];
        for (int i = 0; i < 3; i++) {
            silhouetteLeafNode.positions[1][i][w] = pb[i];
            silhouetteLeafNode.positions[2][i][w] = pc[i];
        }

        if (silhouetteEdge->hasFace(0)) {
            Vector3 n0 = silhouetteEdge->normal(0);
            for (int i = 0; i < 3; i++) {
                silhouetteLeafNode.positions[0][i][w] = n0[i];
            }
        }

        if (silhouetteEdge->hasFace(1)) {
            Vector3 n1 = silhouetteEdge->normal(1);
            for (int i = 0; i < 3; i++) {
                silhouetteLeafNode.positions[3][i][w] = n1[i];
            }
        }
    }
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline void Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::populateSilhouetteLeafNodes()
{
    if (silhouetteTypeSupportsVectorizedQueries) {
        silhouetteLeafNodes.resize(nSilhouetteLeafs);

        for (int i = 0; i < nNodes; i++) {
            const NodeType& node = flatTree[i];
            if (isLeafNode(node)) {
                populateSilhouetteLeafNode(node, silhouetteRefs, silhouetteLeafNodes, WIDTH);
            }
        }
    }
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline Mbvh<WIDTH, DIM,
            PrimitiveType,
            SilhouetteType,
            NodeType,
            LeafNodeType,
            SilhouetteLeafNodeType>::Mbvh(std::vector<PrimitiveType *>& primitives_,
                                          std::vector<SilhouetteType *>& silhouettes_):
nNodes(0),
nLeafs(0),
nSilhouetteLeafs(0),
maxDepth(0),
primitives(primitives_),
silhouettes(silhouettes_),
range(enoki::arange<enoki::Array<int, DIM>>())
{
    primitiveTypeIsAggregate = std::is_base_of<Aggregate<DIM>, PrimitiveType>::value;
    primitiveTypeSupportsVectorizedQueries = std::is_same<PrimitiveType, LineSegment>::value ||
                                             std::is_same<PrimitiveType, Triangle>::value;
    silhouetteTypeSupportsVectorizedQueries = std::is_same<SilhouetteType, SilhouetteVertex>::value ||
                                              std::is_same<SilhouetteType, SilhouetteEdge>::value;
    static_assert(FCPW_MBVH_BRANCHING_FACTOR == 4 || FCPW_MBVH_BRANCHING_FACTOR == 8,
                  "Branching factor must be atleast 4");
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
template<typename BvhNodeType>
inline void Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::initialize(const Bvh<DIM, BvhNodeType, PrimitiveType, SilhouetteType> *bvh)
{
    // clear previous data
    nNodes = 0;
    nLeafs = 0;
    nSilhouetteLeafs = 0;
    maxDepth = 0;
    silhouetteRefs.clear();
    flatTree.clear();
    leafNodes.clear();
    silhouetteLeafNodes.clear();

    // collapse bvh
    silhouetteRefs = bvh->silhouetteRefs;
    collapseBvh(bvh, 0, 0xfffffffc, 0);

    // populate leaf nodes if primitive type is supported
    populateLeafNodes();

    // populate silhouette leaf nodes if primitive type is supported
    populateSilhouetteLeafNodes();
}

template<size_t DIM, typename NodeType, typename SilhouetteType>
inline void computeBoundingCone(const std::vector<SilhouetteType *>& silhouetteRefs,
                                const BoundingBox<DIM>& box, NodeType& node,
                                BoundingCone<DIM>& cone)
{
    // do nothing
}

template<size_t DIM, typename SilhouetteType>
inline void computeBoundingCone(const std::vector<SilhouetteType *>& silhouetteRefs,
                                const BoundingBox<DIM>& box, MsnchNode<DIM>& node,
                                BoundingCone<DIM>& cone)
{
    int silhouetteReferenceOffset = node.silhouetteChild[2];
    int nSilhouetteReferences = node.silhouetteChild[3];

    cone = computeBoundingCone<DIM, SilhouetteType>(silhouetteRefs, box.centroid(),
                                                    nSilhouetteReferences,
                                                    silhouetteReferenceOffset);
}

template<size_t DIM, typename NodeType>
inline void assignBoundingCone(const BoundingCone<DIM>& cone, NodeType& node, int index)
{
    // do nothing
}

template<size_t DIM>
inline void assignBoundingCone(const BoundingCone<DIM>& cone, MsnchNode<DIM>& node, int index)
{
    for (size_t i = 0; i < DIM; i++) {
        node.coneAxis[i][index] = cone.axis[i];
    }

    node.coneHalfAngle[index] = cone.halfAngle;
    node.coneRadius[index] = cone.radius;
}

template<size_t DIM, typename NodeType>
inline void mergeBoundingCones(const BoundingCone<DIM>& coneA, const BoundingCone<DIM>& coneB,
                               const BoundingBox<DIM>& boxA, const BoundingBox<DIM>& boxB,
                               const BoundingBox<DIM>& mergedBox, NodeType& node,
                               BoundingCone<DIM>& cone)
{
    // do nothing
}

template<size_t DIM>
inline void mergeBoundingCones(const BoundingCone<DIM>& coneA, const BoundingCone<DIM>& coneB,
                               const BoundingBox<DIM>& boxA, const BoundingBox<DIM>& boxB,
                               const BoundingBox<DIM>& mergedBox, MsnchNode<DIM>& node,
                               BoundingCone<DIM>& cone)
{
    cone = mergeBoundingCones<DIM>(coneA, coneB,
                                   boxA.centroid(),
                                   boxB.centroid(),
                                   mergedBox.centroid());
}

template<size_t WIDTH,
         size_t DIM,
         typename NodeType,
         typename PrimitiveType,
         typename SilhouetteType>
inline std::pair<BoundingBox<DIM>, BoundingCone<DIM>> refitRecursive(const std::vector<PrimitiveType *>& primitives,
                                                                     const std::vector<SilhouetteType *>& silhouetteRefs,
                                                                     std::vector<NodeType>& flatTree, int nodeIndex)
{
    BoundingBox<DIM> box;
    BoundingCone<DIM> cone;
    cone.halfAngle = -M_PI;
    NodeType& node(flatTree[nodeIndex]);

    if (node.child[0] < 0) { // leaf
        // compute bounding box
        int referenceOffset = node.child[2];
        int nReferences = node.child[3];

        for (int p = 0; p < nReferences; p++) {
            int referenceIndex = referenceOffset + p;
            const PrimitiveType *prim = primitives[referenceIndex];

            box.expandToInclude(prim->boundingBox());
        }

        // compute bounding cone
        computeBoundingCone(silhouetteRefs, box, node, cone);

    } else { // not a leaf
        for (int w = 0; w < FCPW_MBVH_BRANCHING_FACTOR; w++) {
            if (node.child[w] != maxInt) {
                // refit child
                std::pair<BoundingBox<DIM>, BoundingCone<DIM>> childBoxCone =
                    refitRecursive<WIDTH, DIM, NodeType, PrimitiveType, SilhouetteType>(
                        primitives, silhouetteRefs, flatTree, node.child[w]);

                // expand bounding box
                BoundingBox<DIM> currentBox = box;
                BoundingBox<DIM> childBox = childBoxCone.first;
                for (size_t i = 0; i < DIM; i++) {
                    node.boxMin[i][w] = childBox.pMin[i];
                    node.boxMax[i][w] = childBox.pMax[i];
                }
                box.expandToInclude(childBox);

                // expand bounding cone
                BoundingCone<DIM> childCone = childBoxCone.second;
                assignBoundingCone(childCone, node, w);
                mergeBoundingCones(cone, childCone, currentBox, childBox, box, node, cone);
            }
        }
    }

    return std::make_pair(box, cone);
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline void Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::refit()
{
    // refit primitives if they are aggregates
    if (primitiveTypeIsAggregate) {
        for (int p = 0; p < (int)primitives.size(); p++) {
            Aggregate<DIM> *aggregate = reinterpret_cast<Aggregate<DIM> *>(primitives[p]);
            aggregate->refit();
        }
    }

    // update leaf and silhouette leaf nodes if primitive type is supported
    populateLeafNodes();
    populateSilhouetteLeafNodes();

    // update flatTree
    if (nNodes > 0) {
        refitRecursive<WIDTH, DIM, NodeType, PrimitiveType, SilhouetteType>(
            primitives, silhouetteRefs, flatTree, 0);
    }
}

template<typename NodeType>
inline void updateSilhouetteLeafInfo(const NodeType& mbvhNode, size_t WIDTH,
                                     float& nSilhouetteLeafsNotFull)
{
    // do nothing
}

template<size_t DIM>
inline void updateSilhouetteLeafInfo(const MsnchNode<DIM>& mbvhNode, size_t WIDTH,
                                     float& nSilhouetteLeafsNotFull)
{
    if (mbvhNode.silhouetteChild[3]%WIDTH != 0) {
        nSilhouetteLeafsNotFull += 1.0f;
    }
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline void Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::printStats() const
{
    // count not-full nodes
    float nLeafsNotFull = 0.0f;
    float nSilhouetteLeafsNotFull = 0.0f;
    float nNodesNotFull = 0.0f;
    int nInnerNodes = 0;

    for (int i = 0; i < nNodes; i++) {
        const NodeType& node = flatTree[i];

        if (isLeafNode(node)) {
            if (node.child[3]%WIDTH != 0) {
                nLeafsNotFull += 1.0f;
            }

            updateSilhouetteLeafInfo(node, WIDTH, nSilhouetteLeafsNotFull);

        } else {
            nInnerNodes++;
            for (int w = 0; w < FCPW_MBVH_BRANCHING_FACTOR; w++) {
                if (node.child[w] == maxInt) {
                    nNodesNotFull += 1.0f;
                    break;
                }
            }
        }
    }

    std::cout << FCPW_MBVH_BRANCHING_FACTOR << "-BVH stats: "
              << nNodes << " nodes, "
              << nLeafs << " leaves, "
              << nSilhouetteLeafs << " silhouette leaves, "
              << (nNodesNotFull*100/nInnerNodes) << "% nodes, "
              << (nLeafsNotFull*100/nLeafs) << "% leaves not full & "
              << (nSilhouetteLeafs > 0 ? (nSilhouetteLeafsNotFull*100/nSilhouetteLeafs) : 0) << "% silhouette leaves not full, "
              << maxDepth << " max depth"
              << std::endl;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline BoundingBox<DIM> Mbvh<WIDTH, DIM,
                             PrimitiveType,
                             SilhouetteType,
                             NodeType,
                             LeafNodeType,
                             SilhouetteLeafNodeType>::boundingBox() const
{
    BoundingBox<DIM> box;
    if (flatTree.size() == 0) return box;

    enoki::scatter(box.pMin.data(), enoki::hmin_inner(flatTree[0].boxMin), range);
    enoki::scatter(box.pMax.data(), enoki::hmax_inner(flatTree[0].boxMax), range);

    return box;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline Vector<DIM> Mbvh<WIDTH, DIM,
                        PrimitiveType,
                        SilhouetteType,
                        NodeType,
                        LeafNodeType,
                        SilhouetteLeafNodeType>::centroid() const
{
    Vector<DIM> c = Vector<DIM>::Zero();
    int nPrimitives = (int)primitives.size();

    for (int p = 0; p < nPrimitives; p++) {
        c += primitives[p]->centroid();
    }

    return c/nPrimitives;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline float Mbvh<WIDTH, DIM,
                  PrimitiveType,
                  SilhouetteType,
                  NodeType,
                  LeafNodeType,
                  SilhouetteLeafNodeType>::surfaceArea() const
{
    float area = 0.0f;
    for (int p = 0; p < (int)primitives.size(); p++) {
        area += primitives[p]->surfaceArea();
    }

    return area;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline float Mbvh<WIDTH, DIM,
                  PrimitiveType,
                  SilhouetteType,
                  NodeType,
                  LeafNodeType,
                  SilhouetteLeafNodeType>::signedVolume() const
{
    float volume = 0.0f;
    for (int p = 0; p < (int)primitives.size(); p++) {
        volume += primitives[p]->signedVolume();
    }

    return volume;
}

template<size_t WIDTH>
inline void enqueueNodes(const IntP<WIDTH>& child, const FloatP<WIDTH>& tMin,
                         const FloatP<WIDTH>& tMax, const MaskP<WIDTH>& mask, float minDist,
                         float& tMaxMin, int& stackPtr, TraversalStack *subtree)
{
    // enqueue nodes
    int closestIndex = -1;
    for (int w = 0; w < WIDTH; w++) {
        if (mask[w]) {
            stackPtr++;
            subtree[stackPtr].node = child[w];
            subtree[stackPtr].distance = tMin[w];
            tMaxMin = std::min(tMaxMin, tMax[w]);

            if (tMin[w] < minDist) {
                closestIndex = stackPtr;
                minDist = tMin[w];
            }
        }
    }

    // put closest node first
    if (closestIndex != -1) {
        std::swap(subtree[stackPtr], subtree[closestIndex]);
    }
}

inline void sortOrder4(const FloatP<4>& t, int& a, int& b, int& c, int& d)
{
    // source: https://stackoverflow.com/questions/25070577/sort-4-numbers-without-array
    int tmp;
    if (t[a] < t[b]) { tmp = a; a = b; b = tmp; }
    if (t[c] < t[d]) { tmp = c; c = d; d = tmp; }
    if (t[a] < t[c]) { tmp = a; a = c; c = tmp; }
    if (t[b] < t[d]) { tmp = b; b = d; d = tmp; }
    if (t[b] < t[c]) { tmp = b; b = c; c = tmp; }
}

template<>
inline void enqueueNodes<4>(const IntP<4>& child, const FloatP<4>& tMin,
                            const FloatP<4>& tMax, const MaskP<4>& mask, float minDist,
                            float& tMaxMin, int& stackPtr, TraversalStack *subtree)
{
    // sort nodes
    int order[4] = {0, 1, 2, 3};
    sortOrder4(tMin, order[0], order[1], order[2], order[3]);

    // enqueue overlapping nodes in sorted order
    for (int w = 0; w < 4; w++) {
        int W = order[w];

        if (mask[W]) {
            stackPtr++;
            subtree[stackPtr].node = child[W];
            subtree[stackPtr].distance = tMin[W];
            tMaxMin = std::min(tMaxMin, tMax[W]);
        }
    }
}

template<size_t WIDTH, size_t DIM>
struct QueryStub {
    // empty
};

template<size_t WIDTH, size_t DIM, typename NodeType, typename LeafNodeType>
inline bool intersectRayPrimitives(QueryStub<WIDTH, DIM> queryStub, const NodeType& node,
                                   const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                   int aggregateIndex, const enokiVector<DIM>& ro, const enokiVector<DIM>& rd,
                                   float& rtMax, Interaction<DIM>& i, bool checkForOcclusion)
{
    std::cerr << "intersectRayPrimitives(): WIDTH: " << WIDTH << ", DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return 0;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline bool intersectRayPrimitives(QueryStub<WIDTH, 2> queryStub, const NodeType& node,
                                   const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                   int aggregateIndex, const enokiVector2& ro, const enokiVector2& rd,
                                   float& rtMax, Interaction<2>& i, bool checkForOcclusion)
{
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    bool didHit = false;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized intersection query
        FloatP<WIDTH> d;
        Vector2P<WIDTH> pt, n;
        FloatP<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector2P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector2P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        MaskP<WIDTH> mask = intersectWideLineSegment<WIDTH>(pa, pb, ro, rd, rtMax, d,
                                                            pt, n, t, checkForOcclusion);

        // determine closest index
        int closestIndex = -1;
        int W = std::min((int)WIDTH, nReferences - startReference);

        for (int w = 0; w < W; w++) {
            if (mask[w] && d[w] <= rtMax) {
                if (checkForOcclusion) return true;
                closestIndex = w;
                rtMax = d[w];
            }
        }

        // update interaction
        if (closestIndex != -1) {
            didHit = true;
            i.d = d[closestIndex];
            i.p[0] = pt[0][closestIndex];
            i.p[1] = pt[1][closestIndex];
            i.n[0] = n[0][closestIndex];
            i.n[1] = n[1][closestIndex];
            i.uv[0] = t[closestIndex];
            i.primitiveIndex = primitiveIndex[closestIndex];
            i.nodeIndex = nodeIndex;
            i.referenceIndex = referenceOffset + startReference + closestIndex;
            i.objectIndex = aggregateIndex;
        }

        startReference += WIDTH;
    }

    return didHit;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline bool intersectRayPrimitives(QueryStub<WIDTH, 3> queryStub, const NodeType& node,
                                   const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                   int aggregateIndex, const enokiVector3& ro, const enokiVector3& rd,
                                   float& rtMax, Interaction<3>& i, bool checkForOcclusion)
{
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    bool didHit = false;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized intersection query
        FloatP<WIDTH> d;
        Vector3P<WIDTH> pt, n;
        Vector2P<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector3P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector3P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const Vector3P<WIDTH>& pc = leafNodes[leafIndex].positions[2];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        MaskP<WIDTH> mask = intersectWideTriangle<WIDTH>(pa, pb, pc, ro, rd, rtMax, d,
                                                         pt, n, t, checkForOcclusion);

        // determine closest index
        int closestIndex = -1;
        int W = std::min((int)WIDTH, nReferences - startReference);

        for (int w = 0; w < W; w++) {
            if (mask[w] && d[w] <= rtMax) {
                if (checkForOcclusion) return true;
                closestIndex = w;
                rtMax = d[w];
            }
        }

        // update interaction
        if (closestIndex != -1) {
            didHit = true;
            i.d = d[closestIndex];
            i.p[0] = pt[0][closestIndex];
            i.p[1] = pt[1][closestIndex];
            i.p[2] = pt[2][closestIndex];
            i.n[0] = n[0][closestIndex];
            i.n[1] = n[1][closestIndex];
            i.n[2] = n[2][closestIndex];
            i.uv[0] = t[0][closestIndex];
            i.uv[1] = t[1][closestIndex];
            i.primitiveIndex = primitiveIndex[closestIndex];
            i.nodeIndex = nodeIndex;
            i.referenceIndex = referenceOffset + startReference + closestIndex;
            i.objectIndex = aggregateIndex;
        }

        startReference += WIDTH;
    }

    return didHit;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline bool Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                                                            int aggregateIndex, int& nodesVisited,
                                                            bool checkForOcclusion) const
{
    bool didHit = false;
    TraversalStack subtree[FCPW_MBVH_MAX_DEPTH];
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> tMin, tMax;
    enokiVector<DIM> ro = enoki::gather<enokiVector<DIM>>(r.o.data(), range);
    enokiVector<DIM> rd = enoki::gather<enokiVector<DIM>>(r.d.data(), range);
    enokiVector<DIM> rinvD = enoki::gather<enokiVector<DIM>>(r.invD.data(), range);
    QueryStub<WIDTH, DIM> queryStub;

    // push root node
    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    subtree[0].node = rootIndex;
    subtree[0].distance = minFloat;
    int stackPtr = 0;

    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found intersection, continue
        if (currentDist > r.tMax) continue;
        const NodeType& node(flatTree[nodeIndex]);

        if (isLeafNode(node)) {
            if (primitiveTypeSupportsVectorizedQueries) {
                // perform vectorized intersection query
                bool hit = intersectRayPrimitives(queryStub, node, leafNodes, nodeIndex, this->pIndex,
                                                  ro, rd, r.tMax, i, checkForOcclusion);

                nodesVisited++;
                if (hit) {
                    if (checkForOcclusion) return true;
                    didHit = true;
                }

            } else {
                // primitive type does not support vectorized intersection query,
                // perform query to each primitive one by one
                int referenceOffset = node.child[2];
                int nReferences = node.child[3];

                for (int p = 0; p < nReferences; p++) {
                    int referenceIndex = referenceOffset + p;
                    const PrimitiveType *prim = primitives[referenceIndex];
                    nodesVisited++;

                    bool hit = false;
                    if (primitiveTypeIsAggregate) {
                        const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(prim);
                        hit = aggregate->intersectFromNode(r, i, nodeStartIndex, aggregateIndex,
                                                           nodesVisited, checkForOcclusion);

                    } else {
                        const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(prim);
                        hit = geometricPrim->intersect(r, i, checkForOcclusion);
                    }

                    if (hit) {
                        if (checkForOcclusion) {
                            return true;
                        }

                        didHit = true;
                        r.tMax = std::min(r.tMax, i.d);
                        i.nodeIndex = nodeIndex;
                        i.referenceIndex = p;
                        i.objectIndex = this->pIndex;
                    }
                }
            }

        } else {
            // intersect ray with boxes
            MaskP<FCPW_MBVH_BRANCHING_FACTOR> mask = intersectWideBox<FCPW_MBVH_BRANCHING_FACTOR, DIM>(
                                                        node.boxMin, node.boxMax, ro, rinvD, r.tMax, tMin, tMax);

            // enqueue intersecting boxes in sorted order
            nodesVisited++;
            mask &= enoki::neq(node.child, maxInt);
            if (enoki::any(mask)) {
                float stub = 0.0f;
                enqueueNodes<FCPW_MBVH_BRANCHING_FACTOR>(node.child, tMin, tMax, mask,
                                                         r.tMax, stub, stackPtr, subtree);
            }
        }
    }

    return didHit;
}

template<size_t WIDTH, size_t DIM, typename NodeType, typename LeafNodeType>
inline int intersectRayPrimitives(QueryStub<WIDTH, DIM> queryStub, const NodeType& node,
                                  const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                  int aggregateIndex, const enokiVector<DIM>& ro, const enokiVector<DIM>& rd,
                                  float& rtMax, std::vector<Interaction<DIM>>& is, bool checkForOcclusion)
{
    std::cerr << "intersectRayPrimitives(): WIDTH: " << WIDTH << ", DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return 0;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline int intersectRayPrimitives(QueryStub<WIDTH, 2> queryStub, const NodeType& node,
                                  const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                  int aggregateIndex, const enokiVector2& ro, const enokiVector2& rd,
                                  float& rtMax, std::vector<Interaction<2>>& is, bool checkForOcclusion)
{
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    int hits = 0;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized intersection query
        FloatP<WIDTH> d;
        Vector2P<WIDTH> pt, n;
        FloatP<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector2P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector2P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        MaskP<WIDTH> mask = intersectWideLineSegment<WIDTH>(pa, pb, ro, rd, rtMax, d,
                                                            pt, n, t, checkForOcclusion);

        // record interactions
        int endReference = startReference + WIDTH;
        if (endReference > nReferences) endReference = nReferences;

        for (int p = startReference; p < endReference; p++) {
            int w = p - startReference;

            if (mask[w]) {
                hits++;
                auto it = is.emplace(is.end(), Interaction<2>());
                it->d = d[w];
                it->p[0] = pt[0][w];
                it->p[1] = pt[1][w];
                it->n[0] = n[0][w];
                it->n[1] = n[1][w];
                it->uv[0] = t[w];
                it->primitiveIndex = primitiveIndex[w];
                it->nodeIndex = nodeIndex;
                it->referenceIndex = referenceOffset + p;
                it->objectIndex = aggregateIndex;
            }
        }

        startReference += WIDTH;
    }

    return hits;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline int intersectRayPrimitives(QueryStub<WIDTH, 3> queryStub, const NodeType& node,
                                  const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                  int aggregateIndex, const enokiVector3& ro, const enokiVector3& rd,
                                  float& rtMax, std::vector<Interaction<3>>& is, bool checkForOcclusion)
{
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    int hits = 0;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized intersection query
        FloatP<WIDTH> d;
        Vector3P<WIDTH> pt, n;
        Vector2P<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector3P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector3P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const Vector3P<WIDTH>& pc = leafNodes[leafIndex].positions[2];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        MaskP<WIDTH> mask = intersectWideTriangle<WIDTH>(pa, pb, pc, ro, rd, rtMax, d,
                                                         pt, n, t, checkForOcclusion);

        // record interactions
        int endReference = startReference + WIDTH;
        if (endReference > nReferences) endReference = nReferences;

        for (int p = startReference; p < endReference; p++) {
            int w = p - startReference;

            if (mask[w]) {
                hits++;
                auto it = is.emplace(is.end(), Interaction<3>());
                it->d = d[w];
                it->p[0] = pt[0][w];
                it->p[1] = pt[1][w];
                it->p[2] = pt[2][w];
                it->n[0] = n[0][w];
                it->n[1] = n[1][w];
                it->n[2] = n[2][w];
                it->uv[0] = t[0][w];
                it->uv[1] = t[1][w];
                it->primitiveIndex = primitiveIndex[w];
                it->nodeIndex = nodeIndex;
                it->referenceIndex = referenceOffset + p;
                it->objectIndex = aggregateIndex;
            }
        }

        startReference += WIDTH;
    }

    return hits;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline int Mbvh<WIDTH, DIM,
                PrimitiveType,
                SilhouetteType,
                NodeType,
                LeafNodeType,
                SilhouetteLeafNodeType>::intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                                           int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                                                           bool checkForOcclusion, bool recordAllHits) const
{
    int hits = 0;
    if (!recordAllHits) is.resize(1);
    TraversalStack subtree[FCPW_MBVH_MAX_DEPTH];
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> tMin, tMax;
    enokiVector<DIM> ro = enoki::gather<enokiVector<DIM>>(r.o.data(), range);
    enokiVector<DIM> rd = enoki::gather<enokiVector<DIM>>(r.d.data(), range);
    enokiVector<DIM> rinvD = enoki::gather<enokiVector<DIM>>(r.invD.data(), range);
    QueryStub<WIDTH, DIM> queryStub;

    // push root node
    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    subtree[0].node = rootIndex;
    subtree[0].distance = minFloat;
    int stackPtr = 0;

    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found intersection, continue
        if (!recordAllHits && currentDist > r.tMax) continue;
        const NodeType& node(flatTree[nodeIndex]);

        if (isLeafNode(node)) {
            if (primitiveTypeSupportsVectorizedQueries) {
                // perform vectorized intersection query
                if (recordAllHits) {
                    hits += intersectRayPrimitives(queryStub, node, leafNodes, nodeIndex, this->pIndex,
                                                   ro, rd, r.tMax, is, checkForOcclusion);

                } else {
                    hits += intersectRayPrimitives(queryStub, node, leafNodes, nodeIndex, this->pIndex,
                                                   ro, rd, r.tMax, is[0], checkForOcclusion);
                }

                nodesVisited++;
                if (hits > 0 && checkForOcclusion) {
                    return 1;
                }

            } else {
                // primitive type does not support vectorized intersection query,
                // perform query to each primitive one by one
                int referenceOffset = node.child[2];
                int nReferences = node.child[3];

                for (int p = 0; p < nReferences; p++) {
                    int referenceIndex = referenceOffset + p;
                    const PrimitiveType *prim = primitives[referenceIndex];
                    nodesVisited++;

                    int hit = 0;
                    std::vector<Interaction<DIM>> cs;
                    if (primitiveTypeIsAggregate) {
                        const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(prim);
                        hit = aggregate->intersectFromNode(r, cs, nodeStartIndex, aggregateIndex,
                                                           nodesVisited, checkForOcclusion, recordAllHits);

                    } else {
                        const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(prim);
                        hit = geometricPrim->intersect(r, cs, checkForOcclusion, recordAllHits);
                        for (int i = 0; i < (int)cs.size(); i++) {
                            cs[i].nodeIndex = nodeIndex;
                            cs[i].referenceIndex = referenceIndex;
                            cs[i].objectIndex = this->pIndex;
                        }
                    }

                    // keep the closest intersection only
                    if (hit > 0) {
                        if (checkForOcclusion) {
                            return 1;
                        }

                        hits += hit;
                        if (recordAllHits) {
                            is.insert(is.end(), cs.begin(), cs.end());

                        } else {
                            r.tMax = std::min(r.tMax, cs[0].d);
                            is[0] = cs[0];
                        }
                    }
                }
            }

        } else {
            // intersect ray with boxes
            MaskP<FCPW_MBVH_BRANCHING_FACTOR> mask = intersectWideBox<FCPW_MBVH_BRANCHING_FACTOR, DIM>(
                                                        node.boxMin, node.boxMax, ro, rinvD, r.tMax, tMin, tMax);

            // enqueue intersecting boxes in sorted order
            nodesVisited++;
            mask &= enoki::neq(node.child, maxInt);
            if (enoki::any(mask)) {
                float stub = 0.0f;
                enqueueNodes<FCPW_MBVH_BRANCHING_FACTOR>(node.child, tMin, tMax, mask,
                                                         r.tMax, stub, stackPtr, subtree);
            }
        }
    }

    if (hits > 0) {
        // sort by distance and remove duplicates
        if (recordAllHits) {
            std::sort(is.begin(), is.end(), compareInteractions<DIM>);
            is = removeDuplicates<DIM>(is);
            hits = (int)is.size();

        } else {
            hits = 1;
        }

        return hits;
    }

    return 0;
}

template<size_t WIDTH, size_t DIM, typename NodeType, typename LeafNodeType>
inline int intersectSpherePrimitives(QueryStub<WIDTH, DIM> queryStub, const NodeType& node,
                                     const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                     int aggregateIndex, const enokiVector<DIM>& sc, float sr2,
                                     float u, Interaction<DIM>& i, float& totalPrimitiveWeight,
                                     bool isNodeInsideSphere=false)
{
    std::cerr << "intersectSpherePrimitives(): WIDTH: " << WIDTH << ", DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return 0;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline int intersectSpherePrimitives(QueryStub<WIDTH, 2> queryStub, const NodeType& node,
                                     const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                     int aggregateIndex, const enokiVector2& sc, float sr2,
                                     float u, Interaction<2>& i, float& totalPrimitiveWeight,
                                     bool isNodeInsideSphere=false)
{
    Vector2 queryPt(sc[0], sc[1]);
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    int hits = 0;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized closest point query
        Vector2P<WIDTH> pt;
        FloatP<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector2P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector2P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        FloatP<WIDTH> surfaceArea = enoki::norm(pb - pa);
        FloatP<WIDTH> d = isNodeInsideSphere ? 0.0f : findClosestPointWideLineSegment<WIDTH, 2>(pa, pb, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;

        // record interactions
        int endReference = startReference + WIDTH;
        if (endReference > nReferences) endReference = nReferences;

        for (int p = startReference; p < endReference; p++) {
            int w = p - startReference;

            if (d2[w] <= sr2) {
                hits++;
                float weight = surfaceArea[w];
                totalPrimitiveWeight += weight;
                float selectionProb = weight/totalPrimitiveWeight;

                if (u < selectionProb) {
                    u = u/selectionProb; // rescale to [0,1)
                    i.d = weight;
                    i.primitiveIndex = primitiveIndex[w];
                    i.nodeIndex = nodeIndex;
                    i.referenceIndex = referenceOffset + p;
                    i.objectIndex = aggregateIndex;

                } else {
                    u = (u - selectionProb)/(1.0f - selectionProb);
                }
            }
        }

        startReference += WIDTH;
    }

    return hits;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline int intersectSpherePrimitives(QueryStub<WIDTH, 3> queryStub, const NodeType& node,
                                     const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                     int aggregateIndex, const enokiVector3& sc, float sr2,
                                     float u, Interaction<3>& i, float& totalPrimitiveWeight,
                                     bool isNodeInsideSphere=false)
{
    Vector3 queryPt(sc[0], sc[1], sc[2]);
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    int hits = 0;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized closest point query
        Vector3P<WIDTH> pt;
        Vector2P<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector3P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector3P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const Vector3P<WIDTH>& pc = leafNodes[leafIndex].positions[2];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        FloatP<WIDTH> surfaceArea = 0.5f*enoki::norm(enoki::cross(pb - pa, pc - pa));
        FloatP<WIDTH> d = isNodeInsideSphere ? 0.0f : findClosestPointWideTriangle<WIDTH>(pa, pb, pc, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;

        // record interactions
        int endReference = startReference + WIDTH;
        if (endReference > nReferences) endReference = nReferences;

        for (int p = startReference; p < endReference; p++) {
            int w = p - startReference;

            if (d2[w] <= sr2) {
                hits++;
                float weight = surfaceArea[w];
                totalPrimitiveWeight += weight;
                float selectionProb = weight/totalPrimitiveWeight;

                if (u < selectionProb) {
                    u = u/selectionProb; // rescale to [0,1)
                    i.d = weight;
                    i.primitiveIndex = primitiveIndex[w];
                    i.nodeIndex = nodeIndex;
                    i.referenceIndex = referenceOffset + p;
                    i.objectIndex = aggregateIndex;

                } else {
                    u = (u - selectionProb)/(1.0f - selectionProb);
                }
            }
        }

        startReference += WIDTH;
    }

    return hits;
}

template<size_t WIDTH, size_t DIM, typename NodeType, typename LeafNodeType>
inline int intersectSpherePrimitives(QueryStub<WIDTH, DIM> queryStub, const NodeType& node,
                                     const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                     int aggregateIndex, const enokiVector<DIM>& sc, float sr2,
                                     std::vector<Interaction<DIM>>& is, float& totalPrimitiveWeight,
                                     bool isNodeInsideSphere=false)
{
    std::cerr << "intersectSpherePrimitives(): WIDTH: " << WIDTH << ", DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return 0;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline int intersectSpherePrimitives(QueryStub<WIDTH, 2> queryStub, const NodeType& node,
                                     const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                     int aggregateIndex, const enokiVector2& sc, float sr2,
                                     std::vector<Interaction<2>>& is, float& totalPrimitiveWeight,
                                     bool isNodeInsideSphere=false)
{
    Vector2 queryPt(sc[0], sc[1]);
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    int hits = 0;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized closest point query
        Vector2P<WIDTH> pt;
        FloatP<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector2P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector2P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        FloatP<WIDTH> surfaceArea = enoki::norm(pb - pa);
        FloatP<WIDTH> d = isNodeInsideSphere ? 0.0f : findClosestPointWideLineSegment<WIDTH, 2>(pa, pb, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;

        // record interactions
        int endReference = startReference + WIDTH;
        if (endReference > nReferences) endReference = nReferences;

        for (int p = startReference; p < endReference; p++) {
            int w = p - startReference;

            if (d2[w] <= sr2) {
                hits++;
                auto it = is.emplace(is.end(), Interaction<2>());
                it->d = 1.0f;
                it->primitiveIndex = primitiveIndex[w];
                it->nodeIndex = nodeIndex;
                it->referenceIndex = referenceOffset + p;
                it->objectIndex = aggregateIndex;
            }
        }

        startReference += WIDTH;
    }

    return hits;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline int intersectSpherePrimitives(QueryStub<WIDTH, 3> queryStub, const NodeType& node,
                                     const std::vector<LeafNodeType>& leafNodes, int nodeIndex,
                                     int aggregateIndex, const enokiVector3& sc, float sr2,
                                     std::vector<Interaction<3>>& is, float& totalPrimitiveWeight,
                                     bool isNodeInsideSphere=false)
{
    Vector3 queryPt(sc[0], sc[1], sc[2]);
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    int hits = 0;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized closest point query
        Vector3P<WIDTH> pt;
        Vector2P<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector3P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector3P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const Vector3P<WIDTH>& pc = leafNodes[leafIndex].positions[2];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        FloatP<WIDTH> surfaceArea = 0.5f*enoki::norm(enoki::cross(pb - pa, pc - pa));
        FloatP<WIDTH> d = isNodeInsideSphere ? 0.0f : findClosestPointWideTriangle<WIDTH>(pa, pb, pc, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;

        // record interactions
        int endReference = startReference + WIDTH;
        if (endReference > nReferences) endReference = nReferences;

        for (int p = startReference; p < endReference; p++) {
            int w = p - startReference;

            if (d2[w] <= sr2) {
                hits++;
                auto it = is.emplace(is.end(), Interaction<3>());
                it->d = 1.0f;
                it->primitiveIndex = primitiveIndex[w];
                it->nodeIndex = nodeIndex;
                it->referenceIndex = referenceOffset + p;
                it->objectIndex = aggregateIndex;
            }
        }

        startReference += WIDTH;
    }

    return hits;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline int Mbvh<WIDTH, DIM,
                PrimitiveType,
                SilhouetteType,
                NodeType,
                LeafNodeType,
                SilhouetteLeafNodeType>::intersectFromNode(const BoundingSphere<DIM>& s,
                                                           std::vector<Interaction<DIM>>& is,
                                                           int nodeStartIndex, int aggregateIndex,
                                                           int& nodesVisited, bool recordOneHit) const
{
    int hits = 0;
    float totalPrimitiveWeight = 0.0f;
    if (recordOneHit && !primitiveTypeIsAggregate) is.resize(1);
    TraversalStack subtree[FCPW_MBVH_MAX_DEPTH];
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> d2Min, d2Max;
    enokiVector<DIM> sc = enoki::gather<enokiVector<DIM>>(s.c.data(), range);
    QueryStub<WIDTH, DIM> queryStub;

    // push root node
    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    subtree[0].node = rootIndex;
    subtree[0].distance = s.r2;
    int stackPtr = 0;

    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        const NodeType& node(flatTree[nodeIndex]);
        stackPtr--;

        if (isLeafNode(node)) {
            if (primitiveTypeSupportsVectorizedQueries) {
                // perform vectorized intersection query
                if (recordOneHit) {
                    float u = uniformRealRandomNumber();
                    hits += intersectSpherePrimitives(queryStub, node, leafNodes, nodeIndex, this->pIndex,
                                                      sc, s.r2, u, is[0], totalPrimitiveWeight);

                } else {
                    hits += intersectSpherePrimitives(queryStub, node, leafNodes, nodeIndex, this->pIndex,
                                                      sc, s.r2, is, totalPrimitiveWeight);
                }

                nodesVisited++;

            } else {
                // primitive type does not support vectorized intersection query,
                // perform query to each primitive one by one
                int referenceOffset = node.child[2];
                int nReferences = node.child[3];

                for (int p = 0; p < nReferences; p++) {
                    int referenceIndex = referenceOffset + p;
                    const PrimitiveType *prim = primitives[referenceIndex];
                    nodesVisited++;

                    if (primitiveTypeIsAggregate) {
                        std::vector<Interaction<DIM>> cs;
                        const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(prim);
                        hits += aggregate->intersectFromNode(s, cs, nodeStartIndex, aggregateIndex,
                                                             nodesVisited, recordOneHit);
                        is.insert(is.end(), cs.begin(), cs.end());

                    } else {
                        Interaction<DIM> c;
                        const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(prim);
                        bool hit = geometricPrim->intersect(s, c, recordOneHit);

                        if (hit) {
                            hits += 1;
                            c.nodeIndex = nodeIndex;
                            c.referenceIndex = referenceIndex;
                            c.objectIndex = this->pIndex;

                            if (recordOneHit) {
                                totalPrimitiveWeight += c.d;
                                if (uniformRealRandomNumber()*totalPrimitiveWeight < c.d) {
                                    is[0] = c;
                                }

                            } else {
                                is.emplace_back(c);
                            }
                        }
                    }
                }
            }

        } else {
            // overlap sphere with boxes
            MaskP<FCPW_MBVH_BRANCHING_FACTOR> mask = overlapWideBox<FCPW_MBVH_BRANCHING_FACTOR, DIM>(
                                                        node.boxMin, node.boxMax, sc, s.r2, d2Min, d2Max);

            // enqueue overlapping boxes
            nodesVisited++;
            mask &= enoki::neq(node.child, maxInt);
            if (enoki::any(mask)) {
                for (int w = 0; w < FCPW_MBVH_BRANCHING_FACTOR; w++) {
                    if (mask[w]) {
                        stackPtr++;
                        subtree[stackPtr].node = node.child[w];
                    }
                }
            }
        }
    }

    if (hits > 0) {
        if (recordOneHit && !primitiveTypeIsAggregate) {
            if (is[0].primitiveIndex == -1) {
                hits = 0;
                is.clear();

            } else if (totalPrimitiveWeight > 0.0f) {
                is[0].d /= totalPrimitiveWeight;
            }
        }

        return hits;
    }

    return 0;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline int Mbvh<WIDTH, DIM,
                PrimitiveType,
                SilhouetteType,
                NodeType,
                LeafNodeType,
                SilhouetteLeafNodeType>::intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                           const Vector<DIM>& randNums, int nodeStartIndex,
                                                           int aggregateIndex, int& nodesVisited,
                                                           const std::function<float(float)>& branchTraversalWeight) const
{
    int hits = 0;
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> d2Min, d2Max;
    float d2NodeMax = maxFloat;
    float u = randNums[0];
    enokiVector<DIM> sc = enoki::gather<enokiVector<DIM>>(s.c.data(), range);
    QueryStub<WIDTH, DIM> queryStub;

    // push root node
    int nodeIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    float traversalPdf = 1.0f;
    int stackPtr = 0;

    while (stackPtr >= 0) {
        // pop off the next node to work on
        const NodeType& node(flatTree[nodeIndex]);
        stackPtr--;

        if (isLeafNode(node)) {
            float totalPrimitiveWeight = 0.0f;
            if (primitiveTypeSupportsVectorizedQueries) {
                // perform vectorized intersection query
                int hit = intersectSpherePrimitives(queryStub, node, leafNodes, nodeIndex, this->pIndex,
                                                    sc, s.r2, u, i, totalPrimitiveWeight, d2NodeMax <= s.r2);
                nodesVisited++;

                if (hit > 0) {
                    hits += hit;
                    i.d *= traversalPdf;
                }

            } else {
                // primitive type does not support vectorized intersection query,
                // perform query to each primitive one by one
                int referenceOffset = node.child[2];
                int nReferences = node.child[3];

                for (int p = 0; p < nReferences; p++) {
                    int referenceIndex = referenceOffset + p;
                    const PrimitiveType *prim = primitives[referenceIndex];
                    nodesVisited++;

                    int hit = 0;
                    Interaction<DIM> c;
                    if (primitiveTypeIsAggregate) {
                        Vector<DIM> modifiedRandNums;
                        modifiedRandNums[0] = u;
                        for (int i = 1; i < DIM; i++) modifiedRandNums[i] = randNums[i];
                        const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(prim);
                        hit = aggregate->intersectFromNode(s, c, modifiedRandNums, nodeStartIndex, aggregateIndex,
                                                           nodesVisited, branchTraversalWeight);

                    } else {
                        if (d2NodeMax <= s.r2) {
                            hit = 1;
                            c.primitiveIndex = prim->getIndex();
                            c.d = prim->surfaceArea();

                        } else {
                            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(prim);
                            hit = geometricPrim->intersect(s, c, true) ? 1 : 0;
                        }

                        c.nodeIndex = nodeIndex;
                        c.referenceIndex = referenceIndex;
                        c.objectIndex = this->pIndex;
                    }

                    if (hit > 0) {
                        hits += hit;
                        float selectionProb = 0.0f;
                        if (primitiveTypeIsAggregate) {
                            // choose sample uniformly amongst aggregates
                            totalPrimitiveWeight += 1.0f;
                            selectionProb = 1.0f/totalPrimitiveWeight;

                        } else {
                            // choose sample in proportion to surface area amongst primitives
                            totalPrimitiveWeight += c.d;
                            selectionProb = c.d/totalPrimitiveWeight;
                        }

                        if (u < selectionProb) {
                            u = u/selectionProb; // rescale to [0,1)
                            i = c;
                            i.d *= traversalPdf;

                        } else {
                            u = (u - selectionProb)/(1.0f - selectionProb);
                        }
                    }
                }
            }

            if (totalPrimitiveWeight > 0.0f) {
                i.d /= totalPrimitiveWeight;
            }

        } else {
            // overlap sphere with boxes
            MaskP<FCPW_MBVH_BRANCHING_FACTOR> mask = overlapWideBox<FCPW_MBVH_BRANCHING_FACTOR, DIM>(
                                                        node.boxMin, node.boxMax, sc, s.r2, d2Min, d2Max);

            // enqueue overlapping boxes
            nodesVisited++;
            mask &= enoki::neq(node.child, maxInt);
            if (enoki::any(mask)) {
                int selectedIndex = -1;
                float selectedWeight = 0.0f;
                float totalTraversalWeight = 0.0f;
                FloatP<FCPW_MBVH_BRANCHING_FACTOR> r2;
                if (branchTraversalWeight) {
                    VectorP<FCPW_MBVH_BRANCHING_FACTOR, DIM> boxCenter = (node.boxMin + node.boxMax)*0.5f;
                    r2 = enoki::squared_norm(sc - boxCenter);
                }

                for (int w = 0; w < FCPW_MBVH_BRANCHING_FACTOR; w++) {
                    if (mask[w]) {
                        float weight = branchTraversalWeight ? branchTraversalWeight(r2[w]) : 1.0f;
                        totalTraversalWeight += weight;
                        float prob = weight/totalTraversalWeight;

                        if (u < prob) {
                            selectedIndex = w;
                            selectedWeight = weight;
                            u = u/prob;

                        } else {
                            u = (u - prob)/(1.0f - prob);
                        }
                    }
                }

                if (selectedIndex != -1) {
                    stackPtr++;
                    nodeIndex = node.child[selectedIndex];
                    traversalPdf *= selectedWeight/totalTraversalWeight;
                    d2NodeMax = d2Max[selectedIndex];
                }
            }
        }
    }

    if (hits > 0) {
        if (!primitiveTypeIsAggregate) {
            if (i.primitiveIndex == -1) {
                hits = 0;

            } else {
                // sample a point on the selected geometric primitive
                const PrimitiveType *prim = primitives[i.referenceIndex];
                float pdf = i.samplePoint(prim, randNums);
                i.d *= pdf;
            }
        }

        return hits;
    }

    return 0;
}

template<size_t WIDTH, size_t DIM, typename NodeType, typename LeafNodeType>
inline bool findClosestPointPrimitives(QueryStub<WIDTH, DIM> queryStub, const NodeType& node,
                                       const std::vector<LeafNodeType>& leafNodes,
                                       int nodeIndex, int aggregateIndex, const enokiVector<DIM>& sc,
                                       float& sr2, Interaction<DIM>& i)
{
    std::cerr << "findClosestPointPrimitives(): WIDTH: " << WIDTH << ", DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return false;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline bool findClosestPointPrimitives(QueryStub<WIDTH, 2> queryStub, const NodeType& node,
                                       const std::vector<LeafNodeType>& leafNodes,
                                       int nodeIndex, int aggregateIndex, const enokiVector2& sc,
                                       float& sr2, Interaction<2>& i)
{
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    bool found = false;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized closest point query
        Vector2P<WIDTH> pt;
        FloatP<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector2P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector2P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        FloatP<WIDTH> d = findClosestPointWideLineSegment<WIDTH, 2>(pa, pb, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;

        // determine closest index
        int closestIndex = -1;
        int W = std::min((int)WIDTH, nReferences - startReference);

        for (int w = 0; w < W; w++) {
            if (d2[w] <= sr2) {
                closestIndex = w;
                sr2 = d2[w];
            }
        }

        // update interaction
        if (closestIndex != -1) {
            i.d = d[closestIndex];
            i.p[0] = pt[0][closestIndex];
            i.p[1] = pt[1][closestIndex];
            i.uv[0] = t[closestIndex];
            i.primitiveIndex = primitiveIndex[closestIndex];
            i.nodeIndex = nodeIndex;
            i.referenceIndex = referenceOffset + startReference + closestIndex;
            i.objectIndex = aggregateIndex;
            found = true;
        }

        startReference += WIDTH;
    }

    return found;
}

template<size_t WIDTH, typename NodeType, typename LeafNodeType>
inline bool findClosestPointPrimitives(QueryStub<WIDTH, 3> queryStub, const NodeType& node,
                                       const std::vector<LeafNodeType>& leafNodes,
                                       int nodeIndex, int aggregateIndex, const enokiVector3& sc,
                                       float& sr2, Interaction<3>& i)
{
    int leafOffset = -node.child[0] - 1;
    int nLeafs = node.child[1];
    int referenceOffset = node.child[2];
    int nReferences = node.child[3];
    int startReference = 0;
    bool found = false;

    for (int l = 0; l < nLeafs; l++) {
        // perform vectorized closest point query
        Vector3P<WIDTH> pt;
        Vector2P<WIDTH> t;
        int leafIndex = leafOffset + l;
        const Vector3P<WIDTH>& pa = leafNodes[leafIndex].positions[0];
        const Vector3P<WIDTH>& pb = leafNodes[leafIndex].positions[1];
        const Vector3P<WIDTH>& pc = leafNodes[leafIndex].positions[2];
        const IntP<WIDTH>& primitiveIndex = leafNodes[leafIndex].primitiveIndex;
        FloatP<WIDTH> d = findClosestPointWideTriangle<WIDTH>(pa, pb, pc, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;

        // determine closest index
        int closestIndex = -1;
        int W = std::min((int)WIDTH, nReferences - startReference);

        for (int w = 0; w < W; w++) {
            if (d2[w] <= sr2) {
                closestIndex = w;
                sr2 = d2[w];
            }
        }

        // update interaction
        if (closestIndex != -1) {
            i.d = d[closestIndex];
            i.p[0] = pt[0][closestIndex];
            i.p[1] = pt[1][closestIndex];
            i.p[2] = pt[2][closestIndex];
            i.uv[0] = t[0][closestIndex];
            i.uv[1] = t[1][closestIndex];
            i.primitiveIndex = primitiveIndex[closestIndex];
            i.nodeIndex = nodeIndex;
            i.referenceIndex = referenceOffset + startReference + closestIndex;
            i.objectIndex = aggregateIndex;
            found = true;
        }

        startReference += WIDTH;
    }

    return found;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline bool Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                   int nodeStartIndex, int aggregateIndex,
                                                                   int& nodesVisited, bool recordNormal) const
{
    bool notFound = true;
    TraversalStack subtree[FCPW_MBVH_MAX_DEPTH];
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> d2Min, d2Max;
    enokiVector<DIM> sc = enoki::gather<enokiVector<DIM>>(s.c.data(), range);
    QueryStub<WIDTH, DIM> queryStub;

    // push root node
    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    subtree[0].node = rootIndex;
    subtree[0].distance = minFloat;
    int stackPtr = 0;

    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found primitive, continue
        if (currentDist > s.r2) continue;
        const NodeType& node(flatTree[nodeIndex]);

        if (isLeafNode(node)) {
            if (primitiveTypeSupportsVectorizedQueries) {
                // perform vectorized closest point query
                bool found = findClosestPointPrimitives(queryStub, node, leafNodes, nodeIndex,
                                                        this->pIndex, sc, s.r2, i);
                if (found) notFound = false;
                nodesVisited++;

            } else {
                // primitive type does not support vectorized closest point query,
                // perform query to each primitive one by one
                int referenceOffset = node.child[2];
                int nReferences = node.child[3];

                for (int p = 0; p < nReferences; p++) {
                    int referenceIndex = referenceOffset + p;
                    const PrimitiveType *prim = primitives[referenceIndex];
                    nodesVisited++;

                    bool found = false;
                    Interaction<DIM> c;

                    if (primitiveTypeIsAggregate) {
                        const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(prim);
                        found = aggregate->findClosestPointFromNode(s, c, nodeStartIndex, aggregateIndex,
                                                                    nodesVisited, recordNormal);

                    } else {
                        const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(prim);
                        found = geometricPrim->findClosestPoint(s, c);
                        c.nodeIndex = nodeIndex;
                        c.referenceIndex = referenceIndex;
                        c.objectIndex = this->pIndex;
                    }

                    // keep the closest point only
                    if (found) {
                        notFound = false;
                        s.r2 = std::min(s.r2, c.d*c.d);
                        i = c;
                    }
                }
            }

        } else {
            // overlap sphere with boxes
            MaskP<FCPW_MBVH_BRANCHING_FACTOR> mask = overlapWideBox<FCPW_MBVH_BRANCHING_FACTOR, DIM>(
                                                                node.boxMin, node.boxMax, sc, s.r2, d2Min, d2Max);

            // enqueue overlapping boxes in sorted order
            nodesVisited++;
            mask &= enoki::neq(node.child, maxInt);
            if (enoki::any(mask)) {
                enqueueNodes<FCPW_MBVH_BRANCHING_FACTOR>(node.child, d2Min, d2Max, mask,
                                                         s.r2, s.r2, stackPtr, subtree);
            }
        }
    }

    if (!notFound) {
        // compute normal
        if (recordNormal && !primitiveTypeIsAggregate) {
            i.computeNormal(primitives[i.referenceIndex]);
        }

        return true;
    }

    return false;
}

template<size_t WIDTH, size_t DIM, typename SilhouetteLeafNodeType>
inline bool findClosestSilhouettes(QueryStub<WIDTH, DIM> queryStub, const MsnchNode<DIM>& node,
                                   const std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                   int nodeIndex, int aggregateIndex, const enokiVector<DIM>& sc,
                                   float& sr2, Interaction<DIM>& i, bool flipNormalOrientation,
                                   float squaredMinRadius, float precision)
{
    std::cerr << "findClosestSilhouettes(): WIDTH: " << WIDTH << ", DIM: " << DIM << " not supported" << std::endl;
    exit(EXIT_FAILURE);

    return false;
}

template<size_t WIDTH, typename SilhouetteLeafNodeType>
inline bool findClosestSilhouettes(QueryStub<WIDTH, 2> queryStub, const MsnchNode<2>& node,
                                   const std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                   int nodeIndex, int aggregateIndex, const enokiVector2& sc,
                                   float& sr2, Interaction<2>& i, bool flipNormalOrientation,
                                   float squaredMinRadius, float precision)
{
    if (squaredMinRadius >= sr2) return false;

    int silhouetteLeafOffset = -node.silhouetteChild[0] - 1;
    int nSilhouetteLeafs = node.silhouetteChild[1];
    int silhouetteReferenceOffset = node.silhouetteChild[2];
    int nSilhouetteReferences = node.silhouetteChild[3];
    int startReference = 0;
    bool found = false;

    for (int l = 0; l < nSilhouetteLeafs; l++) {
        // perform vectorized closest silhouette query
        int leafIndex = silhouetteLeafOffset + l;
        const Vector2P<WIDTH>& pb = silhouetteLeafNodes[leafIndex].positions[1];
        const Vector2P<WIDTH>& n0 = silhouetteLeafNodes[leafIndex].positions[0];
        const Vector2P<WIDTH>& n1 = silhouetteLeafNodes[leafIndex].positions[2];
        const IntP<WIDTH>& primitiveIndex = silhouetteLeafNodes[leafIndex].primitiveIndex;
        const MaskP<WIDTH>& missingFace = silhouetteLeafNodes[leafIndex].missingFace;
        Vector2P<WIDTH> viewDir = sc - pb;
        FloatP<WIDTH> d = enoki::norm(viewDir);
        FloatP<WIDTH> d2 = d*d;
        if (enoki::all(d2 > sr2)) continue;
        MaskP<WIDTH> isSilhouette = missingFace;
        enoki::masked(isSilhouette, ~missingFace) = isWideSilhouetteVertex(n0, n1, viewDir, d, flipNormalOrientation, precision);

        // determine closest index
        int closestIndex = -1;
        int W = std::min((int)WIDTH, nSilhouetteReferences - startReference);

        for (int w = 0; w < W; w++) {
            if (isSilhouette[w] && d2[w] <= sr2) {
                closestIndex = w;
                sr2 = d2[w];
            }
        }

        // update interaction
        if (closestIndex != -1) {
            i.d = d[closestIndex];
            i.p[0] = pb[0][closestIndex];
            i.p[1] = pb[1][closestIndex];
            i.uv[0] = -1;
            i.primitiveIndex = primitiveIndex[closestIndex];
            i.nodeIndex = nodeIndex;
            i.referenceIndex = silhouetteReferenceOffset + startReference + closestIndex;
            i.objectIndex = aggregateIndex;
            found = true;
        }

        startReference += WIDTH;
    }

    return found;
}

template<size_t WIDTH, typename SilhouetteLeafNodeType>
inline bool findClosestSilhouettes(QueryStub<WIDTH, 3> queryStub, const MsnchNode<3>& node,
                                   const std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                   int nodeIndex, int aggregateIndex, const enokiVector3& sc,
                                   float& sr2, Interaction<3>& i, bool flipNormalOrientation,
                                   float squaredMinRadius, float precision)
{
    if (squaredMinRadius >= sr2) return false;

    int silhouetteLeafOffset = -node.silhouetteChild[0] - 1;
    int nSilhouetteLeafs = node.silhouetteChild[1];
    int silhouetteReferenceOffset = node.silhouetteChild[2];
    int nSilhouetteReferences = node.silhouetteChild[3];
    int startReference = 0;
    bool found = false;

    for (int l = 0; l < nSilhouetteLeafs; l++) {
        // perform vectorized closest silhouette query
        Vector3P<WIDTH> pt;
        FloatP<WIDTH> t;
        int leafIndex = silhouetteLeafOffset + l;
        const Vector3P<WIDTH>& pb = silhouetteLeafNodes[leafIndex].positions[1];
        const Vector3P<WIDTH>& pc = silhouetteLeafNodes[leafIndex].positions[2];
        const Vector3P<WIDTH>& n0 = silhouetteLeafNodes[leafIndex].positions[0];
        const Vector3P<WIDTH>& n1 = silhouetteLeafNodes[leafIndex].positions[3];
        const IntP<WIDTH>& primitiveIndex = silhouetteLeafNodes[leafIndex].primitiveIndex;
        const MaskP<WIDTH>& missingFace = silhouetteLeafNodes[leafIndex].missingFace;
        FloatP<WIDTH> d = findClosestPointWideLineSegment<WIDTH, 3>(pb, pc, sc, pt, t);
        FloatP<WIDTH> d2 = d*d;
        if (enoki::all(d2 > sr2)) continue;
        Vector3P<WIDTH> viewDir = sc - pt;
        MaskP<WIDTH> isSilhouette = missingFace;
        enoki::masked(isSilhouette, ~missingFace) = isWideSilhouetteEdge(pb, pc, n0, n1, viewDir, d, flipNormalOrientation, precision);

        // determine closest index
        int closestIndex = -1;
        int W = std::min((int)WIDTH, nSilhouetteReferences - startReference);

        for (int w = 0; w < W; w++) {
            if (isSilhouette[w] && d2[w] <= sr2) {
                closestIndex = w;
                sr2 = d2[w];
            }
        }

        // update interaction
        if (closestIndex != -1) {
            i.d = d[closestIndex];
            i.p[0] = pt[0][closestIndex];
            i.p[1] = pt[1][closestIndex];
            i.p[2] = pt[2][closestIndex];
            i.uv[0] = t[closestIndex];
            i.uv[1] = -1;
            i.primitiveIndex = primitiveIndex[closestIndex];
            i.nodeIndex = nodeIndex;
            i.referenceIndex = silhouetteReferenceOffset + startReference + closestIndex;
            i.objectIndex = aggregateIndex;
            found = true;
        }

        startReference += WIDTH;
    }

    return found;
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename SilhouetteLeafNodeType>
inline void processSubtreeForClosestSilhouettePoint(QueryStub<WIDTH, DIM> queryStub,
                                                    const std::vector<NodeType>& flatTree,
                                                    const std::vector<PrimitiveType *>& primitives,
                                                    const std::vector<SilhouetteType *>& silhouetteRefs,
                                                    const std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                                    const enokiVector<DIM>& sc, BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                    int nodeStartIndex, int aggregateIndex, int objectIndex,
                                                    bool primitiveTypeIsAggregate, bool silhouetteTypeSupportsVectorizedQueries,
                                                    bool flipNormalOrientation, float squaredMinRadius, float precision,
                                                    bool recordNormal, TraversalStack *subtree,
                                                    FloatP<FCPW_MBVH_BRANCHING_FACTOR>& d2Min,
                                                    bool& notFound, int& nodesVisited)
{
    std::cerr << "Mbvh::processSubtreeForClosestSilhouettePoint() not implemented" << std::endl;
    exit(EXIT_FAILURE);
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename SilhouetteLeafNodeType>
inline void processSubtreeForClosestSilhouettePoint(QueryStub<WIDTH, DIM> queryStub,
                                                    const std::vector<MsnchNode<DIM>>& flatTree,
                                                    const std::vector<PrimitiveType *>& primitives,
                                                    const std::vector<SilhouetteType *>& silhouetteRefs,
                                                    const std::vector<SilhouetteLeafNodeType>& silhouetteLeafNodes,
                                                    const enokiVector<DIM>& sc, BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                    int nodeStartIndex, int aggregateIndex, int objectIndex,
                                                    bool primitiveTypeIsAggregate, bool silhouetteTypeSupportsVectorizedQueries,
                                                    bool flipNormalOrientation, float squaredMinRadius, float precision,
                                                    bool recordNormal, TraversalStack *subtree,
                                                    FloatP<FCPW_MBVH_BRANCHING_FACTOR>& d2Min,
                                                    bool& notFound, int& nodesVisited)
{
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> stubs[2];
    int stackPtr = 0;
    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found primitive, continue
        if (currentDist > s.r2) continue;
        const MsnchNode<DIM>& node(flatTree[nodeIndex]);

        if (node.child[0] < 0) { // is leaf
            if (primitiveTypeIsAggregate) {
                int referenceOffset = node.child[2];
                int nReferences = node.child[3];

                for (int p = 0; p < nReferences; p++) {
                    int referenceIndex = referenceOffset + p;
                    const PrimitiveType *prim = primitives[referenceIndex];
                    nodesVisited++;

                    Interaction<DIM> c;
                    const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(prim);
                    bool found = aggregate->findClosestSilhouettePointFromNode(s, c, nodeStartIndex, aggregateIndex,
                                                                               nodesVisited, flipNormalOrientation,
                                                                               squaredMinRadius, precision, recordNormal);

                    // keep the closest silhouette point
                    if (found) {
                        notFound = false;
                        s.r2 = std::min(s.r2, c.d*c.d);
                        i = c;

                        if (squaredMinRadius >= s.r2) {
                            break;
                        }
                    }
                }

            } else {
                if (silhouetteTypeSupportsVectorizedQueries) {
                    // perform vectorized closest silhouette query
                    nodesVisited++;
                    bool found = findClosestSilhouettes(queryStub, node, silhouetteLeafNodes, nodeIndex, objectIndex,
                                                        sc, s.r2, i, flipNormalOrientation, squaredMinRadius, precision);
                    if (found) {
                        notFound = false;
                        if (squaredMinRadius >= s.r2) break;
                    }

                } else {
                    // silhouette type does not support vectorized closest silhouette
                    // query, perform query to each silhouette one by one
                    int silhouetteReferenceOffset = node.silhouetteChild[2];
                    int nSilhouetteReferences = node.silhouetteChild[3];

                    for (int p = 0; p < nSilhouetteReferences; p++) {
                        int referenceIndex = silhouetteReferenceOffset + p;
                        const SilhouetteType *silhouette = silhouetteRefs[referenceIndex];

                        // skip query if silhouette index is the same as i.primitiveIndex (and object indices match)
                        if (silhouette->getIndex() == i.primitiveIndex && objectIndex == i.objectIndex) continue;
                        nodesVisited++;

                        Interaction<DIM> c;
                        bool found = silhouette->findClosestSilhouettePoint(s, c, flipNormalOrientation,
                                                                            squaredMinRadius, precision);

                        // keep the closest silhouette point
                        if (found) {
                            notFound = false;
                            s.r2 = std::min(s.r2, c.d*c.d);
                            i = c;
                            i.nodeIndex = nodeIndex;
                            i.referenceIndex = referenceIndex;
                            i.objectIndex = objectIndex;

                            if (squaredMinRadius >= s.r2) {
                                break;
                            }
                        }
                    }
                }
            }

        } else { // not a leaf
            // overlap sphere with boxes, and normal and view cones
            MaskP<FCPW_MBVH_BRANCHING_FACTOR> mask = enoki::neq(node.child, maxInt) && node.coneHalfAngle >= 0.0f;
            mask &= overlapWideBox<FCPW_MBVH_BRANCHING_FACTOR, DIM>(node.boxMin, node.boxMax, sc, s.r2, d2Min);
            overlapWideCone<FCPW_MBVH_BRANCHING_FACTOR, DIM>(node.coneAxis, node.coneHalfAngle, node.coneRadius, sc,
                                                             node.boxMin, node.boxMax, d2Min, stubs[0], stubs[1], mask);

            // enqueue overlapping boxes in sorted order
            nodesVisited++;
            if (enoki::any(mask)) {
                float stub = 0.0f;
                enqueueNodes<FCPW_MBVH_BRANCHING_FACTOR>(node.child, d2Min, 0.0f, mask,
                                                         s.r2, stub, stackPtr, subtree);
            }
        }
    }
}

template<size_t WIDTH, size_t DIM,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline bool Mbvh<WIDTH, DIM,
                 PrimitiveType,
                 SilhouetteType,
                 NodeType,
                 LeafNodeType,
                 SilhouetteLeafNodeType>::findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                             int nodeStartIndex, int aggregateIndex,
                                                                             int& nodesVisited, bool flipNormalOrientation,
                                                                             float squaredMinRadius, float precision,
                                                                             bool recordNormal) const
{
    if (squaredMinRadius >= s.r2) return false;

    bool notFound = true;
    FloatP<FCPW_MBVH_BRANCHING_FACTOR> d2Min;
    enokiVector<DIM> sc = enoki::gather<enokiVector<DIM>>(s.c.data(), range);
    TraversalStack subtree[FCPW_MBVH_MAX_DEPTH];
    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    subtree[0].node = rootIndex;
    subtree[0].distance = minFloat;
    QueryStub<WIDTH, DIM> queryStub;

    processSubtreeForClosestSilhouettePoint(queryStub, flatTree, primitives, silhouetteRefs, silhouetteLeafNodes,
                                            sc, s, i, nodeStartIndex, aggregateIndex, this->pIndex,
                                            primitiveTypeIsAggregate, silhouetteTypeSupportsVectorizedQueries,
                                            flipNormalOrientation, squaredMinRadius, precision, recordNormal,
                                            subtree, d2Min, notFound, nodesVisited);

    if (!notFound) {
        // compute normal
        if (recordNormal && !primitiveTypeIsAggregate) {
            i.computeSilhouetteNormal(silhouetteRefs[i.referenceIndex]);
        }

        return true;
    }

    return false;
}

} // namespace fcpw