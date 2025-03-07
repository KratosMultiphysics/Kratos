namespace fcpw {

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline float Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::computeSplitCost(const BoundingBox<DIM>& boxLeft,
                                                                                 const BoundingBox<DIM>& boxRight,
                                                                                 int nReferencesLeft, int nReferencesRight,
                                                                                 int depth) const
{
    float cost = maxFloat;
    if (packLeaves && depth > 0 && ((float)depthGuess/depth) < 1.5f &&
        nReferencesLeft%leafSize != 0 && nReferencesRight%leafSize != 0) {
        return cost;
    }

    if (costHeuristic == CostHeuristic::SurfaceArea) {
        cost = nReferencesLeft*boxLeft.surfaceArea() + nReferencesRight*boxRight.surfaceArea();

    } else if (costHeuristic == CostHeuristic::OverlapSurfaceArea) {
        // set the cost to be negative if the left and right boxes don't overlap at all
        BoundingBox<DIM> boxIntersected = boxLeft.intersect(boxRight);
        cost = (nReferencesLeft/boxRight.surfaceArea() +
                nReferencesRight/boxLeft.surfaceArea())*std::fabs(boxIntersected.surfaceArea());
        if (!boxIntersected.isValid()) cost *= -1;

    } else if (costHeuristic == CostHeuristic::Volume) {
        cost = nReferencesLeft*boxLeft.volume() + nReferencesRight*boxRight.volume();

    } else if (costHeuristic == CostHeuristic::OverlapVolume) {
        // set the cost to be negative if the left and right boxes don't overlap at all
        BoundingBox<DIM> boxIntersected = boxLeft.intersect(boxRight);
        cost = (nReferencesLeft/boxRight.volume() +
                nReferencesRight/boxLeft.volume())*std::fabs(boxIntersected.volume());
        if (!boxIntersected.isValid()) cost *= -1;
    }

    return cost;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline float Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::computeObjectSplit(const BoundingBox<DIM>& nodeBoundingBox,
                                                                                   const BoundingBox<DIM>& nodeCentroidBox,
                                                                                   const std::vector<BoundingBox<DIM>>& referenceBoxes,
                                                                                   const std::vector<Vector<DIM>>& referenceCentroids,
                                                                                   int depth, int nodeStart, int nodeEnd,
                                                                                   int& splitDim, float& splitCoord)
{
    float splitCost = maxFloat;
    splitDim = -1;
    splitCoord = 0.0f;

    if (costHeuristic != CostHeuristic::LongestAxisCenter) {
        Vector<DIM> extent = nodeBoundingBox.extent();

        // find the best split across all dimensions
        for (size_t dim = 0; dim < DIM; dim++) {
            // ignore flat dimension
            if (extent[dim] < 1e-6f) continue;

            // bin references into buckets
            float bucketWidth = extent[dim]/nBuckets;
            for (int b = 0; b < nBuckets; b++) {
                buckets[b].first = BoundingBox<DIM>();
                buckets[b].second = 0;
            }

            for (int p = nodeStart; p < nodeEnd; p++) {
                int bucketIndex = (int)((referenceCentroids[p][dim] - nodeBoundingBox.pMin[dim])/bucketWidth);
                bucketIndex = std::clamp(bucketIndex, 0, nBuckets - 1);
                buckets[bucketIndex].first.expandToInclude(referenceBoxes[p]);
                buckets[bucketIndex].second += 1;
            }

            // sweep right to left to build right bucket bounding boxes
            BoundingBox<DIM> boxRefRight;
            for (int b = nBuckets - 1; b > 0; b--) {
                boxRefRight.expandToInclude(buckets[b].first);
                rightBuckets[b].first = boxRefRight;
                rightBuckets[b].second = buckets[b].second;
                if (b != nBuckets - 1) rightBuckets[b].second += rightBuckets[b + 1].second;
            }

            // evaluate bucket split costs
            BoundingBox<DIM> boxRefLeft;
            int nReferencesLeft = 0;
            for (int b = 1; b < nBuckets; b++) {
                boxRefLeft.expandToInclude(buckets[b - 1].first);
                nReferencesLeft += buckets[b - 1].second;

                if (nReferencesLeft > 0 && rightBuckets[b].second > 0) {
                    float cost = computeSplitCost(boxRefLeft, rightBuckets[b].first,
                                                  nReferencesLeft, rightBuckets[b].second,
                                                  depth);

                    if (cost < splitCost) {
                        splitCost = cost;
                        splitDim = dim;
                        splitCoord = nodeBoundingBox.pMin[dim] + b*bucketWidth;
                    }
                }
            }
        }
    }

    // if no split dimension was chosen, fallback to LongestAxisCenter heuristic
    if (splitDim == -1) {
        splitDim = nodeCentroidBox.maxDimension();
        splitCoord = (nodeCentroidBox.pMin[splitDim] + nodeCentroidBox.pMax[splitDim])*0.5f;
    }

    return splitCost;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline int Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::performObjectSplit(int nodeStart, int nodeEnd, int splitDim, float splitCoord,
                                                                                 std::vector<BoundingBox<DIM>>& referenceBoxes,
                                                                                 std::vector<Vector<DIM>>& referenceCentroids)
{
    int mid = nodeStart;
    for (int i = nodeStart; i < nodeEnd; i++) {
        if (referenceCentroids[i][splitDim] < splitCoord) {
            std::swap(primitives[i], primitives[mid]);
            std::swap(referenceBoxes[i], referenceBoxes[mid]);
            std::swap(referenceCentroids[i], referenceCentroids[mid]);
            mid++;
        }
    }

    // if we get a bad split, just choose the center...
    if (mid == nodeStart || mid == nodeEnd) {
        mid = nodeStart + (nodeEnd - nodeStart)/2;

        // ensure the number of primitives in one branch is a multiple of the leaf size
        if (packLeaves) {
            while ((mid - nodeStart)%leafSize != 0 && mid < nodeEnd) mid++;
            if (mid == nodeEnd) mid = nodeStart + (nodeEnd - nodeStart)/2;
        }
    }

    return mid;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::buildRecursive(std::vector<BoundingBox<DIM>>& referenceBoxes,
                                                                              std::vector<Vector<DIM>>& referenceCentroids,
                                                                              std::vector<NodeType>& buildNodes,
                                                                              int parent, int start, int end, int depth)
{
    const int Untouched    = 0xffffffff;
    const int TouchedTwice = 0xfffffffd;
    maxDepth = std::max(depth, maxDepth);

    // add node to tree
    NodeType node;
    int currentNodeIndex = nNodes;
    int nReferences = end - start;

    nNodes++;

    // calculate the bounding box for this node
    BoundingBox<DIM> bb, bc;
    for (int p = start; p < end; p++) {
        bb.expandToInclude(referenceBoxes[p]);
        bc.expandToInclude(referenceCentroids[p]);
    }

    node.box = bb;

    // if the number of references at this point is less than the leaf
    // size, then this will become a leaf
    if (nReferences <= leafSize || depth == FCPW_BVH_MAX_DEPTH - 2) {
        node.referenceOffset = start;
        node.nReferences = nReferences;
        nLeafs++;

    } else {
        node.secondChildOffset = Untouched;
        node.nReferences = 0;
    }

    buildNodes.emplace_back(node);

    // child touches parent...
    // special case: don't do this for the root
    if (parent != 0xfffffffc) {
        buildNodes[parent].secondChildOffset--;

        // when this is the second touch, this is the right child;
        // the right child sets up the offset for the flat tree
        if (buildNodes[parent].secondChildOffset == TouchedTwice) {
            buildNodes[parent].secondChildOffset = nNodes - 1 - parent;
        }
    }

    // if this is a leaf, no need to subdivide
    if (node.nReferences > 0) return;

    // compute object split
    int splitDim;
    float splitCoord;
    float splitCost = computeObjectSplit(bb, bc, referenceBoxes, referenceCentroids,
                                         depth, start, end, splitDim, splitCoord);

    // partition the list of references on split
    int mid = performObjectSplit(start, end, splitDim, splitCoord, referenceBoxes, referenceCentroids);

    // push left and right children
    buildRecursive(referenceBoxes, referenceCentroids, buildNodes, currentNodeIndex, start, mid, depth + 1);
    buildRecursive(referenceBoxes, referenceCentroids, buildNodes, currentNodeIndex, mid, end, depth + 1);
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::build()
{
    // precompute bounding boxes and centroids
    int nReferences = (int)primitives.size();
    std::vector<BoundingBox<DIM>> referenceBoxes;
    std::vector<Vector<DIM>> referenceCentroids;

    referenceBoxes.resize(nReferences);
    referenceCentroids.resize(nReferences);
    flatTree.reserve(nReferences*2);

    for (int i = 0; i < nReferences; i++) {
        referenceBoxes[i] = primitives[i]->boundingBox();
        referenceCentroids[i] = primitives[i]->centroid();
    }

    // build tree recursively
    buildRecursive(referenceBoxes, referenceCentroids, flatTree, 0xfffffffc, 0, nReferences, 0);

    // clear working set
    buckets.clear();
    rightBuckets.clear();
}

template<size_t DIM, typename SilhouetteType>
inline void computeBoundingConesRecursive(const std::vector<SilhouetteType *>& silhouetteRefs,
                                          const std::vector<Vector<DIM>>& silhouetteNormals,
                                          const std::vector<Vector<DIM>>& silhouetteFaceNormals,
                                          std::vector<SnchNode<DIM>>& flatTree, int start, int end)
{
    BoundingCone<DIM> cone;
    SnchNode<DIM>& node(flatTree[start]);
    Vector<DIM> centroid = node.box.centroid();

    // compute bounding cone axis
    bool anySilhouetteRefs = false;
    bool silhouettesHaveTwoAdjacentFaces = true;
    for (int i = start; i < end; i++) {
        SnchNode<DIM>& childNode(flatTree[i]);

        for (int j = 0; j < childNode.nSilhouetteReferences; j++) { // is leaf if nSilhouetteReferences > 0
            int referenceIndex = childNode.silhouetteReferenceOffset + j;
            const SilhouetteType *silhouette = silhouetteRefs[referenceIndex];

            cone.axis += silhouetteNormals[referenceIndex];
            cone.radius = std::max(cone.radius, (silhouette->centroid() - centroid).norm());
            silhouettesHaveTwoAdjacentFaces = silhouettesHaveTwoAdjacentFaces &&
                                              silhouette->hasFace(0) &&
                                              silhouette->hasFace(1);
            anySilhouetteRefs = true;
        }
    }

    // compute bounding cone angle
    if (!anySilhouetteRefs) {
        node.cone.halfAngle = -M_PI;

    } else if (!silhouettesHaveTwoAdjacentFaces) {
        node.cone.halfAngle = M_PI;

    } else {
        float axisNorm = cone.axis.norm();
        if (axisNorm > epsilon) {
            cone.axis /= axisNorm;
            cone.halfAngle = 0.0f;

            for (int i = start; i < end; i++) {
                SnchNode<DIM>& childNode(flatTree[i]);

                for (int j = 0; j < childNode.nSilhouetteReferences; j++) { // is leaf if nSilhouetteReferences > 0
                    int referenceIndex = childNode.silhouetteReferenceOffset + j;
                    const SilhouetteType *silhouette = silhouetteRefs[referenceIndex];

                    for (int k = 0; k < 2; k++) {
                        const Vector<DIM>& n = silhouetteFaceNormals[2*referenceIndex + k];
                        float angle = std::acos(std::max(-1.0f, std::min(1.0f, cone.axis.dot(n))));
                        cone.halfAngle = std::max(cone.halfAngle, angle);
                    }
                }
            }

            node.cone = cone;
        }
    }

    // recurse on children
    if (node.nReferences == 0) { // not a leaf
        computeBoundingConesRecursive<DIM, SilhouetteType>(silhouetteRefs, silhouetteNormals, silhouetteFaceNormals,
                                                           flatTree, start + 1, start + node.secondChildOffset);
        computeBoundingConesRecursive<DIM, SilhouetteType>(silhouetteRefs, silhouetteNormals, silhouetteFaceNormals,
                                                           flatTree, start + node.secondChildOffset, end);
    }
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::assignGeometricDataToNodes(const std::function<bool(float, int)>& ignoreSilhouette)
{
    // do nothing
}

template<>
inline void Bvh<2, SnchNode<2>, LineSegment, SilhouetteVertex>::assignGeometricDataToNodes(const std::function<bool(float, int)>& ignoreSilhouette)
{
    // collect silhouette references
    Vector2 zero = Vector2::Zero();
    std::vector<Vector2> silhouetteRefNormals, silhouetteRefFaceNormals;

    for (int i = 0; i < (int)flatTree.size(); i++) {
        SnchNode<2>& node(flatTree[i]);
        std::unordered_map<int, bool> seenVertex;
        int start = (int)silhouetteRefs.size();

        for (int j = 0; j < node.nReferences; j++) { // leaf node if nReferences > 0
            int referenceIndex = node.referenceOffset + j;
            LineSegment *lineSegment = primitives[referenceIndex];

            for (int k = 0; k < 2; k++) {
                int vIndex = lineSegment->indices[k];
                SilhouetteVertex *silhouetteVertex = silhouettes[vIndex];

                if (seenVertex.find(vIndex) == seenVertex.end()) {
                    seenVertex[vIndex] = true;
                    bool ignore = false;
                    bool hasFace0 = silhouetteVertex->hasFace(0);
                    bool hasFace1 = silhouetteVertex->hasFace(1);
                    bool hasTwoAdjacentFaces = hasFace0 && hasFace1;
                    Vector2 n0 = hasFace0 ? silhouetteVertex->normal(0, true) : zero;
                    Vector2 n1 = hasFace1 ? silhouetteVertex->normal(1, true) : zero;
                    Vector2 n = silhouetteVertex->normal();

                    if (hasTwoAdjacentFaces && ignoreSilhouette) {
                        float det = n0[0]*n1[1] - n0[1]*n1[0];
                        ignore = ignoreSilhouette(det, lineSegment->getIndex());
                    }

                    if (!ignore) {
                        silhouetteRefs.emplace_back(silhouetteVertex);
                        silhouetteRefNormals.emplace_back(n);
                        silhouetteRefFaceNormals.emplace_back(n0);
                        silhouetteRefFaceNormals.emplace_back(n1);
                    }
                }
            }
        }

        int end = (int)silhouetteRefs.size();
        node.silhouetteReferenceOffset = start;
        node.nSilhouetteReferences = end - start;
    }

    // compute bounding cones recursively
    computeBoundingConesRecursive<2, SilhouetteVertex>(silhouetteRefs, silhouetteRefNormals,
                                                       silhouetteRefFaceNormals, flatTree,
                                                       0, (int)flatTree.size());
}

template<>
inline void Bvh<3, SnchNode<3>, Triangle, SilhouetteEdge>::assignGeometricDataToNodes(const std::function<bool(float, int)>& ignoreSilhouette)
{
    // collect silhouette references
    Vector3 zero = Vector3::Zero();
    std::vector<Vector3> silhouetteRefNormals, silhouetteRefFaceNormals;

    for (int i = 0; i < (int)flatTree.size(); i++) {
        SnchNode<3>& node(flatTree[i]);
        std::unordered_map<int, bool> seenEdge;
        int start = (int)silhouetteRefs.size();

        for (int j = 0; j < node.nReferences; j++) { // leaf node if nReferences > 0
            int referenceIndex = node.referenceOffset + j;
            Triangle *triangle = primitives[referenceIndex];

            for (int k = 0; k < 3; k++) {
                int eIndex = triangle->soup->eIndices[3*triangle->getIndex() + k];
                SilhouetteEdge *silhouetteEdge = silhouettes[eIndex];

                if (seenEdge.find(eIndex) == seenEdge.end()) {
                    seenEdge[eIndex] = true;
                    bool ignore = false;
                    bool hasFace0 = silhouetteEdge->hasFace(0);
                    bool hasFace1 = silhouetteEdge->hasFace(1);
                    bool hasTwoAdjacentFaces = hasFace0 && hasFace1;
                    Vector3 n0 = hasFace0 ? silhouetteEdge->normal(0, true) : zero;
                    Vector3 n1 = hasFace1 ? silhouetteEdge->normal(1, true) : zero;
                    Vector3 n = silhouetteEdge->normal();

                    if (hasTwoAdjacentFaces && ignoreSilhouette) {
                        const Vector3& pa = silhouetteEdge->soup->positions[silhouetteEdge->indices[1]];
                        const Vector3& pb = silhouetteEdge->soup->positions[silhouetteEdge->indices[2]];
                        Vector3 edgeDir = (pb - pa).normalized();
                        float dihedralAngle = std::atan2(edgeDir.dot(n0.cross(n1)), n0.dot(n1));
                        ignore = ignoreSilhouette(dihedralAngle, triangle->getIndex());
                    }

                    if (!ignore) {
                        silhouetteRefs.emplace_back(silhouetteEdge);
                        silhouetteRefNormals.emplace_back(n);
                        silhouetteRefFaceNormals.emplace_back(n0);
                        silhouetteRefFaceNormals.emplace_back(n1);
                    }
                }
            }
        }

        int end = (int)silhouetteRefs.size();
        node.silhouetteReferenceOffset = start;
        node.nSilhouetteReferences = end - start;
    }

    // compute bounding cones recursively
    computeBoundingConesRecursive<3, SilhouetteEdge>(silhouetteRefs, silhouetteRefNormals,
                                                     silhouetteRefFaceNormals, flatTree,
                                                     0, (int)flatTree.size());
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::Bvh(const CostHeuristic& costHeuristic_,
                                                              std::vector<PrimitiveType *>& primitives_,
                                                              std::vector<SilhouetteType *>& silhouettes_,
                                                              SortPositionsFunc<DIM, NodeType, PrimitiveType, SilhouetteType> sortPositions_,
                                                              const std::function<bool(float, int)>& ignoreSilhouette_,
                                                              bool packLeaves_, int leafSize_, int nBuckets_):
costHeuristic(costHeuristic_),
nNodes(0),
nLeafs(0),
leafSize(leafSize_),
nBuckets(nBuckets_),
maxDepth(0),
depthGuess(std::log2(primitives_.size())),
buckets(nBuckets, std::make_pair(BoundingBox<DIM>(), 0)),
rightBuckets(nBuckets, std::make_pair(BoundingBox<DIM>(), 0)),
primitives(primitives_),
silhouettes(silhouettes_),
packLeaves(packLeaves_),
primitiveTypeIsAggregate(std::is_base_of<Aggregate<DIM>, PrimitiveType>::value)
{
    // build bvh
    build();

    // sort positions
    if (sortPositions_) {
        sortPositions_(flatTree, primitives, silhouettes);
    }

    // assigns geometric data (e.g. cones) to nodes
    assignGeometricDataToNodes(ignoreSilhouette_);
}

template<typename NodeType>
inline void mergeBoundingCones(const NodeType& left, const NodeType& right, NodeType& node)
{
    // do nothing
}

template<size_t DIM>
inline void mergeBoundingCones(const SnchNode<DIM>& left, const SnchNode<DIM>& right, SnchNode<DIM>& node)
{
    node.cone = mergeBoundingCones<DIM>(left.cone, right.cone,
                                        left.box.centroid(),
                                        right.box.centroid(),
                                        node.box.centroid());
}

template<size_t DIM, typename SilhouetteType>
inline BoundingCone<DIM> computeBoundingCone(const std::vector<SilhouetteType *>& silhouetteRefs,
                                             const Vector<DIM>& centroid,
                                             int nSilhouetteReferences,
                                             int silhouetteReferenceOffset)
{
    // compute bounding cone axis
    BoundingCone<DIM> cone;
    bool anySilhouetteRefs = false;
    bool silhouettesHaveTwoAdjacentFaces = true;
    for (int p = 0; p < nSilhouetteReferences; p++) {
        int referenceIndex = silhouetteReferenceOffset + p;
        const SilhouetteType *silhouette = silhouetteRefs[referenceIndex];

        cone.axis += silhouette->normal();
        cone.radius = std::max(cone.radius, (silhouette->centroid() - centroid).norm());
        silhouettesHaveTwoAdjacentFaces = silhouettesHaveTwoAdjacentFaces &&
                                          silhouette->hasFace(0) &&
                                          silhouette->hasFace(1);
        anySilhouetteRefs = true;
    }

    // compute bounding cone angle
    if (!anySilhouetteRefs) {
        cone.halfAngle = -M_PI;

    } else if (!silhouettesHaveTwoAdjacentFaces) {
        cone.halfAngle = M_PI;

    } else {
        float axisNorm = cone.axis.norm();
        if (axisNorm > epsilon) {
            cone.axis /= axisNorm;
            cone.halfAngle = 0.0f;

            for (int p = 0; p < nSilhouetteReferences; p++) {
                int referenceIndex = silhouetteReferenceOffset + p;
                const SilhouetteType *silhouette = silhouetteRefs[referenceIndex];

                for (int k = 0; k < 2; k++) {
                    Vector<DIM> n = silhouette->normal(k, true);
                    float angle = std::acos(std::max(-1.0f, std::min(1.0f, cone.axis.dot(n))));
                    cone.halfAngle = std::max(cone.halfAngle, angle);
                }
            }
        }
    }

    return cone;
}

template<typename NodeType, typename SilhouetteType>
inline void computeBoundingCone(const std::vector<SilhouetteType *>& silhouetteRefs, NodeType& node)
{
    // do nothing
}

template<size_t DIM, typename SilhouetteType>
inline void computeBoundingCone(const std::vector<SilhouetteType *>& silhouetteRefs, SnchNode<DIM>& node)
{
    node.cone = computeBoundingCone<DIM, SilhouetteType>(silhouetteRefs, node.box.centroid(),
                                                         node.nSilhouetteReferences,
                                                         node.silhouetteReferenceOffset);
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void refitRecursive(const std::vector<PrimitiveType *>& primitives,
                           const std::vector<SilhouetteType *>& silhouetteRefs,
                           std::vector<NodeType>& flatTree, int nodeIndex)
{
    NodeType& node(flatTree[nodeIndex]);

    if (node.nReferences == 0) { // not a leaf
        refitRecursive<DIM, NodeType, PrimitiveType, SilhouetteType>(
            primitives, silhouetteRefs, flatTree, nodeIndex + 1);
        refitRecursive<DIM, NodeType, PrimitiveType, SilhouetteType>(
            primitives, silhouetteRefs, flatTree, nodeIndex + node.secondChildOffset);

        // merge left and right child bounding boxes
        node.box = flatTree[nodeIndex + 1].box;
        node.box.expandToInclude(flatTree[nodeIndex + node.secondChildOffset].box);

        // merge left and right child bounding cones
        mergeBoundingCones(flatTree[nodeIndex + 1], flatTree[nodeIndex + node.secondChildOffset], node);

    } else { // leaf
        // compute bounding box
        node.box = BoundingBox<DIM>();
        for (int p = 0; p < node.nReferences; p++) {
            int referenceIndex = node.referenceOffset + p;
            const PrimitiveType *prim = primitives[referenceIndex];

            node.box.expandToInclude(prim->boundingBox());
        }

        // compute bounding cone
        computeBoundingCone(silhouetteRefs, node);
    }
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::refit()
{
    // refit primitives if they are aggregates
    if (primitiveTypeIsAggregate) {
        for (int p = 0; p < (int)primitives.size(); p++) {
            Aggregate<DIM> *aggregate = reinterpret_cast<Aggregate<DIM> *>(primitives[p]);
            aggregate->refit();
        }
    }

    // refit
    if (nNodes > 0) {
        refitRecursive<DIM, NodeType, PrimitiveType, SilhouetteType>(
            primitives, silhouetteRefs, flatTree, 0);
    }
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::printStats() const
{
    std::cout << "BVH stats: "
              << nNodes << " nodes, "
              << nLeafs << " leaves, "
              << maxDepth << " max depth, "
              << primitives.size() << " primitives, "
              << silhouettes.size() << " silhouettes and "
              << silhouetteRefs.size() << " silhouetteRefs"
              << std::endl;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline BoundingBox<DIM> Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::boundingBox() const
{
    return flatTree.size() > 0 ? flatTree[0].box : BoundingBox<DIM>();
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline Vector<DIM> Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::centroid() const
{
    Vector<DIM> c = Vector<DIM>::Zero();
    int nPrimitives = (int)primitives.size();

    for (int p = 0; p < nPrimitives; p++) {
        c += primitives[p]->centroid();
    }

    return c/nPrimitives;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline float Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::surfaceArea() const
{
    float area = 0.0f;
    for (int p = 0; p < (int)primitives.size(); p++) {
        area += primitives[p]->surfaceArea();
    }

    return area;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline float Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::signedVolume() const
{
    float volume = 0.0f;
    for (int p = 0; p < (int)primitives.size(); p++) {
        volume += primitives[p]->signedVolume();
    }

    return volume;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline bool Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::processSubtreeForIntersection(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                                                                                             int aggregateIndex, bool checkForOcclusion,
                                                                                             TraversalStack *subtree, float *boxHits,
                                                                                             bool& didHit, int& nodesVisited) const
{
    int stackPtr = 0;
    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found intersection, continue
        if (currentDist > r.tMax) continue;
        const NodeType& node(flatTree[nodeIndex]);

        // is leaf -> intersect
        if (node.nReferences > 0) {
            for (int p = 0; p < node.nReferences; p++) {
                int referenceIndex = node.referenceOffset + p;
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

        } else { // not a leaf
            bool hit0 = flatTree[nodeIndex + 1].box.intersect(r, boxHits[0], boxHits[1]);
            bool hit1 = flatTree[nodeIndex + node.secondChildOffset].box.intersect(r, boxHits[2], boxHits[3]);

            // did we hit both nodes?
            if (hit0 && hit1) {
                // we assume that the left child is a closer hit...
                int closer = nodeIndex + 1;
                int other = nodeIndex + node.secondChildOffset;

                // ... if the right child was actually closer, swap the relavent values
                if (boxHits[2] < boxHits[0]) {
                    std::swap(boxHits[0], boxHits[2]);
                    std::swap(boxHits[1], boxHits[3]);
                    std::swap(closer, other);
                }

                // it's possible that the nearest object is still in the other side, but we'll
                // check the farther-away node later...

                // push the farther first, then the closer
                stackPtr++;
                subtree[stackPtr].node = other;
                subtree[stackPtr].distance = boxHits[2];

                stackPtr++;
                subtree[stackPtr].node = closer;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit0) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + 1;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit1) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + node.secondChildOffset;
                subtree[stackPtr].distance = boxHits[2];
            }

            nodesVisited++;
        }
    }

    return false;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline bool Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                                                                                 int aggregateIndex, int& nodesVisited,
                                                                                 bool checkForOcclusion) const
{
    bool didHit = false;
    TraversalStack subtree[FCPW_BVH_MAX_DEPTH];
    float boxHits[4];

    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    if (flatTree[rootIndex].box.intersect(r, boxHits[0], boxHits[1])) {
        subtree[0].node = rootIndex;
        subtree[0].distance = boxHits[0];
        bool occluded = processSubtreeForIntersection(r, i, nodeStartIndex, aggregateIndex, checkForOcclusion,
                                                      subtree, boxHits, didHit, nodesVisited);
        if (occluded) return true;
    }

    return didHit;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline bool Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::processSubtreeForIntersection(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                                                                             int nodeStartIndex, int aggregateIndex, bool checkForOcclusion,
                                                                                             bool recordAllHits, TraversalStack *subtree,
                                                                                             float *boxHits, int& hits, int& nodesVisited) const
{
    int stackPtr = 0;
    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found intersection, continue
        if (!recordAllHits && currentDist > r.tMax) continue;
        const NodeType& node(flatTree[nodeIndex]);

        // is leaf -> intersect
        if (node.nReferences > 0) {
            for (int p = 0; p < node.nReferences; p++) {
                int referenceIndex = node.referenceOffset + p;
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
                        return true;
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

        } else { // not a leaf
            bool hit0 = flatTree[nodeIndex + 1].box.intersect(r, boxHits[0], boxHits[1]);
            bool hit1 = flatTree[nodeIndex + node.secondChildOffset].box.intersect(r, boxHits[2], boxHits[3]);

            // did we hit both nodes?
            if (hit0 && hit1) {
                // we assume that the left child is a closer hit...
                int closer = nodeIndex + 1;
                int other = nodeIndex + node.secondChildOffset;

                // ... if the right child was actually closer, swap the relavent values
                if (boxHits[2] < boxHits[0]) {
                    std::swap(boxHits[0], boxHits[2]);
                    std::swap(boxHits[1], boxHits[3]);
                    std::swap(closer, other);
                }

                // it's possible that the nearest object is still in the other side, but we'll
                // check the farther-away node later...

                // push the farther first, then the closer
                stackPtr++;
                subtree[stackPtr].node = other;
                subtree[stackPtr].distance = boxHits[2];

                stackPtr++;
                subtree[stackPtr].node = closer;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit0) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + 1;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit1) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + node.secondChildOffset;
                subtree[stackPtr].distance = boxHits[2];
            }

            nodesVisited++;
        }
    }

    return false;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline int Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                                                                int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                                                                                bool checkForOcclusion, bool recordAllHits) const
{
    int hits = 0;
    if (!recordAllHits) is.resize(1);
    TraversalStack subtree[FCPW_BVH_MAX_DEPTH];
    float boxHits[4];

    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    if (flatTree[rootIndex].box.intersect(r, boxHits[0], boxHits[1])) {
        subtree[0].node = rootIndex;
        subtree[0].distance = boxHits[0];
        bool occluded = processSubtreeForIntersection(r, is, nodeStartIndex, aggregateIndex, checkForOcclusion,
                                                      recordAllHits, subtree, boxHits, hits, nodesVisited);
        if (occluded) return 1;
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

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline float Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::processSubtreeForIntersection(const BoundingSphere<DIM>& s, std::vector<Interaction<DIM>>& is,
                                                                                              int nodeStartIndex, int aggregateIndex, bool recordOneHit,
                                                                                              TraversalStack *subtree, float *boxHits, int& hits, int& nodesVisited) const
{
    float totalPrimitiveWeight = 0.0f;
    int stackPtr = 0;

    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        const NodeType& node(flatTree[nodeIndex]);
        stackPtr--;

        // is leaf -> intersect
        if (node.nReferences > 0) {
            for (int p = 0; p < node.nReferences; p++) {
                int referenceIndex = node.referenceOffset + p;
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

        } else { // not a leaf
            bool hit0 = flatTree[nodeIndex + 1].box.overlap(s, boxHits[0], boxHits[1]);
            bool hit1 = flatTree[nodeIndex + node.secondChildOffset].box.overlap(s, boxHits[2], boxHits[3]);

            if (hit0) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + 1;
            }

            if (hit1) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + node.secondChildOffset;
            }

            nodesVisited++;
        }
    }

    return totalPrimitiveWeight;
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline int Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::intersectFromNode(const BoundingSphere<DIM>& s,
                                                                                std::vector<Interaction<DIM>>& is,
                                                                                int nodeStartIndex, int aggregateIndex,
                                                                                int& nodesVisited, bool recordOneHit) const
{
    int hits = 0;
    float totalPrimitiveWeight = 0.0f;
    if (recordOneHit && !primitiveTypeIsAggregate) is.resize(1);
    TraversalStack subtree[FCPW_BVH_MAX_DEPTH];
    float boxHits[4];

    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    if (flatTree[rootIndex].box.overlap(s, boxHits[0], boxHits[1])) {
        subtree[0].node = rootIndex;
        subtree[0].distance = s.r2;
        totalPrimitiveWeight = processSubtreeForIntersection(s, is, nodeStartIndex, aggregateIndex, recordOneHit,
                                                             subtree, boxHits, hits, nodesVisited);
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

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::processSubtreeForIntersection(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                             const Vector<DIM>& randNums, int nodeStartIndex, int aggregateIndex,
                                                                                             const std::function<float(float)>& branchTraversalWeight,
                                                                                             int nodeIndex, float traversalPdf, float *boxHits,
                                                                                             int& hits, int& nodesVisited) const
{
    int stackPtr = 0;
    float d2NodeMax = boxHits[1];
    float u = randNums[0];

    while (stackPtr >= 0) {
        // pop off the next node to work on
        const NodeType& node(flatTree[nodeIndex]);
        stackPtr--;

        // is leaf -> intersect
        if (node.nReferences > 0) {
            float totalPrimitiveWeight = 0.0f;
            for (int p = 0; p < node.nReferences; p++) {
                int referenceIndex = node.referenceOffset + p;
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

            if (totalPrimitiveWeight > 0.0f) {
                i.d /= totalPrimitiveWeight;
            }

        } else { // not a leaf
            const BoundingBox<DIM>& box0(flatTree[nodeIndex + 1].box);
            float weight0 = box0.overlap(s, boxHits[0], boxHits[1]) ? 1.0f : 0.0f;

            const BoundingBox<DIM>& box1(flatTree[nodeIndex + node.secondChildOffset].box);
            float weight1 = box1.overlap(s, boxHits[2], boxHits[3]) ? 1.0f : 0.0f;

            if (branchTraversalWeight) {
                if (weight0 > 0.0f) weight0 *= branchTraversalWeight((s.c - box0.centroid()).squaredNorm());
                if (weight1 > 0.0f) weight1 *= branchTraversalWeight((s.c - box1.centroid()).squaredNorm());
            }

            float totalTraversalWeight = weight0 + weight1;
            if (totalTraversalWeight > 0.0f) {
                stackPtr++;
                float traversalProb0 = weight0/totalTraversalWeight;
                float traversalProb1 = 1.0f - traversalProb0;

                if (u < traversalProb0) {
                    u = u/traversalProb0; // rescale to [0,1)
                    nodeIndex = nodeIndex + 1;
                    traversalPdf *= traversalProb0;
                    d2NodeMax = boxHits[1];

                } else {
                    u = (u - traversalProb0)/traversalProb1; // rescale to [0,1)
                    nodeIndex = nodeIndex + node.secondChildOffset;
                    traversalPdf *= traversalProb1;
                    d2NodeMax = boxHits[3];
                }
            }

            nodesVisited++;
        }
    }
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline int Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                const Vector<DIM>& randNums, int nodeStartIndex,
                                                                                int aggregateIndex, int& nodesVisited,
                                                                                const std::function<float(float)>& branchTraversalWeight) const
{
    int hits = 0;
    float boxHits[4];

    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    if (flatTree[rootIndex].box.overlap(s, boxHits[0], boxHits[1])) {
        processSubtreeForIntersection(s, i, randNums, nodeStartIndex, aggregateIndex, branchTraversalWeight,
                                      rootIndex, 1.0f, boxHits, hits, nodesVisited);
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

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::processSubtreeForClosestPoint(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                             int nodeStartIndex, int aggregateIndex,
                                                                                             bool recordNormal, TraversalStack *subtree,
                                                                                             float *boxHits, bool& notFound, int& nodesVisited) const
{
    int stackPtr = 0;
    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found primitive, continue
        if (currentDist > s.r2) continue;
        const NodeType& node(flatTree[nodeIndex]);

        // is leaf -> compute squared distance
        if (node.nReferences > 0) {
            for (int p = 0; p < node.nReferences; p++) {
                int referenceIndex = node.referenceOffset + p;
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

        } else { // not a leaf
            bool hit0 = flatTree[nodeIndex + 1].box.overlap(s, boxHits[0], boxHits[1]);
            s.r2 = std::min(s.r2, boxHits[1]);

            bool hit1 = flatTree[nodeIndex + node.secondChildOffset].box.overlap(s, boxHits[2], boxHits[3]);
            s.r2 = std::min(s.r2, boxHits[3]);

            // is there overlap with both nodes?
            if (hit0 && hit1) {
                // we assume that the left child is a closer hit...
                int closer = nodeIndex + 1;
                int other = nodeIndex + node.secondChildOffset;

                // ... if the right child was actually closer, swap the relavent values
                if (boxHits[0] == 0.0f && boxHits[2] == 0.0f) {
                    if (boxHits[3] < boxHits[1]) {
                        std::swap(closer, other);
                    }

                } else if (boxHits[2] < boxHits[0]) {
                    std::swap(boxHits[0], boxHits[2]);
                    std::swap(closer, other);
                }

                // it's possible that the nearest object is still in the other side, but we'll
                // check the farther-away node later...

                // push the farther first, then the closer
                stackPtr++;
                subtree[stackPtr].node = other;
                subtree[stackPtr].distance = boxHits[2];

                stackPtr++;
                subtree[stackPtr].node = closer;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit0) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + 1;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit1) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + node.secondChildOffset;
                subtree[stackPtr].distance = boxHits[2];
            }

            nodesVisited++;
        }
    }
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline bool Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                        int nodeStartIndex, int aggregateIndex,
                                                                                        int& nodesVisited, bool recordNormal) const
{
    bool notFound = true;
    TraversalStack subtree[FCPW_BVH_MAX_DEPTH];
    float boxHits[4];

    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    if (flatTree[rootIndex].box.overlap(s, boxHits[0], boxHits[1])) {
        s.r2 = std::min(s.r2, boxHits[1]);
        subtree[0].node = rootIndex;
        subtree[0].distance = boxHits[0];
        processSubtreeForClosestPoint(s, i, nodeStartIndex, aggregateIndex, recordNormal,
                                      subtree, boxHits, notFound, nodesVisited);
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

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void processSubtreeForClosestSilhouettePoint(const std::vector<NodeType>& flatTree,
                                                    const std::vector<PrimitiveType *>& primitives,
                                                    const std::vector<SilhouetteType *>& silhouetteRefs,
                                                    BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                    int nodeStartIndex, int aggregateIndex, int objectIndex,
                                                    bool primitiveTypeIsAggregate, bool flipNormalOrientation,
                                                    float squaredMinRadius, float precision, bool recordNormal,
                                                    TraversalStack *subtree, float *boxHits, bool& notFound,
                                                    int& nodesVisited)
{
    std::cerr << "Bvh::processSubtreeForClosestSilhouettePoint() not implemented" << std::endl;
    exit(EXIT_FAILURE);
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline void processSubtreeForClosestSilhouettePoint(const std::vector<SnchNode<DIM>>& flatTree,
                                                    const std::vector<PrimitiveType *>& primitives,
                                                    const std::vector<SilhouetteType *>& silhouetteRefs,
                                                    BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                    int nodeStartIndex, int aggregateIndex, int objectIndex,
                                                    bool primitiveTypeIsAggregate, bool flipNormalOrientation,
                                                    float squaredMinRadius, float precision, bool recordNormal,
                                                    TraversalStack *subtree, float *boxHits, bool& notFound,
                                                    int& nodesVisited)
{
    float stubs[2];
    int stackPtr = 0;
    while (stackPtr >= 0) {
        // pop off the next node to work on
        int nodeIndex = subtree[stackPtr].node;
        float currentDist = subtree[stackPtr].distance;
        stackPtr--;

        // if this node is further than the closest found silhouette, continue
        if (currentDist > s.r2) continue;
        const SnchNode<DIM>& node(flatTree[nodeIndex]);

        // is leaf -> compute silhouette distance
        if (node.nReferences > 0) {
            if (primitiveTypeIsAggregate) {
                for (int p = 0; p < node.nReferences; p++) {
                    int referenceIndex = node.referenceOffset + p;
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
                for (int p = 0; p < node.nSilhouetteReferences; p++) {
                    int referenceIndex = node.silhouetteReferenceOffset + p;
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

        } else { // not a leaf
            const SnchNode<DIM>& node0(flatTree[nodeIndex + 1]);
            bool hit0 = node0.cone.isValid() && node0.box.overlap(s, boxHits[0]) &&
                        node0.cone.overlap(s.c, node0.box, boxHits[0], stubs[0], stubs[1]);

            const SnchNode<DIM>& node1(flatTree[nodeIndex + node.secondChildOffset]);
            bool hit1 = node1.cone.isValid() && node1.box.overlap(s, boxHits[1]) &&
                        node1.cone.overlap(s.c, node1.box, boxHits[1], stubs[0], stubs[1]);

            // is there overlap with both nodes?
            if (hit0 && hit1) {
                // we assume that the left child is a closer hit...
                int closer = nodeIndex + 1;
                int other = nodeIndex + node.secondChildOffset;

                // ... if the right child was actually closer, swap the relavent values
                if (boxHits[1] < boxHits[0]) {
                    std::swap(boxHits[0], boxHits[1]);
                    std::swap(closer, other);
                }

                // it's possible that the nearest object is still in the other side, but we'll
                // check the farther-away node later...

                // push the farther first, then the closer
                stackPtr++;
                subtree[stackPtr].node = other;
                subtree[stackPtr].distance = boxHits[1];

                stackPtr++;
                subtree[stackPtr].node = closer;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit0) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + 1;
                subtree[stackPtr].distance = boxHits[0];

            } else if (hit1) {
                stackPtr++;
                subtree[stackPtr].node = nodeIndex + node.secondChildOffset;
                subtree[stackPtr].distance = boxHits[1];
            }

            nodesVisited++;
        }
    }
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline bool Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>::findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                                  int nodeStartIndex, int aggregateIndex,
                                                                                                  int& nodesVisited, bool flipNormalOrientation,
                                                                                                  float squaredMinRadius, float precision,
                                                                                                  bool recordNormal) const
{
    if (squaredMinRadius >= s.r2) return false;

    bool notFound = true;
    TraversalStack subtree[FCPW_BVH_MAX_DEPTH];
    float boxHits[2];

    int rootIndex = aggregateIndex == this->pIndex ? nodeStartIndex : 0;
    if (flatTree[rootIndex].box.overlap(s, boxHits[0])) {
        subtree[0].node = rootIndex;
        subtree[0].distance = boxHits[0];
        processSubtreeForClosestSilhouettePoint(flatTree, primitives, silhouetteRefs, s, i, 
                                                nodeStartIndex, aggregateIndex, this->pIndex,
                                                primitiveTypeIsAggregate, flipNormalOrientation,
                                                squaredMinRadius, precision, recordNormal,
                                                subtree, boxHits, notFound, nodesVisited);
    }

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