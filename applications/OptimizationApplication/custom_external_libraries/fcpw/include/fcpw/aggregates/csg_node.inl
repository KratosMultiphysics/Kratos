namespace fcpw {

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline void CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::computeBoundingBox()
{
    if (operation == BooleanOperation::Intersection) {
        // use the child bounding box with the smaller extent; this is not the tightest fit box
        BoundingBox<DIM> boxLeft = left->boundingBox();
        BoundingBox<DIM> boxRight = right->boundingBox();
        box.expandToInclude(boxLeft.extent().squaredNorm() <
                            boxRight.extent().squaredNorm() ?
                            boxLeft : boxRight);

    } else if (operation == BooleanOperation::Difference) {
        // use the bounding box of the left child (i.e., the object that is subtracted from);
        // this is not the tightest fit box
        box.expandToInclude(left->boundingBox());

    } else {
        // this is the tightest fit box for the union and none operations
        BoundingBox<DIM> boxLeft = left->boundingBox();
        BoundingBox<DIM> boxRight = right->boundingBox();
        box.expandToInclude(boxLeft);
        box.expandToInclude(boxRight);
    }
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::CsgNode(std::unique_ptr<PrimitiveTypeLeft> left_,
                                                                    std::unique_ptr<PrimitiveTypeRight> right_,
                                                                    const BooleanOperation& operation_):
left(std::move(left_)),
right(std::move(right_)),
operation(operation_),
leftPrimitiveTypeIsAggregate(std::is_base_of<Aggregate<DIM>, PrimitiveTypeLeft>::value),
rightPrimitiveTypeIsAggregate(std::is_base_of<Aggregate<DIM>, PrimitiveTypeRight>::value)
{
    // compute bounding box
    computeBoundingBox();
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline void CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::refit()
{
    if (leftPrimitiveTypeIsAggregate) left->refit();
    if (rightPrimitiveTypeIsAggregate) right->refit();
    computeBoundingBox();
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline BoundingBox<DIM> CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::boundingBox() const
{
    return box;
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline Vector<DIM> CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::centroid() const
{
    return box.centroid();
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline float CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::surfaceArea() const
{
    // NOTE: this is an overestimate
    return left->surfaceArea() + right->surfaceArea();
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline float CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::signedVolume() const
{
    // NOTE: these are overestimates
    float boxVolume = box.volume();
    if (boxVolume == 0.0f) boxVolume = maxFloat;

    if (operation == BooleanOperation::Intersection) {
        return std::min(boxVolume, std::min(left->signedVolume(), right->signedVolume()));

    } else if (operation == BooleanOperation::Difference) {
        return std::min(boxVolume, left->signedVolume());
    }

    return std::min(boxVolume, left->signedVolume() + right->signedVolume());
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline void CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::computeInteractions(const std::vector<Interaction<DIM>>& isLeft,
                                                                                     const std::vector<Interaction<DIM>>& isRight,
                                                                                     std::vector<Interaction<DIM>>& is) const
{
    int nLeft = 0;
    int nRight = 0;
    int hitsLeft = (int)isLeft.size();
    int hitsRight = (int)isRight.size();
    bool isLeftIntervalStart = hitsLeft%2 == 0;
    bool isRightIntervalStart = hitsRight%2 == (operation == BooleanOperation::Difference ? 1 : 0);
    int counter = 0;
    if (!isLeftIntervalStart) counter++;
    if (!isRightIntervalStart) counter++;

    auto addInteraction = [](const BooleanOperation& operation, int before, int after) -> bool {
        if (operation == BooleanOperation::Intersection || operation == BooleanOperation::Difference) {
            return (before == 1 && after == 2) || (before == 2 && after == 1);
        }

        // operation is union
        return (before == 0 && after == 1) || (before == 1 && after == 0);
    };

    // traverse the left & right interaction lists, appending interactions based on the operation
    while (nLeft != hitsLeft || nRight != hitsRight) {
        if (operation == BooleanOperation::Intersection && (nLeft == hitsLeft || nRight == hitsRight)) break;
        if (operation == BooleanOperation::Difference && nLeft == hitsLeft) break;

        int counterBefore = counter;
        if (nRight == hitsRight || (nLeft != hitsLeft && isLeft[nLeft].d < isRight[nRight].d)) {
            // left interaction is closer than right interaction
            counter += isLeftIntervalStart ? 1 : -1;
            isLeftIntervalStart = !isLeftIntervalStart;
            if (addInteraction(operation, counterBefore, counter)) is.emplace_back(isLeft[nLeft]);
            nLeft++;

        } else {
            // right interaction is closer than left interaction
            counter += isRightIntervalStart ? 1 : -1;
            isRightIntervalStart = !isRightIntervalStart;
            if (addInteraction(operation, counterBefore, counter)) {
                is.emplace_back(isRight[nRight]);
                if (operation == BooleanOperation::Difference) is[is.size() - 1].n *= -1; // flip normal if operation is difference
            }
            nRight++;
        }
    }
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline bool CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                                                                                   int aggregateIndex, int& nodesVisited,
                                                                                   bool checkForOcclusion) const
{
    std::cerr << "CsgNode::intersectFromNode(const Ray<DIM>&) not implemented" << std::endl;
    exit(EXIT_FAILURE);

    return false;
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline int CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                                                                  int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                                                                                  bool checkForOcclusion, bool recordAllHits) const
{
    // TODO: optimize for checkForOcclusion == true
    int hits = 0;
    is.clear();
    float tMin, tMax;
    nodesVisited++;

    if (box.intersect(r, tMin, tMax)) {
        // perform intersection query for left child
        int hitsLeft = 0;
        Ray<DIM> rLeft = r;
        std::vector<Interaction<DIM>> isLeft;
        if (leftPrimitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(left.get());
            hitsLeft = aggregate->intersectFromNode(rLeft, isLeft, nodeStartIndex, aggregateIndex,
                                                    nodesVisited, false, true);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(left.get());
            hitsLeft = geometricPrim->intersect(rLeft, isLeft, false, true);
        }

        // return if no intersections for the left child were found and
        // the operation is intersection or difference
        if (hitsLeft == 0 && (operation == BooleanOperation::Intersection ||
                              operation == BooleanOperation::Difference)) return 0;

        // perform intersection query for right child
        int hitsRight = 0;
        Ray<DIM> rRight = r;
        std::vector<Interaction<DIM>> isRight;
        if (rightPrimitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(right.get());
            hitsRight = aggregate->intersectFromNode(rRight, isRight, nodeStartIndex, aggregateIndex,
                                                     nodesVisited, false, true);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(right.get());
            hitsRight = geometricPrim->intersect(rRight, isRight, false, true);
        }

        // return if no intersections were found for both children
        if (hitsLeft == 0 && hitsRight == 0) return 0;

        if (hitsLeft > 0 && hitsRight > 0) {
            // determine interactions based on the operation
            if (operation == BooleanOperation::None) {
                // merge the left and right sorted interaction lists
                is.resize(isLeft.size() + isRight.size());
                std::merge(isLeft.begin(), isLeft.end(),
                           isRight.begin(), isRight.end(),
                           is.begin(), compareInteractions<DIM>);

            } else {
                computeInteractions(isLeft, isRight, is);
            }

        } else if (hitsLeft > 0) {
            // return if no intersections for the right child were found and the operation
            // is intersection
            if (operation == BooleanOperation::Intersection) return 0;

            // set the interactions to the left child's interactions for the
            // difference, union and none operations
            is = isLeft;

        } else if (hitsRight > 0) {
            // set the interactions to the right child's interactions for the
            // union and none operations
            is = isRight;
        }

        // shrink ray's tMax if possible
        if (!recordAllHits) {
            r.tMax = is[0].d; // list is already sorted
            is.resize(1);
            hits = 1;

        } else {
            hits = (int)is.size();
        }
    }

    return hits;
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline int CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::intersectFromNode(const BoundingSphere<DIM>& s,
                                                                                  std::vector<Interaction<DIM>>& is,
                                                                                  int nodeStartIndex, int aggregateIndex,
                                                                                  int& nodesVisited, bool recordOneHit) const
{
    std::cerr << "CsgNode::intersectFromNode(const BoundingSphere<DIM>&) not implemented" << std::endl;
    exit(EXIT_FAILURE);

    return 0;
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline int CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                  const Vector<DIM>& randNums, int nodeStartIndex,
                                                                                  int aggregateIndex, int& nodesVisited,
                                                                                  const std::function<float(float)>& branchTraversalWeight) const
{
    std::cerr << "CsgNode::intersectFromNode(const BoundingSphere<DIM>&) not implemented" << std::endl;
    exit(EXIT_FAILURE);

    return 0;
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline bool CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                          int nodeStartIndex, int aggregateIndex,
                                                                                          int& nodesVisited, bool recordNormal) const
{
    bool notFound = true;
    float d2Min, d2Max;
    nodesVisited++;

    if (box.overlap(s, d2Min, d2Max)) {
        // perform closest point query on left child
        bool foundLeft = false;
        Interaction<DIM> iLeft;
        BoundingSphere<DIM> sLeft = s;
        if (leftPrimitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(left.get());
            foundLeft = aggregate->findClosestPointFromNode(sLeft, iLeft, nodeStartIndex, aggregateIndex,
                                                            nodesVisited, true);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(left.get());
            foundLeft = geometricPrim->findClosestPoint(sLeft, iLeft);

            // compute normal
            if (recordNormal && foundLeft) {
                iLeft.computeNormal(left.get());
            }
        }

        // return if no closest point for the left child is found and
        // the operation is intersection or difference
        if (!foundLeft && (operation == BooleanOperation::Intersection ||
                           operation == BooleanOperation::Difference)) return false;

        // perform closest point query on right child
        bool foundRight = false;
        Interaction<DIM> iRight;
        BoundingSphere<DIM> sRight = s;
        if (rightPrimitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(right.get());
            foundRight = aggregate->findClosestPointFromNode(sRight, iRight, nodeStartIndex, aggregateIndex,
                                                             nodesVisited, true);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(right.get());
            foundRight = geometricPrim->findClosestPoint(sRight, iRight);

            // compute normal
            if (recordNormal && foundRight) {
                iRight.computeNormal(right.get());
            }
        }

        // return if no closest point was found to both children
        if (!foundLeft && !foundRight) return false;

        if (foundLeft && foundRight) {
            // compute signed distances
            float sdLeft = iLeft.signedDistance(s.c);
            float sdRight = iRight.signedDistance(s.c);
            DistanceInfo info = iLeft.distanceInfo == DistanceInfo::Exact &&
                                iRight.distanceInfo == DistanceInfo::Exact ?
                                DistanceInfo::Exact : DistanceInfo::Bounded;

            // determine which interaction to set and whether the distance info is
            // exact or bounded based on the operation
            if (operation == BooleanOperation::Union) {
                i = sdLeft < sdRight ? iLeft : iRight; // min(sdLeft, sdRight)
                i.distanceInfo = info == DistanceInfo::Exact &&
                                 sdLeft > 0 && sdRight > 0 ?
                                 DistanceInfo::Exact : DistanceInfo::Bounded;

            } else if (operation == BooleanOperation::Intersection) {
                i = sdLeft > sdRight ? iLeft : iRight; // max(sdLeft, sdRight)
                i.distanceInfo = info == DistanceInfo::Exact &&
                                 sdLeft < 0 && sdRight < 0 ?
                                 DistanceInfo::Exact : DistanceInfo::Bounded;

            } else if (operation == BooleanOperation::Difference) {
                iRight.n *= -1; // flip normal of right child
                iRight.sign *= -1; // flip sign of right child
                i = sdLeft > -sdRight ? iLeft : iRight; // max(sdLeft, -sdRight)
                i.distanceInfo = info == DistanceInfo::Exact &&
                                 sdLeft < 0 && sdRight > 0 ?
                                 DistanceInfo::Exact : DistanceInfo::Bounded;

            } else {
                // set the closer of the two interactions
                i = iLeft.d < iRight.d ? iLeft : iRight;
            }

        } else if (foundLeft) {
            // return if no closest point was found to the right child and the operation
            // is intersection
            if (operation == BooleanOperation::Intersection) return false;

            // set the interaction to the left child's interaction for the
            // difference, union and none operations
            i = iLeft;

        } else if (foundRight) {
            // set the interaction to the right child's interaction for the
            // union and none operations
            i = iRight;
        }

        // shrink sphere radius if possible
        s.r2 = std::min(s.r2, i.d*i.d);
        notFound = false;
    }

    return !notFound;
}

template<size_t DIM, typename PrimitiveTypeLeft, typename PrimitiveTypeRight>
inline bool CsgNode<DIM, PrimitiveTypeLeft, PrimitiveTypeRight>::findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                                    int nodeStartIndex, int aggregateIndex,
                                                                                                    int& nodesVisited, bool flipNormalOrientation,
                                                                                                    float squaredMinRadius, float precision,
                                                                                                    bool recordNormal) const
{
    std::cerr << "CsgNode::findClosestSilhouettePointFromNode() not implemented" << std::endl;
    exit(EXIT_FAILURE);

    return false;
}

} // namespace fcpw