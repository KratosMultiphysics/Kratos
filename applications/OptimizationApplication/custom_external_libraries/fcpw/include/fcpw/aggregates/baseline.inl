namespace fcpw {

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline Baseline<DIM, PrimitiveType, SilhouetteType>::Baseline(std::vector<PrimitiveType *>& primitives_,
                                                              std::vector<SilhouetteType *>& silhouettes_):
primitives(primitives_),
silhouettes(silhouettes_),
primitiveTypeIsAggregate(std::is_base_of<Aggregate<DIM>, PrimitiveType>::value)
{

}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline void Baseline<DIM, PrimitiveType, SilhouetteType>::refit()
{
    // refit primitives if they are aggregates
    if (primitiveTypeIsAggregate) {
        for (int p = 0; p < (int)primitives.size(); p++) {
            Aggregate<DIM> *aggregate = reinterpret_cast<Aggregate<DIM> *>(primitives[p]);
            aggregate->refit();
        }
    }
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline BoundingBox<DIM> Baseline<DIM, PrimitiveType, SilhouetteType>::boundingBox() const
{
    BoundingBox<DIM> bb;
    for (int p = 0; p < (int)primitives.size(); p++) {
        bb.expandToInclude(primitives[p]->boundingBox());
    }

    return bb;
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline Vector<DIM> Baseline<DIM, PrimitiveType, SilhouetteType>::centroid() const
{
    Vector<DIM> c = Vector<DIM>::Zero();
    int nPrimitives = (int)primitives.size();

    for (int p = 0; p < nPrimitives; p++) {
        c += primitives[p]->centroid();
    }

    return c/nPrimitives;
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline float Baseline<DIM, PrimitiveType, SilhouetteType>::surfaceArea() const
{
    float area = 0.0f;
    for (int p = 0; p < (int)primitives.size(); p++) {
        area += primitives[p]->surfaceArea();
    }

    return area;
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline float Baseline<DIM, PrimitiveType, SilhouetteType>::signedVolume() const
{
    float volume = 0.0f;
    for (int p = 0; p < (int)primitives.size(); p++) {
        volume += primitives[p]->signedVolume();
    }

    return volume;
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline bool Baseline<DIM, PrimitiveType, SilhouetteType>::intersectFromNode(Ray<DIM>& r, Interaction<DIM>& i, int nodeStartIndex,
                                                                            int aggregateIndex, int& nodesVisited,
                                                                            bool checkForOcclusion) const
{
    // find closest hit
    bool didHit = false;
    for (int p = 0; p < (int)primitives.size(); p++) {
        nodesVisited++;

        bool hit = false;
        if (primitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(primitives[p]);
            hit = aggregate->intersectFromNode(r, i, nodeStartIndex, aggregateIndex,
                                               nodesVisited, checkForOcclusion);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(primitives[p]);
            hit = geometricPrim->intersect(r, i, checkForOcclusion);
        }

        if (hit) {
            if (checkForOcclusion) {
                return true;
            }

            didHit = true;
            r.tMax = std::min(r.tMax, i.d);
            i.referenceIndex = p;
            i.objectIndex = this->pIndex;
        }
    }

    return didHit;
}

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline int Baseline<DIM, PrimitiveType, SilhouetteType>::intersectFromNode(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                                                           int nodeStartIndex, int aggregateIndex, int& nodesVisited,
                                                                           bool checkForOcclusion, bool recordAllHits) const
{
    int hits = 0;
    if (!recordAllHits) is.resize(1);

    // find closest hit
    for (int p = 0; p < (int)primitives.size(); p++) {
        nodesVisited++;

        int hit = 0;
        std::vector<Interaction<DIM>> cs;
        if (primitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(primitives[p]);
            hit = aggregate->intersectFromNode(r, cs, nodeStartIndex, aggregateIndex,
                                               nodesVisited, checkForOcclusion, recordAllHits);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(primitives[p]);
            hit = geometricPrim->intersect(r, cs, checkForOcclusion, recordAllHits);
            for (int i = 0; i < (int)cs.size(); i++) {
                cs[i].referenceIndex = p;
                cs[i].objectIndex = this->pIndex;
            }
        }

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

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline int Baseline<DIM, PrimitiveType, SilhouetteType>::intersectFromNode(const BoundingSphere<DIM>& s,
                                                                           std::vector<Interaction<DIM>>& is,
                                                                           int nodeStartIndex, int aggregateIndex,
                                                                           int& nodesVisited, bool recordOneHit) const
{
    int hits = 0;
    float totalPrimitiveWeight = 0.0f;
    if (recordOneHit && !primitiveTypeIsAggregate) is.resize(1);

    for (int p = 0; p < (int)primitives.size(); p++) {
        nodesVisited++;

        if (primitiveTypeIsAggregate) {
            std::vector<Interaction<DIM>> cs;
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(primitives[p]);
            hits += aggregate->intersectFromNode(s, cs, nodeStartIndex, aggregateIndex,
                                                 nodesVisited, recordOneHit);
            is.insert(is.end(), cs.begin(), cs.end());

        } else {
            Interaction<DIM> c;
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(primitives[p]);
            bool hit = geometricPrim->intersect(s, c, recordOneHit);

            if (hit) {
                hits += 1;
                c.referenceIndex = p;
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

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline int Baseline<DIM, PrimitiveType, SilhouetteType>::intersectFromNode(const BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                           const Vector<DIM>& randNums, int nodeStartIndex,
                                                                           int aggregateIndex, int& nodesVisited,
                                                                           const std::function<float(float)>& branchTraversalWeight) const
{
    int hits = 0;
    float totalPrimitiveWeight = 0.0f;
    float u = randNums[0];

    for (int p = 0; p < (int)primitives.size(); p++) {
        nodesVisited++;

        int hit = 0;
        Interaction<DIM> c;
        if (primitiveTypeIsAggregate) {
            Vector<DIM> modifiedRandNums;
            modifiedRandNums[0] = u;
            for (int i = 1; i < DIM; i++) modifiedRandNums[i] = randNums[i];
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(primitives[p]);
            hit = aggregate->intersectFromNode(s, c, modifiedRandNums, nodeStartIndex, aggregateIndex,
                                               nodesVisited, branchTraversalWeight);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(primitives[p]);
            hit = geometricPrim->intersect(s, c, true) ? 1 : 0;
            c.referenceIndex = p;
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

            } else {
                u = (u - selectionProb)/(1.0f - selectionProb);
            }
        }
    }

    if (hits > 0) {
        if (totalPrimitiveWeight > 0.0f) {
            i.d /= totalPrimitiveWeight;
        }

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

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline bool Baseline<DIM, PrimitiveType, SilhouetteType>::findClosestPointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                   int nodeStartIndex, int aggregateIndex,
                                                                                   int& nodesVisited, bool recordNormal) const
{
    // find closest point
    bool notFound = true;
    for (int p = 0; p < (int)primitives.size(); p++) {
        nodesVisited++;

        bool found = false;
        Interaction<DIM> c;
        if (primitiveTypeIsAggregate) {
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(primitives[p]);
            found = aggregate->findClosestPointFromNode(s, c, nodeStartIndex, aggregateIndex,
                                                        nodesVisited, recordNormal);

        } else {
            const GeometricPrimitive<DIM> *geometricPrim = reinterpret_cast<const GeometricPrimitive<DIM> *>(primitives[p]);
            found = geometricPrim->findClosestPoint(s, c);
            c.referenceIndex = p;
            c.objectIndex = this->pIndex;
        }

        // keep the closest point only
        if (found) {
            notFound = false;
            s.r2 = std::min(s.r2, c.d*c.d);
            i = c;
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

template<size_t DIM, typename PrimitiveType, typename SilhouetteType>
inline bool Baseline<DIM, PrimitiveType, SilhouetteType>::findClosestSilhouettePointFromNode(BoundingSphere<DIM>& s, Interaction<DIM>& i,
                                                                                             int nodeStartIndex, int aggregateIndex,
                                                                                             int& nodesVisited, bool flipNormalOrientation,
                                                                                             float squaredMinRadius, float precision,
                                                                                             bool recordNormal) const
{
    if (squaredMinRadius >= s.r2) return false;

    bool notFound = true;
    if (primitiveTypeIsAggregate) {
        for (int p = 0; p < (int)primitives.size(); p++) {
            nodesVisited++;
            Interaction<DIM> c;
            const Aggregate<DIM> *aggregate = reinterpret_cast<const Aggregate<DIM> *>(primitives[p]);
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
        for (int p = 0; p < (int)silhouettes.size(); p++) {
            nodesVisited++;
            Interaction<DIM> c;
            bool found = silhouettes[p]->findClosestSilhouettePoint(s, c, flipNormalOrientation,
                                                                    squaredMinRadius, precision);

            // keep the closest silhouette point
            if (found) {
                notFound = false;
                s.r2 = std::min(s.r2, c.d*c.d);
                i = c;
                i.referenceIndex = p;
                i.objectIndex = this->pIndex;

                if (squaredMinRadius >= s.r2) {
                    break;
                }
            }
        }
    }

    if (!notFound) {
        // compute normal
        if (recordNormal && !primitiveTypeIsAggregate) {
            i.computeSilhouetteNormal(silhouettes[i.referenceIndex]);
        }

        return true;
    }

    return false;
}

} // namespace fcpw