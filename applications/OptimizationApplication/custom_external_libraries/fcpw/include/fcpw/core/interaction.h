#pragma once

#include <fcpw/core/core.h>

namespace fcpw {

enum class DistanceInfo {
    Exact,
    Bounded
};

template<size_t DIM>
struct Interaction {
    // constructor
    Interaction(): d(maxFloat), sign(0), primitiveIndex(-1), nodeIndex(-1), referenceIndex(-1),
                   objectIndex(-1), p(Vector<DIM>::Zero()), n(Vector<DIM>::Zero()),
                   uv(Vector<DIM - 1>::Zero()), distanceInfo(DistanceInfo::Exact) {}

    // comparison operators
    bool operator==(const Interaction<DIM>& i) const {
        bool distancesMatch = std::fabs(d - i.d) < 1e-6f;
        if (distanceInfo == DistanceInfo::Bounded) return distancesMatch;

        return distancesMatch && (p - i.p).squaredNorm() < 1e-6f;
    }

    bool operator!=(const Interaction<DIM>& i) const {
        return !(*this == i);
    }

    // returns signed distance
    float signedDistance(const Vector<DIM>& x) const {
        return sign == 0 ? ((x - p).dot(n) > 0.0f ? 1.0f : -1.0f)*d : sign*d;
    }

    // samples a random point on geometric primitive and returns sampling pdf
    float samplePoint(const Primitive<DIM> *primitive, const Vector<DIM>& randNums) {
        return static_cast<const GeometricPrimitive<DIM> *>(primitive)->samplePoint(randNums, uv, p, n);
    }

    // computes normal from geometric primitive if unspecified
    void computeNormal(const Primitive<DIM> *primitive) {
        n = static_cast<const GeometricPrimitive<DIM> *>(primitive)->normal(uv);
    }

    // computes silhouette normal
    void computeSilhouetteNormal(const SilhouettePrimitive<DIM> *silhouette) {
        n = silhouette->normal();
    }

    // applies transform
    void applyTransform(const Transform<DIM>& t,
                        const Transform<DIM>& tInv,
                        const Vector<DIM>& query,
                        bool overwriteDistance = true) {
        p = t*p;
        if (overwriteDistance) d = (p - query).norm();
        n = Transform<DIM>(tInv.matrix().transpose())*n;
        n.normalize();
    }

    // members
    float d;
    int sign; // sign bit used for difference ops
    int primitiveIndex; // index of primitive in polygon soup
    int nodeIndex; // index of aggregate node containing intersected or closest point
    int referenceIndex; // reference index of primitive for internal library use
    int objectIndex; // index of object containing the primitive
    Vector<DIM> p, n;
    Vector<DIM - 1> uv;
    DistanceInfo distanceInfo;
};

template<size_t DIM>
inline bool compareInteractions(const Interaction<DIM>& i, const Interaction<DIM>& j) {
    return i.d < j.d;
}

template<size_t DIM>
inline std::vector<Interaction<DIM>> removeDuplicates(const std::vector<Interaction<DIM>>& is) {
    int N = (int)is.size();
    std::vector<bool> isDuplicate(N, false);
    std::vector<Interaction<DIM>> cs;

    for (int i = 0; i < N - 1; i++) {
        if (is[i] == is[i + 1]) isDuplicate[i + 1] = true;
    }

    for (int i = 0; i < N; i++) {
        if (!isDuplicate[i]) cs.emplace_back(is[i]);
    }

    return cs;
}

} // namespace fcpw