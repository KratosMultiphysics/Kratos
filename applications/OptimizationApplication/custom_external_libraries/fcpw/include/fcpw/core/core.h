#pragma once

// global includes
#define _USE_MATH_DEFINES
#define NOMINMAX
#include <algorithm>
#include <math.h>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <type_traits>
#ifdef FCPW_USE_ENOKI
    #include <enoki/array.h>
#endif
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace fcpw {

#ifdef FCPW_USE_ENOKI
    template<size_t WIDTH>
    using IntP = enoki::Packet<int, WIDTH>;
    template<size_t WIDTH>
    using FloatP = enoki::Packet<float, WIDTH>;
    template<size_t WIDTH>
    using MaskP = enoki::mask_t<FloatP<WIDTH>>;
    template<size_t DIM>
    using enokiVector = enoki::Array<float, DIM>;
    using enokiVector2 = enoki::Array<float, 2>;
    using enokiVector3 = enoki::Array<float, 3>;
    template<size_t WIDTH, size_t DIM>
    using VectorP = enoki::Array<FloatP<WIDTH>, DIM>;
    template<size_t WIDTH>
    using Vector1P = VectorP<WIDTH, 1>;
    template<size_t WIDTH>
    using Vector2P = VectorP<WIDTH, 2>;
    template<size_t WIDTH>
    using Vector3P = VectorP<WIDTH, 3>;
#endif

template<size_t DIM>
using Vector = Eigen::Matrix<float, DIM, 1>;
using Vector1 = Vector<1>;
using Vector2 = Vector<2>;
using Vector3 = Vector<3>;

using Vector2i = Eigen::Vector2i;
using Vector3i = Eigen::Vector3i;

template<size_t DIM>
using Transform = Eigen::Transform<float, DIM, Eigen::Affine>;

template<size_t DIM>
struct Ray;
template<size_t DIM>
struct BoundingSphere;
template<size_t DIM>
struct BoundingBox;
template<size_t DIM>
struct BoundingCone;
template<size_t DIM>
struct Interaction;
template<size_t DIM>
class Primitive;
template<size_t DIM>
class GeometricPrimitive;
template<size_t DIM>
class SilhouettePrimitive;
template<size_t DIM>
class Aggregate;

static const float minFloat = std::numeric_limits<float>::lowest();
static const float maxFloat = std::numeric_limits<float>::max();
static const int minInt = std::numeric_limits<int>::min();
static const int maxInt = std::numeric_limits<int>::max();
static const float epsilon = std::numeric_limits<float>::epsilon();
static const float oneMinusEpsilon = 1.0f - epsilon;

inline bool inRange(float val, float low, float high) {
    return val >= low && val <= high;
}

inline float uniformRealRandomNumber(float a=0.0f, float b=1.0f)
{
    thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<float> distribution(a, b);

    return distribution(generator);
}

template<size_t DIM>
inline Vector<DIM> uniformRealRandomVector(float a=0.0f, float b=1.0f)
{
    Vector<DIM> v;
    for (size_t i = 0; i < DIM; i++) {
        v[i] = uniformRealRandomNumber(a, b);
    }

    return v;
}

template<size_t DIM>
inline Vector<DIM> rotate(const Vector<DIM>& u, const Vector<DIM>& v, float theta)
{
    std::cerr << "rotate<DIM>() not implemented" << std::endl;
    exit(EXIT_FAILURE);

    return Vector<DIM>::Zero();
}

template<>
inline Vector2 rotate<2>(const Vector2& u, const Vector2& v, float theta)
{
    float det = u[0]*v[1] - u[1]*v[0];
    theta *= std::copysign(1.0f, det);
    float cosTheta = std::cos(theta);
    float sinTheta = std::sin(theta);

    return Vector2(cosTheta*u[0] - sinTheta*u[1],
                   sinTheta*u[0] + cosTheta*u[1]);
}

template<>
inline Vector3 rotate<3>(const Vector3& u, const Vector3& v, float theta)
{
    float cosTheta = std::cos(theta);
    float sinTheta = std::sin(theta);
    Vector3 w = u.cross(v).normalized();
    Vector3 oneMinusCosThetaW = (1.0f - cosTheta)*w;

    Eigen::Matrix3f R; 
    R << cosTheta + oneMinusCosThetaW[0]*w[0],
         oneMinusCosThetaW[1]*w[0] - sinTheta*w[2],
         oneMinusCosThetaW[2]*w[0] + sinTheta*w[1],
         oneMinusCosThetaW[0]*w[1] + sinTheta*w[2],
         cosTheta + oneMinusCosThetaW[1]*w[1],
         oneMinusCosThetaW[2]*w[1] - sinTheta*w[0],
         oneMinusCosThetaW[0]*w[2] - sinTheta*w[1],
         oneMinusCosThetaW[1]*w[2] + sinTheta*w[0],
         cosTheta + oneMinusCosThetaW[2]*w[2];

    return R*u;
}

} // namespace fcpw