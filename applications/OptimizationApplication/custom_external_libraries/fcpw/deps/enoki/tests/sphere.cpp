/*
    tests/sphere.cpp -- a simple sphere ray tracer

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#if defined(_MSC_VER)
#  pragma warning(disable: 4723) /// potential divide by 0
#endif

#include <enoki/dynamic.h>
#include <fstream>
#include "ray.h"
#include <chrono>

// -----------------------------------------------------------------------
//! @{ \name Convenient type aliases
// -----------------------------------------------------------------------

/* Floats and packets of floats */
using FloatP    = Array<float>;
using FloatX    = DynamicArray<FloatP>;

/* 2D vectors and static/dynamic packets of 2D vectors */
using Vector2f  = Array<float, 2>;
using Vector2fP = Array<FloatP, 2>;
using Vector2fX = Array<FloatX, 2>;

/* 3D vectors and static/dynamic packets of 3D vectors */
using Vector3f  = Array<float, 3>;
using Vector3fP = Array<FloatP, 3>;
using Vector3fX = Array<FloatX, 3>;

/* rays and static/dynamic packets of rays */
using Ray3f     = Ray<Vector3f>;
using Ray3fP    = Ray<Vector3fP>;
using Ray3fX    = Ray<Vector3fX>;

/* Aliases to create types that are compatible with other type */
template <typename T> using vector3f_t = Array<value_t<T>, 3>;
template <typename T> using ray3f_t    = Ray<vector3f_t<T>>;

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Computational kernels of the ray tracer
// -----------------------------------------------------------------------

/// "Sensor": trace rays for a given X and Y coordinate
template <typename Vector2f> ray3f_t<Vector2f> make_rays(const Vector2f &p) {
    using Vector3f = vector3f_t<Vector2f>;
    using Ray3f = ray3f_t<Vector2f>;

    return Ray3f(Vector3f(p.x(), p.y(), -1.f),
                 Vector3f(0.f, 0.f, 1.f));
}

/// "Shape": intersect against sphere
template <typename Ray> ENOKI_INLINE typename Ray::Vector intersect_rays(const Ray &r) {
    /* Coefficients of quadratic */
    auto a = dot(r.d, r.d);
    auto b = 2.f * dot(r.o, r.d);
    auto c = dot(r.o, r.o) - 1.f;

    /* Solve quadratic equation */
    auto discrim = b*b - 4.f*a*c;
    auto t = (-b + sqrt(discrim)) / (2.f * a);

    return select(discrim >= 0.f, r(t), 0.f);
}

/// "Shader": directional illumination
template <typename Vector3f> ENOKI_INLINE typename Vector3f::Value shade_hits(Vector3f n) {
    return 0.2f + max(dot(n, Vector3f(-1.f, -1.f, 2.f)), 0.f) * 90.f;
}

/// All three kernels combined into one
template <typename Vector2> ENOKI_INLINE typename Vector2::Value combined(Vector2 n) {
    return shade_hits(intersect_rays(make_rays(n)));
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Wrappers which execute the above kernels for dynamic arrays
// -----------------------------------------------------------------------

Ray3fX make_rays_dynamic(const Vector2fX &p) {
    return vectorize([](auto &&p) { return make_rays<Vector2fP>(p); }, p);
}

Vector3fX intersect_rays_dynamic(const Ray3fX &r) {
    return vectorize([](auto &&r) { return intersect_rays<Ray3fP>(r); }, r);
}

FloatX shade_hits_dynamic(const Vector3fX &n) {
    return vectorize([](auto &&n) { return shade_hits<Vector3fP>(n); }, n);
}

FloatX combined_dynamic(const Vector2fX &p) {
    return vectorize([](auto &&p) { return combined<Vector2fP>(p); }, p);
}

//! @}
// -----------------------------------------------------------------------

auto clk() { return std::chrono::high_resolution_clock::now(); }

template <typename T> float clkdiff(T a, T b) {
    return std::chrono::duration<float>(b - a).count() * 1000;
}

void write_image(const std::string &filename, const FloatX &image) {
    std::ofstream os(filename);
    os << "P3\n1024 1024\n255\n";
    for (float v : image)
        os << (int) v << " " << (int) v << " " << (int) v << "\n";
}

int main(int /* argc */, char ** /* argv */) {
    auto idx = linspace<FloatX>(-1.2f, 1.2f, 1024);
    auto grid = meshgrid(idx, idx);

    /* benchmark1 */ {
        auto time_start = clk();
        Ray3fX    rays  = make_rays_dynamic(grid);
        Vector3fX hits  = intersect_rays_dynamic(rays);
        FloatX image = shade_hits_dynamic(hits);
        auto time_end = clk();
        std::cerr << "Separate kernels: " << clkdiff(time_start, time_end) << " ms " << std::endl;
        write_image("sphere1.ppm", image);
    }

    /* benchmark2 */ {
        auto time_start = clk();
        FloatX image = combined_dynamic(grid);
        auto time_end = clk();
        std::cerr << "Combined kernels: " << clkdiff(time_start, time_end) << " ms " << std::endl;
        write_image("sphere2.ppm", image);
    }

    return 0;
}
