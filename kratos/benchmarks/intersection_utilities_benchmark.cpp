//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "geometries/point.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{

// Sample data for benchmarking
Point::Pointer p_point_1(make_shared<Point>( 0.0, 0.0, 0.0));
Point::Pointer p_point_2(make_shared<Point>( 1.0, 0.0, 0.0));
Point::Pointer p_point_3(make_shared<Point>( 0.0, 1.0, 0.0));

Triangle3D3<Point> triangle(p_point_1, p_point_3, p_point_2);
Point line_point_1(1.0, 0.25, 0.25);
Point line_point_2(-1.0, 0.25, 0.25);
Point intersection_point(0.0, 0.0, 0.0);
Point triangle_test_point(0.2, 0.2, 0.0);

static void BM_ComputeTriangleLineIntersection(benchmark::State& state) {
    for (auto _ : state) {
        const int intersection_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triangle,
            line_point_1.Coordinates(),
            line_point_2.Coordinates(),
            intersection_point.Coordinates());

        benchmark::DoNotOptimize(intersection_id);
        benchmark::DoNotOptimize(intersection_point);
    }
}

static void BM_PointInTriangle(benchmark::State& state) {
    for (auto _ : state) {
        const bool is_inside = IntersectionUtilities::PointInTriangle(
            triangle[0].Coordinates(),
            triangle[1].Coordinates(),
            triangle[2].Coordinates(),
            triangle_test_point.Coordinates());

        benchmark::DoNotOptimize(is_inside);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_ComputeTriangleLineIntersection);
BENCHMARK(BM_PointInTriangle);

}  // namespace Kratos

BENCHMARK_MAIN();