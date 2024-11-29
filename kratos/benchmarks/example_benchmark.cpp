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
#include "utilities/geometry_utilities/nearest_point_utilities.h"

namespace Kratos
{

// Sample data for benchmarking
Point::Pointer p_point_1(make_shared<Point>( 0.0, 0.0, 0.0));
Point::Pointer p_point_2(make_shared<Point>( 1.0, 0.0, 0.0));
Point::Pointer p_point_3(make_shared<Point>( 0.0, 1.0, 0.0));

Triangle3D3<Point> triangle(p_point_1, p_point_3, p_point_2);
Point nearest_point(0.0, 0.0, 0.0);
Point inside_point(0.2, 0.1, 0.00);

static void BM_TriangleNearestPoint(benchmark::State& state) {
    for (auto _ : state) {
        NearestPointUtilities::TriangleNearestPoint(inside_point, triangle, nearest_point);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_TriangleNearestPoint);

}  // namespace Kratos

BENCHMARK_MAIN();