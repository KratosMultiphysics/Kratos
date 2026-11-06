//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo Crescenzio
//

// System includes

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "geometries/quadrilateral_2d_4.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "utilities/quadrature_points_utility.h"

// Application includes
#include "mpm_application_variables.h"
#include "custom_utilities/brute_force_material_point_locator.h"
#include "custom_elements/mpm_updated_lagrangian.hpp"

namespace Kratos
{

static void BM_MPElementLocator(benchmark::State& state) {

    const int mp_per_dimension = state.range(0);

    const double xmin = -1.0;
    const double xmax = +1.0;
    const double ymin = -1.0;
    const double ymax = +1.0;

    const double dx = (xmax-xmin)/(mp_per_dimension+1);
    const double dy = (ymax-ymin)/(mp_per_dimension+1);

    Model model;
    ModelPart& mpm_model_part = model.CreateModelPart("MPMModelPart");

    auto p_node_1 = mpm_model_part.CreateNewNode(1, xmin, ymin, 0.0);
    auto p_node_2 = mpm_model_part.CreateNewNode(2, xmax, ymin, 0.0);
    auto p_node_3 = mpm_model_part.CreateNewNode(3, xmin, ymax, 0.0);
    auto p_node_4 = mpm_model_part.CreateNewNode(4, xmax, ymax, 0.0);
    auto p_geometry = Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_1, p_node_2, p_node_4, p_node_3);

    const auto& process_info = mpm_model_part.GetProcessInfo();
    auto p_prop = mpm_model_part.CreateNewProperties(0);

    for (int i=1; i<=mp_per_dimension; ++i) {
        for (int j=1; j<=mp_per_dimension; ++j) {
            double xcoord = xmin + dx*i;
            double ycoord = ymin + dy*j;
            array_1d<double, 3> mp_coord{ xcoord, ycoord, 0.0 };
            auto p_quad = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(p_geometry, mp_coord, 1.0);
            auto p_elem = Kratos::make_intrusive<MPMUpdatedLagrangian>(mp_per_dimension*(i-1)+j, p_quad, p_prop);
            mpm_model_part.AddElement(p_elem);
            p_elem->SetValuesOnIntegrationPoints(MP_COORD, { mp_coord }, process_info);
        }
    }

    auto locator = BruteForceMaterialPointLocator(mpm_model_part);

    const Point point{ 0.0, 0.0, 0.0 };
    constexpr double tolerance = 1e-6;
    // const int expected_id = (mp_per_dimension*mp_per_dimension-1)/2 + 1;

    for (auto _ : state) {
        locator.FindElement(point, tolerance);
    }
}

// Register benchmarks and provide input size as command-line option
BENCHMARK(BM_MPElementLocator)->Arg(11)->Arg(101)->Arg(1001);

}  // namespace Kratos.

BENCHMARK_MAIN();
