// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_elements/transient_Pw_element.h"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "geometries/tetrahedra_3d_4.h"
#include "test_setup_utilities/element_setup_utilities.hpp"

#include <benchmark/benchmark.h>

namespace
{

using namespace Kratos;

std::shared_ptr<Properties> CreatePropertiesForTransientPwElementBenchmark()
{
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    p_properties->SetValue(POROSITY, 1.000000e-01);
    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    p_properties->SetValue(BULK_MODULUS_FLUID, 200.0);
    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_YY, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_ZZ, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_XY, 0.000000e+00);
    p_properties->SetValue(PERMEABILITY_YZ, 0.000000e+00);
    p_properties->SetValue(PERMEABILITY_ZX, 0.000000e+00);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_properties->SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);
    p_properties->SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED");

    return p_properties;
}

auto CreateTransientPwElement3D4NForBenchmark(const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 1.0, 1.0, 1.0));

    auto p_element = make_intrusive<TransientPwElement<3, 4>>(
        1, std::make_shared<Tetrahedra3D4<Node>>(nodes), rProperties,
        std::make_unique<ThreeDimensionalStressState>(), nullptr);

    const auto solution_step_variables = Geo::ConstVariableDataRefs{
        std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)};
    Testing::ElementSetupUtilities::AddVariablesToNodes(nodes, solution_step_variables,
                                                         Geo::ConstVariableRefs{std::cref(WATER_PRESSURE)});
    for (auto& r_node : nodes) {
        r_node.SetBufferSize(2);
    }

    return p_element;
}

void SetTransientPwElement3D4NValues(const Element::Pointer& rElement)
{
    const auto gravity = array_1d<double, 3>{0.0, -10.0, 0.0};
    for (int counter = 0; auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = counter * 1.0e5;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE)   = counter * 5.0e5;
        ++counter;
    }
}

} // namespace

namespace Kratos
{

void benchmarkTransientPwElement3D4NLocalSystemCalculation(benchmark::State& rState)
{
    const auto p_properties = CreatePropertiesForTransientPwElementBenchmark();
    auto       p_element    = CreateTransientPwElement3D4NForBenchmark(p_properties);
    SetTransientPwElement3D4NValues(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    for (auto _ : rState) {
        auto left_hand_side  = Matrix{};
        auto right_hand_side = Vector{};
        p_element->CalculateLocalSystem(left_hand_side, right_hand_side, dummy_process_info);
        benchmark::DoNotOptimize(left_hand_side);
        benchmark::DoNotOptimize(right_hand_side);
    }
}

BENCHMARK(benchmarkTransientPwElement3D4NLocalSystemCalculation);

} // namespace Kratos

BENCHMARK_MAIN();
