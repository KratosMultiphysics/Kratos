// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "test_setup_utilities/element_setup_utilities.hpp"

#include <benchmark/benchmark.h>

namespace
{

using namespace Kratos;

std::shared_ptr<Properties> CreatePropertiesForUPwDiffOrderElementBenchmark()
{
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<GeoIncrementalLinearElasticLaw>(
                                                 std::make_unique<PlaneStrain>()));
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    p_properties->SetValue(DENSITY_SOLID, 2.650000e+03);
    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    p_properties->SetValue(POROSITY, 1.000000e-01);
    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    p_properties->SetValue(BULK_MODULUS_FLUID, 200.0); // small to get a significant value for the compressibility term
    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_YY, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_XY, 0.000000e+00);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_properties->SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);
    p_properties->SetValue(IGNORE_UNDRAINED, false);

    return p_properties;
}

void SetSolutionStepValuesForGeneralCheck(const Element::Pointer& rElement)
{
    const auto zero_values = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto gravity     = array_1d<double, 3>{0.0, -9.81, 0.0};

    for (int counter = 0; auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(VELOCITY)            = zero_values;
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = counter * 1.0e5;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE)   = counter * 5.0e5;
        ++counter;
    }
    rElement->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{-0.015, 0.0, 0.0};
    rElement->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.015, 0.00, 0.0};
    rElement->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.015, 0.0};
}

auto CreateSmallStrainUPwDiffOrderElementWithUPwDofs(const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 0.0, -1.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 0.0, -0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(5, 0.5, -0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(6, 0.5, 0.05, 0.0));

    auto result = Testing::ElementSetupUtilities::Create2D6NDiffOrderElement(nodes, rProperties);
    const auto solution_step_variables = Geo::ConstVariableDataRefs{
        std::cref(WATER_PRESSURE),     std::cref(DT_WATER_PRESSURE), std::cref(DISPLACEMENT),
        std::cref(VELOCITY),           std::cref(ACCELERATION),      std::cref(VOLUME_ACCELERATION),
        std::cref(HYDRAULIC_DISCHARGE)};
    const auto degrees_of_freedom =
        Geo::ConstVariableRefs{std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT_X),
                               std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_Z)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(result, solution_step_variables, degrees_of_freedom);

    for (auto& r_node : nodes) {
        r_node.SetBufferSize(2);
    }
    return result;
}

} // namespace

namespace Kratos
{

void benchmarkUPwDiffOrderLocalSystemCalculation(benchmark::State& rState)
{
    const auto p_properties = CreatePropertiesForUPwDiffOrderElementBenchmark();
    auto       p_element    = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    for (auto _ : rState) {
        auto left_hand_side  = Matrix{};
        auto right_hand_side = Vector{};
        p_element->CalculateLocalSystem(left_hand_side, right_hand_side, dummy_process_info);
    }
}

void benchmarkUPwDiffOrderRHSCalculation(benchmark::State& rState)
{
    const auto p_properties = CreatePropertiesForUPwDiffOrderElementBenchmark();
    auto       p_element    = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    for (auto _ : rState) {
        auto right_hand_side = Vector{};
        p_element->CalculateRightHandSide(right_hand_side, dummy_process_info);
    }
}

void benchmarkUPwDiffOrderLHSCalculation(benchmark::State& rState)
{
    const auto p_properties = CreatePropertiesForUPwDiffOrderElementBenchmark();
    auto       p_element    = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    for (auto _ : rState) {
        auto left_hand_side = Matrix{};
        p_element->CalculateLeftHandSide(left_hand_side, dummy_process_info);
    }
}

BENCHMARK(benchmarkUPwDiffOrderLocalSystemCalculation);
BENCHMARK(benchmarkUPwDiffOrderRHSCalculation);
BENCHMARK(benchmarkUPwDiffOrderLHSCalculation);

} // namespace Kratos

BENCHMARK_MAIN();