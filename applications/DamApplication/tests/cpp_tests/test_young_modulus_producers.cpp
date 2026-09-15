// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Migration contract of the Young's modulus producers to the standard nodal
// YOUNG_MODULUS + DatabaseAccessor mechanism. Each test verifies the full path
//
//     producer
//         -> nodal historical YOUNG_MODULUS
//         -> DatabaseAccessor
//         -> Properties::GetValue(YOUNG_MODULUS, geometry, N, process_info)
//         -> expected integration-point Young's modulus
//
// while keeping the historical "NODAL_YOUNG_MODULUS" variable_name accepted for
// backward compatibility.
//
#include <array>
#include <cmath>
#include <string>
#include <vector>

// Project includes
#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/table.h"
#include "includes/variables.h"
#include "geometries/line_2d_2.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_utilities/nodal_young_modulus_utilities.h"
#include "custom_processes/dam_random_fields_variable_process.hpp"
#include "custom_processes/dam_chemo_mechanical_aging_young_process.hpp"
#include "custom_processes/dam_azenha_heat_source_process.hpp"

namespace Kratos
{
namespace Testing
{
namespace
{

/// Builds a small two-node model part carrying the nodal field variables and a
/// single property, and returns the geometry used to query the accessor path.
ModelPart& BuildYoungModulusModel(Model& rModel, Geometry<Node>::PointsArrayType& rPoints)
{
    ModelPart& r_mp = rModel.CreateModelPart("M", 2);
    for (auto& v : std::vector<const VariableData*>{&YOUNG_MODULUS, &NODAL_YOUNG_MODULUS,
            &ALPHA_HEAT_SOURCE, &TEMPERATURE, &HEAT_FLUX}) {
        r_mp.AddNodalSolutionStepVariable(*v);
    }
    r_mp.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_mp.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_prop = r_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = 2.0e7; // constant property fallback
    (*p_prop)[POISSON_RATIO] = 0.2;
    rPoints.clear();
    rPoints.push_back(r_mp.pGetNode(1));
    rPoints.push_back(r_mp.pGetNode(2));
    return r_mp;
}

/// Evaluates nodal YOUNG_MODULUS through the Properties::GetValue accessor path
/// with uniform shape functions (midpoint of a two-node line).
double InterpolatedYoungViaAccessor(Properties& rProps, Geometry<Node>& rGeometry, const ProcessInfo& rPi)
{
    std::vector<double> n_data = {0.5, 0.5};
    Vector N(n_data.size());
    for (std::size_t i = 0; i < n_data.size(); ++i) N[i] = n_data[i];
    return rProps.GetValue(YOUNG_MODULUS, rGeometry, N, rPi);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RandomFieldsLegacyNameProducesNodalYoungModulus, KratosDamFastSuite)
{
    Model model;
    Geometry<Node>::PointsArrayType pts;
    ModelPart& r_mp = BuildYoungModulusModel(model, pts);

    // Per-node-ID field (as produced by the gstools wrappers).
    Table<double, double> table;
    table.PushBack(1.0, 10.0);
    table.PushBack(2.0, 20.0);

    Parameters params(R"({"model_part_name":"M","variable_name":"NODAL_YOUNG_MODULUS"})");
    DamRandomFieldsVariableProcess process(r_mp, table, params);
    process.ExecuteBeforeSolutionLoop();

    // The legacy caller produced nodal historical YOUNG_MODULUS...
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(1)->FastGetSolutionStepValue(YOUNG_MODULUS), 10.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(2)->FastGetSolutionStepValue(YOUNG_MODULUS), 20.0, 1.0e-12);
    // ... and must not touch the legacy variable.
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(1)->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS), 0.0, 1.0e-12);

    // The DatabaseAccessor is installed and the interpolation reaches the law.
    auto& r_props = *r_mp.pGetProperties(1);
    KRATOS_EXPECT_TRUE(r_props.HasAccessor(YOUNG_MODULUS));
    KRATOS_EXPECT_NEAR(InterpolatedYoungViaAccessor(r_props, *Geometry<Node>::Pointer(new Line2D2<Node>(pts)), r_mp.GetProcessInfo()), 15.0, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(ChemoMechanicalAgingLegacyNameProducesNodalYoungModulus, KratosDamFastSuite)
{
    Model model;
    Geometry<Node>::PointsArrayType pts;
    ModelPart& r_mp = BuildYoungModulusModel(model, pts);

    r_mp.GetProcessInfo()[TIME] = 10.0 * 31536000.0; // 10 years in seconds

    Parameters params(R"({
        "model_part_name" : "M",
        "variable_name" : "NODAL_YOUNG_MODULUS",
        "initial_elastic_modulus" : 30.0e9,
        "initial_porosity" : 0.2,
        "max_chemical_porosity" : 0.32,
        "chemical_characteristic_aging_time" : 100.0,
        "max_mechanical_damage" : 0.32,
        "damage_characteristic_aging_time" : 100.0
    })");
    DamChemoMechanicalAgingYoungProcess process(r_mp, params);
    process.ExecuteBeforeSolutionLoop();

    const double time_years = 10.0;
    const double sound_concrete = 30.0e9 * std::sqrt(1.0 + 0.0805 * std::log(time_years));
    const double chemical_porosity = 0.32 * (1.0 - std::exp(-time_years / 100.0));
    const double damage_mechanical = 0.32 * (1.0 - std::exp(-time_years / 100.0));
    const double expected = ((1.0 - 0.2 - chemical_porosity) * (1.0 - damage_mechanical) * sound_concrete) / (1.0 - 0.2);

    KRATOS_EXPECT_NEAR(r_mp.pGetNode(1)->FastGetSolutionStepValue(YOUNG_MODULUS), expected, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(2)->FastGetSolutionStepValue(YOUNG_MODULUS), expected, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(1)->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS), 0.0, 1.0e-12);

    auto& r_props = *r_mp.pGetProperties(1);
    KRATOS_EXPECT_TRUE(r_props.HasAccessor(YOUNG_MODULUS));
    KRATOS_EXPECT_NEAR(InterpolatedYoungViaAccessor(r_props, *Geometry<Node>::Pointer(new Line2D2<Node>(pts)), r_mp.GetProcessInfo()), expected, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(AzenhaHeatFluxAgingWritesNodalYoungModulus, KratosDamFastSuite)
{
    Model model;
    Geometry<Node>::PointsArrayType pts;
    ModelPart& r_mp = BuildYoungModulusModel(model, pts);

    // Aging branch: ExecuteBeforeSolutionLoop -> ExecuteInitializeAging writes
    // E(alpha_initial) = sqrt(alpha_initial)*E_inf.
    Parameters params(R"({
        "model_part_name" : "M",
        "variable_name" : "HEAT_FLUX",
        "activation_energy" : 100.0,
        "gas_constant" : 8.314,
        "constant_rate" : 1.0,
        "alpha_initial" : 0.5,
        "q_total" : 1.0,
        "aging" : true,
        "young_inf" : 2.0e8,
        "A" : 0.0,
        "B" : 0.0,
        "C" : 0.0,
        "D" : 0.0
    })");
    DamAzenhaHeatFluxProcess process(r_mp, params);
    process.ExecuteBeforeSolutionLoop();

    const double expected = std::sqrt(0.5) * 2.0e8;
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(1)->FastGetSolutionStepValue(YOUNG_MODULUS), expected, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(2)->FastGetSolutionStepValue(YOUNG_MODULUS), expected, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_mp.pGetNode(1)->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS), 0.0, 1.0e-12);

    auto& r_props = *r_mp.pGetProperties(1);
    KRATOS_EXPECT_TRUE(r_props.HasAccessor(YOUNG_MODULUS));
    KRATOS_EXPECT_NEAR(InterpolatedYoungViaAccessor(r_props, *Geometry<Node>::Pointer(new Line2D2<Node>(pts)), r_mp.GetProcessInfo()), expected, 1.0e-12);
}

} // namespace Testing
} // namespace Kratos