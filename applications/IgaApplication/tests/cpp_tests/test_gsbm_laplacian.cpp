//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes
#include <array>
#include <cmath>
#include <vector>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "iga_application_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/kratos_components.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "spaces/default_spaces.h"

namespace Kratos::Testing
{

namespace
{

using IndexType = std::size_t;

/// Appends a closed polygonal loop of LineCondition2D2N to the skin
void CreateSkinLoop(
    ModelPart& rSkinModelPart,
    const std::vector<std::array<double, 2>>& rVertices)
{
    if (!rSkinModelPart.HasProperties(0)) {
        rSkinModelPart.CreateNewProperties(0);
    }
    Properties::Pointer p_properties = rSkinModelPart.pGetProperties(0);

    const IndexType first_node_id = rSkinModelPart.NumberOfNodes() + 1;
    const IndexType first_condition_id = rSkinModelPart.NumberOfConditions() + 1;
    const IndexType number_of_vertices = rVertices.size();
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        rSkinModelPart.CreateNewNode(
            first_node_id + i, rVertices[i][0], rVertices[i][1], 0.0
        );
    }
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        const std::vector<IndexType> node_ids{
            first_node_id + i, first_node_id + (i + 1) % number_of_vertices
        };
        rSkinModelPart.CreateNewCondition(
            "LineCondition2D2N", first_condition_id + i, node_ids, p_properties
        );
    }
}

} // namespace

// Cut patch test: the domain is the rectangle [0.25, 0.8] x [0.2, 0.8] on the
// patch [0, 1]^2 with 8 x 8 knot spans of size 0.125 and degrees (2, 3). The
// volume fractions are 0, 0.16, 0.4 and 1, with alpha = 0 on the border spans,
// and the facets cover every case: inside (both conditions off), void (ghost
// penalty only), cut (both on, d != 0) and the knot-aligned side x = 0.25
// (Dirichlet only, d = 0). The quadratic u = 1 + x^2 - 2y^2 + 3xy, with
// f = -Lap(u) = 2, lies in the spline space, satisfies the facet terms for any
// alpha, and the shifted boundary condition exactly, since the Taylor
// expansion has order p >= 2. It is thus reproduced on the whole patch.
KRATOS_TEST_CASE_IN_SUITE(GsbmLaplacianCutPatchTest, KratosIgaFastSuite)
{
    using SparseSpaceType = TDefaultSparseSpace<double>;
    using LocalSpaceType = TDefaultDenseSpace<double>;
    using LinearSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>;
    using SchemeType = ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>;
    using StrategyType = ResidualBasedLinearStrategy<
        SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>
    >;

    Model model;
    ModelPart& r_iga_model_part = model.CreateModelPart("IgaModelPart");
    r_iga_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_iga_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
    r_iga_model_part.AddNodalSolutionStepVariable(HEAT_FLUX);

    auto p_settings = Kratos::make_shared<ConvectionDiffusionSettings>();
    p_settings->SetUnknownVariable(TEMPERATURE);
    p_settings->SetDiffusionVariable(CONDUCTIVITY);
    p_settings->SetVolumeSourceVariable(HEAT_FLUX);
    p_settings->SetGradientVariable(TEMPERATURE_GRADIENT);
    r_iga_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_settings);

    NurbsGeometryModeler nurbs_modeler(model, Parameters(R"({
        "model_part_name"     : "IgaModelPart",
        "geometry_name"       : "background",
        "lower_point_xyz"     : [0.0, 0.0, 0.0],
        "upper_point_xyz"     : [1.0, 1.0, 0.0],
        "lower_point_uvw"     : [0.0, 0.0, 0.0],
        "upper_point_uvw"     : [1.0, 1.0, 0.0],
        "polynomial_order"    : [2, 3],
        "number_of_knot_spans": [8, 8]
    })"));
    nurbs_modeler.SetupGeometryModel();

    // Elements, Dirichlet conditions on the skeleton, and ghost penalty
    // conditions sharing its facets
    const Modeler& r_registered_modeler = KratosComponents<Modeler>::Get("IgaModelerGsbm");
    const auto p_gsbm_modeler = r_registered_modeler.Create(model, Parameters(R"({
        "analysis_model_part_name": "IgaModelPart",
        "geometry_name"           : "background",
        "element_condition_list"  : [{
            "type"                            : "element",
            "name"                            : "GsbmLaplacianElement",
            "iga_model_part"                  : "Domain",
            "shape_function_derivatives_order": 2
        }, {
            "type"                            : "condition",
            "name"                            : "GsbmLaplacianConditionDirichlet",
            "iga_model_part"                  : "Dirichlet",
            "shape_function_derivatives_order": 4,
            "gsbm_parameters"                 : {}
        }, {
            "type"                            : "condition",
            "name"                            : "GsbmLaplacianConditionGhostPenalty",
            "iga_model_part"                  : "GhostPenalty",
            "gsbm_parameters"                 : {"skeleton_model_part": "Dirichlet"}
        }]
    })"));
    p_gsbm_modeler->SetupModelPart();

    CreateSkinLoop(
        model.CreateModelPart("Skin"),
        {{0.25, 0.2}, {0.8, 0.2}, {0.8, 0.8}, {0.25, 0.8}}
    );
    VolumeFractionGsbmProcess gsbm_process(model, Parameters(R"({
        "model_part_name"                   : "IgaModelPart",
        "geometry_name"                     : "background",
        "elements_model_part_name"          : "IgaModelPart.Domain",
        "skeleton_model_part_name"          : "IgaModelPart.Dirichlet",
        "skin_model_part_outer_initial_name": "Skin"
    })"));
    gsbm_process.ExecuteInitialize();

    // Every facet case occurs: inside, void, cut, and a full jump of alpha
    std::array<IndexType, 4> number_of_facets{0, 0, 0, 0};
    for (const auto& r_facet : model.GetModelPart("IgaModelPart.Dirichlet").Conditions()) {
        const auto& r_geometry = r_facet.GetGeometry();
        const double alpha_minus = r_geometry.GetGeometryPart(0).GetValue(VOLUME_FRACTION);
        const double alpha_plus = r_geometry.GetGeometryPart(1).GetValue(VOLUME_FRACTION);
        if (alpha_minus == 1.0 && alpha_plus == 1.0) {
            ++number_of_facets[0];
        } else if (alpha_minus == 0.0 && alpha_plus == 0.0) {
            ++number_of_facets[1];
        } else if (std::abs(alpha_minus - alpha_plus) < 1.0) {
            ++number_of_facets[2];
        } else {
            ++number_of_facets[3];
        }
    }
    for (const IndexType number : number_of_facets) {
        KRATOS_EXPECT_GT(number, 0);
    }

    const auto exact_solution = [](const array_1d<double, 3>& rPoint) {
        const double x = rPoint[0];
        const double y = rPoint[1];
        return 1.0 + x * x - 2.0 * y * y + 3.0 * x * y;
    };

    // Unit conductivity, the paper's chi = 20 and kappa = 0.5, source f = 2,
    // and u_D at the closest skin points
    auto p_properties = r_iga_model_part.CreateNewProperties(1);
    p_properties->SetValue(CONDUCTIVITY, 1.0);
    p_properties->SetValue(PENALTY_FACTOR, 20.0);
    p_properties->SetValue(GHOST_PENALTY_FACTOR, 0.5);
    for (auto& r_element : r_iga_model_part.Elements()) {
        r_element.SetProperties(p_properties);
        r_element.SetValue(HEAT_FLUX, 2.0);
    }
    for (auto& r_condition : r_iga_model_part.Conditions()) {
        r_condition.SetProperties(p_properties);
    }
    for (auto& r_node : model.GetModelPart("Skin.projection_data").Nodes()) {
        r_node.SetValue(TEMPERATURE, exact_solution(r_node.Coordinates()));
    }
    for (auto& r_node : r_iga_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }

    StrategyType strategy(
        r_iga_model_part,
        Kratos::make_shared<SchemeType>(),
        Kratos::make_shared<LinearSolverType>()
    );
    strategy.Solve();

    // Solution at every Gauss point of the patch, inside and outside the domain
    const double tolerance = 1.0e-10;
    for (const auto& r_element : r_iga_model_part.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        double solution = 0.0;
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            solution += r_N(0, i) * r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
        }
        KRATOS_EXPECT_NEAR(solution, exact_solution(r_geometry.Center().Coordinates()), tolerance);
    }
}

} // namespace Kratos::Testing
