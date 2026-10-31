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
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "geometries/nurbs_curve_geometry.h"
#include "iga_application_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/global_pointer_variables.h"
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

/// Adds a counterclockwise ellipse to the skin, as one rational quadratic
/// curve with 9 control points: the affine image of the exact circle
void CreateNurbsEllipse(
    ModelPart& rSkinModelPart,
    const array_1d<double, 3>& rCenter,
    const std::array<double, 2>& rSemiAxes,
    const std::string& rLayer)
{
    const std::array<std::array<double, 2>, 9> directions{{
        {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}
    }};
    PointerVector<Node> control_points;
    for (const auto& r_direction : directions) {
        control_points.push_back(Kratos::make_intrusive<Node>(0,
            rCenter[0] + rSemiAxes[0] * r_direction[0],
            rCenter[1] + rSemiAxes[1] * r_direction[1], 0.0
        ));
    }

    const double w = std::sqrt(2.0) / 2.0;
    const std::array<double, 10> knot_values{0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1};
    const std::array<double, 9> weight_values{1, w, 1, w, 1, w, 1, w, 1};
    Vector knots(knot_values.size());
    std::copy(knot_values.begin(), knot_values.end(), knots.begin());
    Vector weights(weight_values.size());
    std::copy(weight_values.begin(), weight_values.end(), weights.begin());

    auto p_curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Node>>>(
        control_points, 2, knots, weights
    );
    p_curve->SetId(rSkinModelPart.NumberOfGeometries());
    p_curve->SetValue(IDENTIFIER, rLayer);
    rSkinModelPart.AddGeometry(p_curve);
}

/// \partial_x^A \partial_y^B u of u = \sum_{i, j <= p} cos(1 + i + 3j) x^i y^j
double ComputePolynomialDerivative(
    const int Degree,
    const int A,
    const int B,
    const array_1d<double, 3>& rPoint)
{
    double value = 0.0;
    for (int i = A; i <= Degree; ++i) {
        for (int j = B; j <= Degree; ++j) {
            // cos(1 + i + 3j) i! / (i - A)! j! / (j - B)! x^(i-A) y^(j-B)
            double term = std::cos(1.0 + i + 3.0 * j);
            for (int l = 0; l < A; ++l) {
                term *= i - l;
            }
            for (int l = 0; l < B; ++l) {
                term *= j - l;
            }
            value += term * std::pow(rPoint[0], i - A) * std::pow(rPoint[1], j - B);
        }
    }
    return value;
}

/// The cut patch of GsbmLaplacianCutPatchTest, with its volume fractions
ModelPart& CreateCutPatchModelPart(Model& rModel)
{
    ModelPart& r_iga_model_part = rModel.CreateModelPart("IgaModelPart");
    r_iga_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_iga_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
    r_iga_model_part.AddNodalSolutionStepVariable(HEAT_FLUX);

    auto p_settings = Kratos::make_shared<ConvectionDiffusionSettings>();
    p_settings->SetUnknownVariable(TEMPERATURE);
    p_settings->SetDiffusionVariable(CONDUCTIVITY);
    p_settings->SetVolumeSourceVariable(HEAT_FLUX);
    p_settings->SetGradientVariable(TEMPERATURE_GRADIENT);
    p_settings->SetDensityVariable(DENSITY);
    p_settings->SetSpecificHeatVariable(SPECIFIC_HEAT);
    r_iga_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_settings);

    NurbsGeometryModeler nurbs_modeler(rModel, Parameters(R"({
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
    const auto p_gsbm_modeler = r_registered_modeler.Create(rModel, Parameters(R"({
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
        rModel.CreateModelPart("Skin"),
        {{0.25, 0.2}, {0.8, 0.2}, {0.8, 0.8}, {0.25, 0.8}}
    );
    VolumeFractionGsbmProcess gsbm_process(rModel, Parameters(R"({
        "model_part_name"                   : "IgaModelPart",
        "geometry_name"                     : "background",
        "elements_model_part_name"          : "IgaModelPart.Domain",
        "skeleton_model_part_name"          : "IgaModelPart.Dirichlet",
        "skin_model_part_outer_initial_name": "Skin"
    })"));
    gsbm_process.ExecuteInitialize();

    return r_iga_model_part;
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
// expansion has order p >= 2. It is thus reproduced on the whole patch, with
// the default skew-symmetric method without penalty and with the symmetric one
// with the penalty chi = 20 of the paper.
KRATOS_TEST_CASE_IN_SUITE(GsbmLaplacianCutPatchTest, KratosIgaFastSuite)
{
    using SparseSpaceType = TDefaultSparseSpace<double>;
    using LocalSpaceType = TDefaultDenseSpace<double>;
    using LinearSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>;
    using SchemeType = ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>;
    using StrategyType = ResidualBasedLinearStrategy<
        SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>
    >;

    for (const bool symmetric : {false, true}) {
        Model model;
        ModelPart& r_iga_model_part = CreateCutPatchModelPart(model);

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

        // Unit conductivity, kappa = 0.5, the defaults or the symmetric
        // method with chi = 20, source f = 2, and u_D at the closest points
        auto p_properties = r_iga_model_part.CreateNewProperties(1);
        p_properties->SetValue(CONDUCTIVITY, 1.0);
        p_properties->SetValue(GHOST_PENALTY_FACTOR, 0.5);
        if (symmetric) {
            p_properties->SetValue(ADJOINT_FACTOR, 1.0);
            p_properties->SetValue(PENALTY_FACTOR, 20.0);
        }
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

        // Solution at every Gauss point of the patch, inside and outside Omega
        const double tolerance = 1.0e-10;
        for (const auto& r_element : r_iga_model_part.Elements()) {
            const auto& r_geometry = r_element.GetGeometry();
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();
            double solution = 0.0;
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                solution += r_N(0, i) * r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
            }
            KRATOS_EXPECT_NEAR(
                solution, exact_solution(r_geometry.Center().Coordinates()), tolerance
            );
        }
    }
}

// Transient patch test: the cut patch of GsbmLaplacianCutPatchTest with
// rho c = 1 and u = q + t r, with q the quadratic of that test, r = x - y and
// f = d_t u - Lap(u) = r + 2. After a stationary solve at t = 0 (f = 2), which
// gives q, two steps of BDF1 and one of BDF2 reproduce u exactly, since both
// schemes are exact for solutions linear in time.
KRATOS_TEST_CASE_IN_SUITE(GsbmLaplacianTransientPatchTest, KratosIgaFastSuite)
{
    using SparseSpaceType = TDefaultSparseSpace<double>;
    using LocalSpaceType = TDefaultDenseSpace<double>;
    using LinearSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>;
    using SchemeType = ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>;
    using StrategyType = ResidualBasedLinearStrategy<
        SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>
    >;

    Model model;
    ModelPart& r_iga_model_part = CreateCutPatchModelPart(model);
    r_iga_model_part.SetBufferSize(3);

    const auto exact_solution = [](const array_1d<double, 3>& rPoint, const double Time) {
        const double x = rPoint[0];
        const double y = rPoint[1];
        return 1.0 + x * x - 2.0 * y * y + 3.0 * x * y + Time * (x - y);
    };

    // Unit conductivity and heat capacity, kappa = 0.5 and the default method
    auto p_properties = r_iga_model_part.CreateNewProperties(1);
    p_properties->SetValue(CONDUCTIVITY, 1.0);
    p_properties->SetValue(DENSITY, 1.0);
    p_properties->SetValue(SPECIFIC_HEAT, 1.0);
    p_properties->SetValue(GHOST_PENALTY_FACTOR, 0.5);
    for (auto& r_element : r_iga_model_part.Elements()) {
        r_element.SetProperties(p_properties);
    }
    for (auto& r_condition : r_iga_model_part.Conditions()) {
        r_condition.SetProperties(p_properties);
    }
    for (auto& r_node : r_iga_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }

    StrategyType strategy(
        r_iga_model_part,
        Kratos::make_shared<SchemeType>(),
        Kratos::make_shared<LinearSolverType>()
    );

    // Stationary at t = 0, then BDF1, BDF1 and BDF2 with dt = 0.1
    const double delta_time = 0.1;
    const std::vector<std::vector<double>> bdf_coefficients{
        {}, {1.0 / delta_time, -1.0 / delta_time}, {1.0 / delta_time, -1.0 / delta_time},
        {1.5 / delta_time, -2.0 / delta_time, 0.5 / delta_time}
    };
    for (IndexType step = 0; step < bdf_coefficients.size(); ++step) {
        const double time = step * delta_time;
        if (step > 0) {
            r_iga_model_part.CloneTimeStep(time);
            const auto& r_step_coefficients = bdf_coefficients[step];
            Vector coefficients(r_step_coefficients.size());
            std::copy(r_step_coefficients.begin(), r_step_coefficients.end(), coefficients.begin());
            r_iga_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, coefficients);
        }

        // Source f = 2 at t = 0 and 2 + x - y after, u_D at the closest points
        for (auto& r_element : r_iga_model_part.Elements()) {
            const auto center = r_element.GetGeometry().Center();
            r_element.SetValue(HEAT_FLUX, 2.0 + (step > 0 ? center[0] - center[1] : 0.0));
        }
        for (auto& r_node : model.GetModelPart("Skin.projection_data").Nodes()) {
            r_node.SetValue(TEMPERATURE, exact_solution(r_node.Coordinates(), time));
        }

        strategy.Solve();

        // Solution at every Gauss point of the patch, inside and outside Omega
        const double tolerance = 1.0e-10;
        for (const auto& r_element : r_iga_model_part.Elements()) {
            const auto& r_geometry = r_element.GetGeometry();
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();
            double solution = 0.0;
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                solution += r_N(0, i) * r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
            }
            KRATOS_EXPECT_NEAR(
                solution, exact_solution(r_geometry.Center().Coordinates(), time), tolerance
            );
        }
    }
}

// Enhanced shift patch test: the domain is the ellipse centred at (0.49, 0.51)
// with semi-axes (0.3, 0.2), an exact NURBS curve, minus a triangular polyline
// hole, on the patch [0, 1]^2 with 8 x 8 knot spans of degree p. The
// polynomial u = \sum_{i, j <= p} cos(1 + i + 3j) x^i y^j, with f = -Lap(u),
// lies in the spline space and has total degree 2p: the shift of order p
// misses its mixed terms, while SHIFT_ORDER = 2p, the enhanced shift, is its
// full Taylor series, so u is reproduced exactly only with the latter (up to
// round-off, which grows about 100-fold per degree: p = 5 reaches 1e-6).
KRATOS_TEST_CASE_IN_SUITE(GsbmLaplacianEnhancedShiftPatchTest, KratosIgaFastSuite)
{
    using SparseSpaceType = TDefaultSparseSpace<double>;
    using LocalSpaceType = TDefaultDenseSpace<double>;
    using LinearSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>;
    using SchemeType = ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>;
    using StrategyType = ResidualBasedLinearStrategy<
        SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>
    >;

    const array_1d<double, 3> center{0.49, 0.51, 0.0};
    const std::array<double, 2> semi_axes{0.3, 0.2};

    for (int p = 2; p <= 4; ++p) {
        for (const bool enhanced : {false, true}) {
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

            Parameters nurbs_parameters(R"({
                "model_part_name"     : "IgaModelPart",
                "geometry_name"       : "background",
                "lower_point_xyz"     : [0.0, 0.0, 0.0],
                "upper_point_xyz"     : [1.0, 1.0, 0.0],
                "lower_point_uvw"     : [0.0, 0.0, 0.0],
                "upper_point_uvw"     : [1.0, 1.0, 0.0],
                "polynomial_order"    : [0, 0],
                "number_of_knot_spans": [8, 8]
            })");
            nurbs_parameters["polynomial_order"][0].SetInt(p);
            nurbs_parameters["polynomial_order"][1].SetInt(p);
            NurbsGeometryModeler nurbs_modeler(model, nurbs_parameters);
            nurbs_modeler.SetupGeometryModel();

            // Elements, Dirichlet conditions on the skeleton, and ghost penalty
            // conditions sharing its facets
            Parameters gsbm_parameters(R"({
                "analysis_model_part_name": "IgaModelPart",
                "geometry_name": "background",
                "element_condition_list": [{
                    "type": "element",
                    "name": "GsbmLaplacianElement",
                    "iga_model_part": "Domain",
                    "shape_function_derivatives_order": 2
                }, {
                    "type": "condition",
                    "name": "GsbmLaplacianConditionDirichlet",
                    "iga_model_part": "Dirichlet",
                    "shape_function_derivatives_order": 0,
                    "gsbm_parameters": {}
                }, {
                    "type": "condition",
                    "name": "GsbmLaplacianConditionGhostPenalty",
                    "iga_model_part": "GhostPenalty",
                    "gsbm_parameters": {"skeleton_model_part": "Dirichlet"}
                }]
            })");
            gsbm_parameters["element_condition_list"][1]["shape_function_derivatives_order"]
                .SetInt(p + 1);
            const Modeler& r_registered_modeler = KratosComponents<Modeler>::Get("IgaModelerGsbm");
            r_registered_modeler.Create(model, gsbm_parameters)->SetupModelPart();

            CreateNurbsEllipse(model.CreateModelPart("Skin"), center, semi_axes, "ellipse");
            CreateSkinLoop(
                model.CreateModelPart("Hole"), {{0.40, 0.44}, {0.58, 0.47}, {0.47, 0.60}}
            );
            VolumeFractionGsbmProcess gsbm_process(model, Parameters(R"({
                "model_part_name"                   : "IgaModelPart",
                "geometry_name"                     : "background",
                "elements_model_part_name"          : "IgaModelPart.Domain",
                "skeleton_model_part_name"          : "IgaModelPart.Dirichlet",
                "skin_model_part_outer_initial_name": "Skin",
                "skin_model_part_inner_initial_name": "Hole"
            })"));
            gsbm_process.ExecuteInitialize();

            // \partial_x^A \partial_y^B of u and f = -Lap(u)
            const auto u = [p](const int A, const int B, const array_1d<double, 3>& rPoint) {
                return ComputePolynomialDerivative(p, A, B, rPoint);
            };
            const auto f = [&u](const int A, const int B, const array_1d<double, 3>& rPoint) {
                return -u(A + 2, B, rPoint) - u(A, B + 2, rPoint);
            };

            // Unit conductivity, kappa = 0.5, the default skew-symmetric
            // method without penalty, the shift of order p or 2p, and f at
            // the Gauss points
            auto p_properties = r_iga_model_part.CreateNewProperties(1);
            p_properties->SetValue(CONDUCTIVITY, 1.0);
            p_properties->SetValue(GHOST_PENALTY_FACTOR, 0.5);
            p_properties->SetValue(SHIFT_ORDER, enhanced ? 2 * p : p);
            for (auto& r_element : r_iga_model_part.Elements()) {
                r_element.SetProperties(p_properties);
                const auto& r_center = r_element.GetGeometry().Center().Coordinates();
                r_element.SetValue(HEAT_FLUX, f(0, 0, r_center));
            }
            for (auto& r_condition : r_iga_model_part.Conditions()) {
                r_condition.SetProperties(p_properties);
            }

            // g = u at the closest points of both skins
            for (const std::string name :
                 {"Skin.ellipse.projection_data", "Hole.projection_data"}) {
                ModelPart& r_projection_data = model.GetModelPart(name);
                KRATOS_EXPECT_GT(r_projection_data.NumberOfNodes(), 0);
                for (auto& r_node : r_projection_data.Nodes()) {
                    r_node.SetValue(TEMPERATURE, u(0, 0, r_node.Coordinates()));
                }
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

            // Largest error at the Gauss points of the patch, and of the spans
            // the domain cuts or fills: the void spans, held by the ghost
            // penalty only, amplify round-off
            double max_error = 0.0;
            double max_error_domain = 0.0;
            for (const auto& r_element : r_iga_model_part.Elements()) {
                const auto& r_geometry = r_element.GetGeometry();
                const Matrix& r_N = r_geometry.ShapeFunctionsValues();
                double solution = 0.0;
                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    solution += r_N(0, i) * r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
                }
                const double error =
                    std::abs(solution - u(0, 0, r_geometry.Center().Coordinates()));
                max_error = std::max(max_error, error);
                if (r_element.GetValue(VOLUME_FRACTION) > 0.0) {
                    max_error_domain = std::max(max_error_domain, error);
                }
            }
            if (enhanced) {
                KRATOS_EXPECT_NEAR(max_error_domain, 0.0, 1.0e-7);
                KRATOS_EXPECT_NEAR(max_error, 0.0, 1.0e-4);
            } else {
                KRATOS_EXPECT_GT(max_error_domain, 1.0e-6);
            }
        }
    }
}

} // namespace Kratos::Testing
