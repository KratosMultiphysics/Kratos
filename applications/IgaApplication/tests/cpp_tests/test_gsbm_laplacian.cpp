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
#include <utility>
#include <vector>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_conditions/gsbm_laplacian_condition_dirichlet_stabilized.h"
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
#include "utilities/math_utils.h"

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

/// \partial_x^A \partial_y^B u of u = \sum_{i+j<=p} cos(1 + i + 3j) x^i y^j
double ComputePolynomialDerivative(
    const int Degree,
    const int A,
    const int B,
    const array_1d<double, 3>& rPoint)
{
    double value = 0.0;
    for (int i = A; i <= Degree; ++i) {
        for (int j = B; i + j <= Degree; ++j) {
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

/**
 * @brief (v . \nabla)^NV (w . \nabla)^NW u of the polynomial u of
 * ComputePolynomialDerivative, by the binomial expansion of both powers.
 */
double ContractPolynomialDerivatives(
    const int Degree,
    const array_1d<double, 3>& rV,
    const int NV,
    const array_1d<double, 3>& rW,
    const int NW,
    const array_1d<double, 3>& rPoint)
{
    const auto binomial = [](const int N, const int K) {
        return MathUtils<double>::Factorial(N)
            / (MathUtils<double>::Factorial(K) * MathUtils<double>::Factorial(N - K));
    };
    double value = 0.0;
    for (int a = 0; a <= NV; ++a) {
        for (int c = 0; c <= NW; ++c) {
            // a and c derivatives in x, the rest in y
            const double factor = binomial(NV, a) * std::pow(rV[0], a) * std::pow(rV[1], NV - a)
                * binomial(NW, c) * std::pow(rW[0], c) * std::pow(rW[1], NW - c);
            value += factor * ComputePolynomialDerivative(
                Degree, a + c, NV + NW - a - c, rPoint
            );
        }
    }
    return value;
}

/// The stabilized condition with its shifted shape functions public, to test
/// them
class StabilizedConditionForTesting : public GsbmLaplacianConditionDirichletStabilized
{
public:
    using GsbmLaplacianConditionDirichletStabilized::
        GsbmLaplacianConditionDirichletStabilized;
    using GsbmLaplacianConditionDirichletStabilized::CalculateShiftedShapeFunctions;
};

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

// Stabilized patch test: the domain is the ellipse centred at (0.49, 0.51)
// with semi-axes (0.3, 0.2), an exact NURBS curve, minus a triangular
// polyline hole, whose corners some facet points project onto, on the patch
// [0, 1]^2 with 8 x 8 knot spans of degree p. The polynomial
// u = \sum_{i+j<=p} cos(1 + i + 3j) x^i y^j, with f = -Lap(u), lies in the
// spline space and its series from x_\Gamma is exact, with \rho^m(u) = 0, so
// it is reproduced with and without \rho^m. Since \rho^m(u) = 0 hides D_k,
// each cut side also checks the shifted u_h against, with j = floor(k / 2),
// u(\tilde{x}) - \sum_k (-1)^j t^k / k! \partial_\perp^(2j) \partial_e^(k-2j) u
KRATOS_TEST_CASE_IN_SUITE(GsbmLaplacianStabilizedPatchTest, KratosIgaFastSuite)
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

    // (p, \rho^m added): \rho^m at every degree, the series alone at p = 3
    const std::vector<std::pair<int, bool>> cases{
        {2, true}, {3, true}, {4, true}, {5, true}, {3, false}
    };
    for (const auto& [p, use_residual_terms] : cases) {
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

        // Elements, stabilized Dirichlet conditions on the skeleton, and
        // ghost penalty conditions sharing its facets
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
                "name": "GsbmLaplacianConditionDirichletStabilized",
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
        CreateSkinLoop(model.CreateModelPart("Hole"), {{0.40, 0.44}, {0.58, 0.47}, {0.47, 0.60}});
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

        // Unit conductivity, the paper's chi = 20 and kappa = 0.5, and the
        // source f at the Gauss points
        auto p_properties = r_iga_model_part.CreateNewProperties(1);
        p_properties->SetValue(CONDUCTIVITY, 1.0);
        p_properties->SetValue(PENALTY_FACTOR, 20.0);
        p_properties->SetValue(GHOST_PENALTY_FACTOR, 0.5);
        p_properties->SetValue(USE_RESIDUAL_TERMS, use_residual_terms);
        for (auto& r_element : r_iga_model_part.Elements()) {
            r_element.SetProperties(p_properties);
            r_element.SetValue(HEAT_FLUX, f(0, 0, r_element.GetGeometry().Center().Coordinates()));
        }
        for (auto& r_condition : r_iga_model_part.Conditions()) {
            r_condition.SetProperties(p_properties);
        }

        // At the closest points of both skins: g = u, and f with its
        // derivatives of orders 0 to p - 2
        for (const std::string name : {"Skin.ellipse.projection_data", "Hole.projection_data"}) {
            ModelPart& r_projection_data = model.GetModelPart(name);
            KRATOS_EXPECT_GT(r_projection_data.NumberOfNodes(), 0);
            for (auto& r_node : r_projection_data.Nodes()) {
                const array_1d<double, 3>& r_point = r_node.Coordinates();
                Vector f_derivatives(p * (p - 1) / 2);
                IndexType index = 0;
                for (int n = 0; n <= p - 2; ++n) {
                    for (int k = 0; k <= n; ++k) {
                        f_derivatives[index++] = f(n - k, k, r_point);
                    }
                }
                r_node.SetValue(TEMPERATURE, u(0, 0, r_point));
                r_node.SetValue(VOLUME_SOURCE_DERIVATIVES, f_derivatives);
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

        // Largest error at the Gauss points of the patch, and of the spans the
        // domain cuts or fills: the void spans, held by the ghost penalty only,
        // amplify round-off, to about 4e-6 at p = 5, while the domain's spans
        // stay below 3e-9
        double max_error = 0.0;
        double max_error_domain = 0.0;
        for (const auto& r_element : r_iga_model_part.Elements()) {
            const auto& r_geometry = r_element.GetGeometry();
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();
            double solution = 0.0;
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                solution += r_N(0, i) * r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
            }
            const double error = std::abs(solution - u(0, 0, r_geometry.Center().Coordinates()));
            max_error = std::max(max_error, error);
            if (r_element.GetValue(VOLUME_FRACTION) > 0.0) {
                max_error_domain = std::max(max_error_domain, error);
            }
        }
        KRATOS_EXPECT_NEAR(max_error_domain, 0.0, 1.0e-7);
        KRATOS_EXPECT_NEAR(max_error, 0.0, 1.0e-4);

        if (!use_residual_terms) {
            continue;
        }

        // Shifted u_h on the domain's cut sides, u_h = u
        double max_difference = 0.0;
        for (const auto& r_condition : model.GetModelPart("IgaModelPart.Dirichlet").Conditions()) {
            if (!r_condition.Has(NEIGHBOUR_NODES)) {
                continue;
            }
            const StabilizedConditionForTesting condition(
                0, r_condition.pGetGeometry(), p_properties
            );
            const Node& r_projection_node = r_condition.GetValue(NEIGHBOUR_NODES)[0];
            const auto& r_minus = r_condition.GetGeometry().GetGeometryPart(0);
            const array_1d<double, 3> point = r_minus.Center().Coordinates();
            const array_1d<double, 3> distance = r_projection_node.Coordinates() - point;
            array_1d<double, 3> distance_perp = ZeroVector(3);
            distance_perp[0] = -distance[1];
            distance_perp[1] = distance[0];

            // u(\tilde{x}) - \sum_k (-1)^j / k! (d_\perp . \nabla)^(2j)
            // (-d . \nabla)^(k-2j) u at x_\Gamma, j = floor(k / 2)
            double reference = u(0, 0, point);
            for (int k = 1; k <= p; ++k) {
                const int j = k / 2;
                const double sign = j % 2 == 0 ? 1.0 : -1.0;
                reference -= sign / MathUtils<double>::Factorial(k)
                    * ContractPolynomialDerivatives(
                        p, distance_perp, 2 * j, -distance, k - 2 * j,
                        r_projection_node.Coordinates()
                    );
            }

            for (IndexType side = 0; side < 2; ++side) {
                const auto& r_side = r_condition.GetGeometry().GetGeometryPart(side);
                if (r_side.GetValue(VOLUME_FRACTION) == 0.0) {
                    continue;  // void spans: u_h only up to the ghost penalty's round-off
                }
                const Vector shifted = condition.CalculateShiftedShapeFunctions(
                    r_side, r_projection_node, distance, p
                );
                double shifted_value = 0.0;
                for (IndexType i = 0; i < r_side.size(); ++i) {
                    shifted_value += shifted[i] * r_side[i].FastGetSolutionStepValue(TEMPERATURE);
                }
                max_difference = std::max(max_difference, std::abs(shifted_value - reference));
            }
        }
        KRATOS_EXPECT_NEAR(max_difference, 0.0, 1.0e-7);  // a wrong D_k gives 1e-3
    }
}

} // namespace Kratos::Testing
