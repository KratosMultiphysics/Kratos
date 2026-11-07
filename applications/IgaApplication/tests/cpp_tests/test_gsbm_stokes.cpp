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

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "geometries/nurbs_surface_geometry.h"
#include "iga_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/kratos_components.h"

namespace Kratos::Testing
{

namespace
{

using IndexType = std::size_t;

/// Patches VelocityX, VelocityY and Pressure of the given degrees on
/// [0, 1] x [0, 0.75], with 4 x 3 knot spans of size 0.25, and
/// GsbmStokesConditionGhostPenalty on the facets of the interior knot lines
ModelPart& CreateStokesGhostPenaltyModelPart(
    Model& rModel,
    const std::array<std::array<int, 2>, 3>& rDegrees)
{
    ModelPart& r_iga_model_part = rModel.CreateModelPart("IgaModelPart");
    r_iga_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_iga_model_part.AddNodalSolutionStepVariable(PRESSURE);

    const std::array<std::string, 3> names{"VelocityX", "VelocityY", "Pressure"};
    for (IndexType part = 0; part < 3; ++part) {
        Parameters parameters(R"({
            "model_part_name"     : "IgaModelPart",
            "geometry_name"       : "",
            "lower_point_xyz"     : [0.0, 0.0, 0.0],
            "upper_point_xyz"     : [1.0, 0.75, 0.0],
            "lower_point_uvw"     : [0.0, 0.0, 0.0],
            "upper_point_uvw"     : [1.0, 0.75, 0.0],
            "polynomial_order"    : [0, 0],
            "number_of_knot_spans": [4, 3]
        })");
        parameters["geometry_name"].SetString(names[part]);
        parameters["polynomial_order"][0].SetInt(rDegrees[part][0]);
        parameters["polynomial_order"][1].SetInt(rDegrees[part][1]);
        NurbsGeometryModeler(rModel, parameters).SetupGeometryModel();
    }

    const Modeler& r_registered_modeler = KratosComponents<Modeler>::Get("IgaModelerGsbm");
    r_registered_modeler.Create(rModel, Parameters(R"({
        "analysis_model_part_name": "IgaModelPart",
        "geometry_name"           : "Pressure",
        "element_condition_list"  : [{
            "type"                            : "condition",
            "name"                            : "GsbmStokesConditionGhostPenalty",
            "iga_model_part"                  : "GhostPenalty",
            "shape_function_derivatives_order": 4,
            "gsbm_parameters"                 : {},
            "geometry_names"                  : ["VelocityX", "VelocityY", "Pressure"]
        }]
    })"))->SetupModelPart();

    for (auto& r_node : r_iga_model_part.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
    }

    return r_iga_model_part;
}

/// Sets the volume fractions of the - and + sides of every facet
void SetFacetVolumeFractions(
    ModelPart& rFacets,
    const double AlphaMinus,
    const double AlphaPlus)
{
    for (auto& r_facet : rFacets.Conditions()) {
        r_facet.GetGeometry().GetGeometryPart(0).SetValue(VOLUME_FRACTION, AlphaMinus);
        r_facet.GetGeometry().GetGeometryPart(1).SetValue(VOLUME_FRACTION, AlphaPlus);
    }
}

/// Field Factor (x - 0.5)_+^a (y - 0.3)^b on a patch of degrees (a, b), whose
/// coefficients are the blossoms at the knots t, s of u and v:
/// \prod_m (t_{i+m} - 0.5)_+ \prod_m (s_{j+m} - 0.3)
void SetTruncatedPowerCoefficients(
    ModelPart& rModelPart,
    const std::string& rPatchName,
    const Variable<double>& rVariable,
    const double Factor)
{
    using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;
    auto& r_surface = dynamic_cast<NurbsSurfaceType&>(rModelPart.GetGeometry(rPatchName));
    const Vector& r_knots_u = r_surface.KnotsU();
    const Vector& r_knots_v = r_surface.KnotsV();
    const IndexType number_of_points_u = r_surface.NumberOfControlPointsU();

    for (IndexType j = 0; j < r_surface.NumberOfControlPointsV(); ++j) {
        for (IndexType i = 0; i < number_of_points_u; ++i) {
            double coefficient = Factor;
            for (IndexType m = 0; m < r_surface.PolynomialDegreeU(); ++m) {
                coefficient *= std::max(r_knots_u[i + m] - 0.5, 0.0);
            }
            for (IndexType m = 0; m < r_surface.PolynomialDegreeV(); ++m) {
                coefficient *= r_knots_v[j + m] - 0.3;
            }
            r_surface[j * number_of_points_u + i].FastGetSolutionStepValue(rVariable) =
                coefficient;
        }
    }
}

/// I_n = \int_0^{0.75} (y - 0.3)^n dy
double IntegrateShiftedPower(const int Exponent)
{
    return (std::pow(0.45, Exponent + 1) - std::pow(-0.3, Exponent + 1)) / (Exponent + 1);
}

} // namespace

// Ghost penalty on the divergence-conforming patches (3, 2), (2, 3) and (2, 2)
// with alpha = 0.5 on every facet, for u_x = (x - 0.5)_+^3 (y - 0.3)^2,
// u_y = 2 (x - 0.5)_+^2 (y - 0.3)^3 and p = 3 (x - 0.5)_+^2 (y - 0.3)^2. They
// are polynomials across every knot line but x = 0.5, where D_n^k jumps by k!
// times the factor in y, with k the degree in x. The facet rule integrates the
// squared jumps exactly, so the energies summed over the facets are those of
// the formulas: the velocity ghost, and either the pressure ghost s(q, p) or
// the mixed one +s(q, div u), -s(div v, p). With alpha = 1 on both sides, or a
// full jump of alpha, the facets are not penalized.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionGhostPenaltyJumps, KratosIgaFastSuite)
{
    const double kappa_u = 0.7;
    const double kappa_p = 0.3;
    const double viscosity = 0.02;
    const double grad_div_factor = 0.5;
    const double h = 0.25;
    const double i_4 = IntegrateShiftedPower(4);
    const double i_6 = IntegrateShiftedPower(6);
    const double tolerance = 1.0e-10;

    for (const bool on_divergence : {false, true}) {
        Model model;
        ModelPart& r_iga_model_part = CreateStokesGhostPenaltyModelPart(
            model, {{{3, 2}, {2, 3}, {2, 2}}}
        );
        ModelPart& r_facets = model.GetModelPart("IgaModelPart.GhostPenalty");
        const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();

        auto p_properties = r_iga_model_part.CreateNewProperties(1);
        p_properties->SetValue(GHOST_PENALTY_FACTOR, kappa_u);
        p_properties->SetValue(PRESSURE_GHOST_PENALTY_FACTOR, kappa_p);
        p_properties->SetValue(DYNAMIC_VISCOSITY, viscosity);
        p_properties->SetValue(GRAD_DIV_FACTOR, grad_div_factor);
        p_properties->SetValue(PRESSURE_GHOST_PENALTY_ON_DIVERGENCE, on_divergence);
        for (auto& r_facet : r_facets.Conditions()) {
            r_facet.SetProperties(p_properties);
        }

        SetTruncatedPowerCoefficients(r_iga_model_part, "VelocityX", VELOCITY_X, 1.0);
        SetTruncatedPowerCoefficients(r_iga_model_part, "VelocityY", VELOCITY_Y, 2.0);
        SetTruncatedPowerCoefficients(r_iga_model_part, "Pressure", PRESSURE, 3.0);
        SetFacetVolumeFractions(r_facets, 0.5, 0.5);

        // Energies of the velocity (u) and pressure (p) parts of the dof values
        // with the local matrices K, summed over the facets
        double energy_uu = 0.0;
        double energy_pp = 0.0;
        double energy_pu = 0.0;
        double energy_up = 0.0;
        for (auto& r_facet : r_facets.Conditions()) {
            KRATOS_EXPECT_EQ(r_facet.Check(r_process_info), 0);

            Matrix left_hand_side;
            Vector right_hand_side;
            r_facet.CalculateLocalSystem(left_hand_side, right_hand_side, r_process_info);
            Condition::DofsVectorType dofs;
            r_facet.GetDofList(dofs, r_process_info);
            KRATOS_EXPECT_EQ(dofs.size(), left_hand_side.size1());

            Vector velocity = ZeroVector(dofs.size());
            Vector pressure = ZeroVector(dofs.size());
            for (IndexType i = 0; i < dofs.size(); ++i) {
                Vector& r_part = dofs[i]->GetVariable() == PRESSURE ? pressure : velocity;
                r_part[i] = dofs[i]->GetSolutionStepValue();
            }

            // Residual form, in the order of the dofs: -K (u + p)
            const Vector residual = -prod(left_hand_side, velocity + pressure);
            KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, residual, 1.0e-12);

            energy_uu += inner_prod(velocity, prod(left_hand_side, velocity));
            energy_pp += inner_prod(pressure, prod(left_hand_side, pressure));
            energy_pu += inner_prod(pressure, prod(left_hand_side, velocity));
            energy_up += inner_prod(velocity, prod(left_hand_side, pressure));
        }

        // Jumps 6 (y - 0.3)^2 of D_n^3 u_x and 4 (y - 0.3)^3 of D_n^2 u_y
        KRATOS_EXPECT_RELATIVE_NEAR(
            energy_uu,
            kappa_u * viscosity * (std::pow(h, 5) * i_4 / 7.0 + 4.0 * std::pow(h, 3) * i_6 / 5.0),
            tolerance
        );

        if (!on_divergence) {
            // Jump 6 (y - 0.3)^2 of D_n^2 p, in s(p, p) / (\mu + \beta_div)
            KRATOS_EXPECT_RELATIVE_NEAR(
                energy_pp,
                9.0 * kappa_p * std::pow(h, 5) * i_4 / (5.0 * (viscosity + grad_div_factor)),
                tolerance
            );
            KRATOS_EXPECT_NEAR(energy_pu, 0.0, 1.0e-14);
            KRATOS_EXPECT_NEAR(energy_up, 0.0, 1.0e-14);
        } else {
            // Jump 18 (y - 0.3)^2 of D_n^2 div u, in +s(p, div u), -s(div u, p)
            const double mixed_energy = 5.4 * kappa_p * std::pow(h, 5) * i_4;
            KRATOS_EXPECT_NEAR(energy_pp, 0.0, 1.0e-14);
            KRATOS_EXPECT_RELATIVE_NEAR(energy_pu, mixed_energy, tolerance);
            KRATOS_EXPECT_RELATIVE_NEAR(energy_up, -mixed_energy, tolerance);
        }

        // Not penalized: alpha = 1 on both sides, or a full jump of alpha
        const std::array<std::pair<double, double>, 2> unpenalized{{{1.0, 1.0}, {0.0, 1.0}}};
        for (const auto& [alpha_minus, alpha_plus] : unpenalized) {
            SetFacetVolumeFractions(r_facets, alpha_minus, alpha_plus);
            for (auto& r_facet : r_facets.Conditions()) {
                Matrix left_hand_side;
                Vector right_hand_side;
                r_facet.CalculateLocalSystem(left_hand_side, right_hand_side, r_process_info);
                KRATOS_EXPECT_DOUBLE_EQ(norm_frobenius(left_hand_side), 0.0);
            }
        }
    }
}

// The mixed pressure ghost penalty needs divergence-conforming spaces, which
// the Taylor-Hood patches (2, 2), (2, 2) and (1, 1) are not, while the
// standard one takes any spaces
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionGhostPenaltyCheck, KratosIgaFastSuite)
{
    Model model;
    ModelPart& r_iga_model_part = CreateStokesGhostPenaltyModelPart(
        model, {{{2, 2}, {2, 2}, {1, 1}}}
    );
    const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();

    auto p_properties = r_iga_model_part.CreateNewProperties(1);
    p_properties->SetValue(GHOST_PENALTY_FACTOR, 1.0);
    p_properties->SetValue(PRESSURE_GHOST_PENALTY_FACTOR, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0);
    auto& r_facet = model.GetModelPart("IgaModelPart.GhostPenalty").Conditions().front();
    r_facet.SetProperties(p_properties);

    KRATOS_EXPECT_EQ(r_facet.Check(r_process_info), 0);
    p_properties->SetValue(PRESSURE_GHOST_PENALTY_ON_DIVERGENCE, true);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        r_facet.Check(r_process_info), "needs divergence-conforming spaces"
    );
}

} // namespace Kratos::Testing
