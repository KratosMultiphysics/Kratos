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
#include <map>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "custom_utilities/shifted_boundary_utilities.h"
#include "geometries/nurbs_surface_geometry.h"
#include "iga_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/constitutive_law.h"
#include "includes/global_pointer_variables.h"
#include "includes/kratos_components.h"
#include "utilities/math_utils.h"

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

/// Field Factor (x - a)^p (y - b)^q on a patch of degrees (p, q), with
/// (x - a)_+^p if TruncatedInX, whose coefficients are the blossoms at the
/// knots t, s of u and v: \prod_m (t_{i+m} - a)_(+) \prod_m (s_{j+m} - b)
void SetProductCoefficients(
    ModelPart& rModelPart,
    const std::string& rPatchName,
    const Variable<double>& rVariable,
    const double Factor,
    const std::array<double, 2>& rShifts,
    const bool TruncatedInX)
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
                const double factor = r_knots_u[i + m] - rShifts[0];
                coefficient *= TruncatedInX ? std::max(factor, 0.0) : factor;
            }
            for (IndexType m = 0; m < r_surface.PolynomialDegreeV(); ++m) {
                coefficient *= r_knots_v[j + m] - rShifts[1];
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

/// Viscosity \mu of the tests on a cut patch
constexpr double StokesCutTestViscosity = 0.7;

/// Linear viscous law \sigma = D \varepsilon, D = \mu diag(2, 2, 1) in Voigt
/// [xx, yy, xy]: Newtonian2DLaw without its volumetric part, which the IGA
/// C++ tests can't load
class GsbmStokesTestViscousLaw : public ConstitutiveLaw
{
public:
    explicit GsbmStokesTestViscousLaw(const double Viscosity)
        : mViscosity(Viscosity)
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GsbmStokesTestViscousLaw>(*this);
    }

    SizeType WorkingSpaceDimension() override
    {
        return 2;
    }

    SizeType GetStrainSize() const override
    {
        return 3;
    }

    void CalculateMaterialResponseCauchy(Parameters& rValues) override
    {
        Matrix& r_D = rValues.GetConstitutiveMatrix();
        r_D = ZeroMatrix(3, 3);
        r_D(0, 0) = 2.0 * mViscosity;
        r_D(1, 1) = 2.0 * mViscosity;
        r_D(2, 2) = mViscosity;
        noalias(rValues.GetStressVector()) = prod(r_D, rValues.GetStrainVector());
    }

private:
    double mViscosity;
};

/// \sigma = 2 \mu \varepsilon(u) - p I of the exact fields of the cut patch
/// tests: u = curl \psi, \psi = (x - 0.4)^2 (y - 0.6)^2, and
/// p = (x - 0.3)(y - 0.7)
BoundedMatrix<double, 2, 2> ComputeStokesCutTestStress(const array_1d<double, 3>& rPoint)
{
    const double x = rPoint[0] - 0.4;
    const double y = rPoint[1] - 0.6;
    const double pressure = (rPoint[0] - 0.3) * (rPoint[1] - 0.7);
    const double mu = StokesCutTestViscosity;

    BoundedMatrix<double, 2, 2> stress;
    stress(0, 0) = 8.0 * mu * x * y - pressure;
    stress(1, 1) = -8.0 * mu * x * y - pressure;
    stress(0, 1) = 2.0 * mu * (x * x - y * y);
    stress(1, 0) = stress(0, 1);
    return stress;
}

/// The divergence-conforming patches (2, 1), (1, 2) and (1, 1) on [0, 1]^2
/// with 6 x 6 knot spans, with GsbmStokesElement, the given condition on the
/// skeleton and GsbmStokesConditionGhostPenalty sharing it, for the fluid
/// inside a 64-gon of radius 0.3 at (0.5, 0.5), so the spans along the patch
/// border are void. The test law, the full Taylor series SHIFT_ORDER = 3, and
/// the exact fields u_x = 2 (x - 0.4)^2 (y - 0.6),
/// u_y = -2 (x - 0.4)(y - 0.6)^2 and p = (x - 0.3)(y - 0.7)
ModelPart& CreateStokesCutModelPart(Model& rModel, const std::string& rConditionName)
{
    ModelPart& r_iga_model_part = rModel.CreateModelPart("IgaModelPart");
    r_iga_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_iga_model_part.AddNodalSolutionStepVariable(PRESSURE);

    const std::array<std::pair<std::string, std::array<int, 2>>, 3> patches{{
        {"VelocityX", {2, 1}}, {"VelocityY", {1, 2}}, {"Pressure", {1, 1}}
    }};
    for (const auto& [r_name, r_degrees] : patches) {
        Parameters parameters(R"({
            "model_part_name"     : "IgaModelPart",
            "geometry_name"       : "",
            "lower_point_xyz"     : [0.0, 0.0, 0.0],
            "upper_point_xyz"     : [1.0, 1.0, 0.0],
            "lower_point_uvw"     : [0.0, 0.0, 0.0],
            "upper_point_uvw"     : [1.0, 1.0, 0.0],
            "polynomial_order"    : [0, 0],
            "number_of_knot_spans": [6, 6]
        })");
        parameters["geometry_name"].SetString(r_name);
        parameters["polynomial_order"][0].SetInt(r_degrees[0]);
        parameters["polynomial_order"][1].SetInt(r_degrees[1]);
        NurbsGeometryModeler(rModel, parameters).SetupGeometryModel();
    }

    const Modeler& r_registered_modeler = KratosComponents<Modeler>::Get("IgaModelerGsbm");
    Parameters modeler_parameters(R"({
        "analysis_model_part_name": "IgaModelPart",
        "geometry_name"           : "Pressure",
        "element_condition_list"  : [{
            "type"                            : "element",
            "name"                            : "GsbmStokesElement",
            "iga_model_part"                  : "Domain",
            "shape_function_derivatives_order": 3,
            "geometry_names"                  : ["VelocityX", "VelocityY", "Pressure"]
        }, {
            "type"                            : "condition",
            "name"                            : "",
            "iga_model_part"                  : "Boundary",
            "shape_function_derivatives_order": 3,
            "gsbm_parameters"                 : {},
            "geometry_names"                  : ["VelocityX", "VelocityY", "Pressure"]
        }, {
            "type"                            : "condition",
            "name"                            : "GsbmStokesConditionGhostPenalty",
            "iga_model_part"                  : "GhostPenalty",
            "gsbm_parameters"                 : {"skeleton_model_part": "Boundary"}
        }]
    })");
    modeler_parameters["element_condition_list"][1]["name"].SetString(rConditionName);
    r_registered_modeler.Create(rModel, modeler_parameters)->SetupModelPart();

    // 64-gon counterclockwise around the fluid, turned off the knot lines
    ModelPart& r_skin = rModel.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(0);
    const IndexType number_of_vertices = 64;
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        const double angle = 2.0 * Globals::Pi * i / number_of_vertices + 0.1;
        r_skin.CreateNewNode(
            i + 1, 0.5 + 0.3 * std::cos(angle), 0.5 + 0.3 * std::sin(angle), 0.0
        );
    }
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        const std::vector<IndexType> node_ids{i + 1, (i + 1) % number_of_vertices + 1};
        r_skin.CreateNewCondition("LineCondition2D2N", i + 1, node_ids, p_skin_properties);
    }
    VolumeFractionGsbmProcess(rModel, Parameters(R"({
        "model_part_name"                   : "IgaModelPart",
        "geometry_name"                     : "Pressure",
        "elements_model_part_name"          : "IgaModelPart.Domain",
        "skeleton_model_part_name"          : "IgaModelPart.Boundary",
        "skin_model_part_outer_initial_name": "Skin"
    })")).ExecuteInitialize();

    for (auto& r_node : r_iga_model_part.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
    }

    auto p_properties = r_iga_model_part.CreateNewProperties(1);
    p_properties->SetValue(
        CONSTITUTIVE_LAW, Kratos::make_shared<GsbmStokesTestViscousLaw>(StokesCutTestViscosity)
    );
    p_properties->SetValue(SHIFT_ORDER, 3);
    p_properties->SetValue(DYNAMIC_VISCOSITY, StokesCutTestViscosity);
    p_properties->SetValue(GHOST_PENALTY_FACTOR, 0.5);
    p_properties->SetValue(PRESSURE_GHOST_PENALTY_FACTOR, 0.1);
    const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
    for (auto& r_element : rModel.GetModelPart("IgaModelPart.Domain").Elements()) {
        r_element.SetProperties(p_properties);
        r_element.Initialize(r_process_info);
    }
    for (auto& r_condition : rModel.GetModelPart("IgaModelPart.Boundary").Conditions()) {
        r_condition.SetProperties(p_properties);
        r_condition.Initialize(r_process_info);
    }

    for (auto& r_condition : rModel.GetModelPart("IgaModelPart.GhostPenalty").Conditions()) {
        r_condition.SetProperties(p_properties);
    }

    SetProductCoefficients(r_iga_model_part, "VelocityX", VELOCITY_X, 2.0, {0.4, 0.6}, false);
    SetProductCoefficients(r_iga_model_part, "VelocityY", VELOCITY_Y, -2.0, {0.4, 0.6}, false);
    SetProductCoefficients(r_iga_model_part, "Pressure", PRESSURE, 1.0, {0.3, 0.7}, false);

    return r_iga_model_part;
}

/// Data t_N = \sigma n + Perturbation of the exact fields at the projection
/// nodes: PRESSURE = p and NORMAL_STRESS = 2 \mu \varepsilon n + Perturbation
void SetNeumannTestData(Model& rModel, const array_1d<double, 3>& rPerturbation)
{
    for (auto& r_node : rModel.GetModelPart("Skin.projection_data").Nodes()) {
        const BoundedMatrix<double, 2, 2> stress = ComputeStokesCutTestStress(r_node.Coordinates());
        const double pressure = (r_node.X() - 0.3) * (r_node.Y() - 0.7);
        const array_1d<double, 3>& r_normal = r_node.GetValue(NORMAL);

        Vector normal_stress(2);
        for (IndexType c = 0; c < 2; ++c) {
            normal_stress[c] = stress(c, 0) * r_normal[0] + stress(c, 1) * r_normal[1]
                + pressure * r_normal[c] + rPerturbation[c];
        }
        r_node.SetValue(PRESSURE, pressure);
        r_node.SetValue(NORMAL_STRESS, normal_stress);
    }
}

/// Velocity u of the exact fields of the cut patch tests
array_1d<double, 3> ComputeStokesCutTestVelocity(const array_1d<double, 3>& rPoint)
{
    const double x = rPoint[0] - 0.4;
    const double y = rPoint[1] - 0.6;

    array_1d<double, 3> velocity = ZeroVector(3);
    velocity[0] = 2.0 * x * x * y;
    velocity[1] = -2.0 * x * y * y;
    return velocity;
}

/// Imposed velocity g = Scale u + Perturbation of the exact fields at the
/// projection nodes
void SetDirichletTestData(
    Model& rModel,
    const array_1d<double, 3>& rPerturbation,
    const double Scale = 1.0)
{
    for (auto& r_node : rModel.GetModelPart("Skin.projection_data").Nodes()) {
        const array_1d<double, 3> velocity =
            Scale * ComputeStokesCutTestVelocity(r_node.Coordinates()) + rPerturbation;
        r_node.SetValue(VELOCITY, velocity);
    }
}

/// Body force MassDensity u + Scale f, with f = -\mu \Delta u + \nabla p of
/// the exact fields, at the Gauss point of each element: f_x =
/// -4 \mu (y - 0.6) + (y - 0.7) and f_y = 4 \mu (x - 0.4) + (x - 0.3)
void SetStokesCutTestBodyForce(
    Model& rModel,
    const double Scale = 1.0,
    const double MassDensity = 0.0)
{
    const double mu = StokesCutTestViscosity;
    for (auto& r_element : rModel.GetModelPart("IgaModelPart.Domain").Elements()) {
        const array_1d<double, 3> point =
            r_element.GetGeometry().GetGeometryPart(0).Center().Coordinates();
        array_1d<double, 3> body_force = MassDensity * ComputeStokesCutTestVelocity(point);
        body_force[0] += Scale * (-4.0 * mu * (point[1] - 0.6) + (point[1] - 0.7));
        body_force[1] += Scale * (4.0 * mu * (point[0] - 0.4) + (point[0] - 0.3));
        r_element.SetValue(BODY_FORCE, body_force);
    }
}

/// Calls a function on the elements, and on the boundary and ghost conditions
template<class TFunction>
void ForEachStokesCutTestEntity(Model& rModel, const TFunction& rFunction)
{
    for (auto& r_element : rModel.GetModelPart("IgaModelPart.Domain").Elements()) {
        rFunction(r_element);
    }
    for (const std::string name : {"IgaModelPart.Boundary", "IgaModelPart.GhostPenalty"}) {
        for (auto& r_condition : rModel.GetModelPart(name).Conditions()) {
            rFunction(r_condition);
        }
    }
}

/// A row per dof of the entities, in the order they are met
std::map<Dof<double>*, IndexType> GetStokesCutTestRows(
    Model& rModel,
    const ProcessInfo& rProcessInfo)
{
    std::map<Dof<double>*, IndexType> rows;
    ForEachStokesCutTestEntity(rModel, [&](auto& rEntity) {
        Element::DofsVectorType dofs;
        rEntity.GetDofList(dofs, rProcessInfo);
        for (Dof<double>* p_dof : dofs) {
            rows.emplace(p_dof, rows.size());
        }
    });
    return rows;
}

/// Pressure dof of the corner (0, 0) of the patch, in the void
Dof<double>* GetStokesCutTestPinnedDof(Model& rModel)
{
    using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;
    auto& r_pressure = dynamic_cast<NurbsSurfaceType&>(
        rModel.GetModelPart("IgaModelPart").GetGeometry("Pressure")
    );
    return r_pressure[0].pGetDof(PRESSURE);
}

/// Solves the dense system of the entities for the increment of the dofs, with
/// the pinned pressure dof set to PinnedValue, as all-Dirichlet leaves the
/// constant pressure free, and updates them
void SolveStokesCutTestIncrement(
    Model& rModel,
    const std::map<Dof<double>*, IndexType>& rRows,
    const ProcessInfo& rProcessInfo,
    const double PinnedValue)
{
    const std::size_t number_of_dofs = rRows.size();
    Matrix system = ZeroMatrix(number_of_dofs, number_of_dofs);
    Vector load = ZeroVector(number_of_dofs);
    ForEachStokesCutTestEntity(rModel, [&](auto& rEntity) {
        Matrix left_hand_side;
        Vector right_hand_side;
        rEntity.CalculateLocalSystem(left_hand_side, right_hand_side, rProcessInfo);
        Element::DofsVectorType dofs;
        rEntity.GetDofList(dofs, rProcessInfo);
        for (IndexType i = 0; i < dofs.size(); ++i) {
            load[rRows.at(dofs[i])] += right_hand_side[i];
            for (IndexType j = 0; j < dofs.size(); ++j) {
                system(rRows.at(dofs[i]), rRows.at(dofs[j])) += left_hand_side(i, j);
            }
        }
    });

    Dof<double>* p_pinned_dof = GetStokesCutTestPinnedDof(rModel);
    const IndexType pinned = rRows.at(p_pinned_dof);
    noalias(row(system, pinned)) = ZeroVector(number_of_dofs);
    system(pinned, pinned) = 1.0;
    load[pinned] = PinnedValue - p_pinned_dof->GetSolutionStepValue();

    Vector increment;
    MathUtils<double>::Solve(system, increment, load);
    for (const auto& [p_dof, row_index] : rRows) {
        p_dof->GetSolutionStepValue() += increment[row_index];
    }
}

/// Adds the right hand side of an entity to the residual of each of its dofs
template<class TEntity>
void AddStokesCutTestResidual(
    std::map<const Dof<double>*, double>& rResidual,
    TEntity& rEntity,
    const ProcessInfo& rProcessInfo)
{
    Matrix left_hand_side;
    Vector right_hand_side;
    rEntity.CalculateLocalSystem(left_hand_side, right_hand_side, rProcessInfo);
    Element::DofsVectorType dofs;
    rEntity.GetDofList(dofs, rProcessInfo);
    for (IndexType i = 0; i < dofs.size(); ++i) {
        rResidual[dofs[i]] += right_hand_side[i];
    }
}

/// Largest absolute value of the residuals
double GetLargestStokesCutTestResidual(const std::map<const Dof<double>*, double>& rResidual)
{
    double largest = 0.0;
    for (const auto& r_entry : rResidual) {
        largest = std::max(largest, std::abs(r_entry.second));
    }
    return largest;
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

        SetProductCoefficients(r_iga_model_part, "VelocityX", VELOCITY_X, 1.0, {0.5, 0.3}, true);
        SetProductCoefficients(r_iga_model_part, "VelocityY", VELOCITY_Y, 2.0, {0.5, 0.3}, true);
        SetProductCoefficients(r_iga_model_part, "Pressure", PRESSURE, 3.0, {0.5, 0.3}, true);
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

// Shifted traction on each cut facet of CreateStokesCutModelPart, with
// the data PRESSURE = p and NORMAL_STRESS = 2 \mu \varepsilon n + \delta at
// the closest point. The full Taylor series is exact for the exact fields, so
// the right hand side of each velocity row is w j N_i (\sigma(\tilde x) n-
// + (n- . n) \delta)_c, with j = \alpha- on the - side and -\alpha+ on the +
// side. The pressure rows have no terms.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionNeumannShiftedTraction, KratosIgaFastSuite)
{
    Model model;
    ModelPart& r_iga_model_part = CreateStokesCutModelPart(model, "GsbmStokesConditionNeumann");
    const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();

    array_1d<double, 3> perturbation = ZeroVector(3);
    perturbation[0] = 0.3;
    perturbation[1] = -0.2;
    SetNeumannTestData(model, perturbation);

    IndexType number_of_cut_facets = 0;
    for (auto& r_facet : model.GetModelPart("IgaModelPart.Boundary").Conditions()) {
        const auto& r_geometry = r_facet.GetGeometry();
        const std::array<double, 2> alpha{
            r_geometry.GetGeometryPart(0).GetValue(VOLUME_FRACTION),
            r_geometry.GetGeometryPart(1).GetValue(VOLUME_FRACTION)
        };
        if (alpha[0] == alpha[1]) {
            continue;
        }
        ++number_of_cut_facets;
        KRATOS_EXPECT_EQ(r_facet.Check(r_process_info), 0);

        Matrix left_hand_side;
        Vector right_hand_side;
        r_facet.CalculateLocalSystem(left_hand_side, right_hand_side, r_process_info);

        // \sigma(\tilde x) n- + (n- . n) \delta
        const auto& r_reference = r_geometry.GetGeometryPart(0).GetGeometryPart(0);
        const double weight =
            r_reference.IntegrationPoints()[0].Weight() * r_reference.DeterminantOfJacobian(0);
        const array_1d<double, 3> normal = r_reference.UnitNormal(0);
        const array_1d<double, 3>& r_skin_normal =
            r_facet.GetValue(NEIGHBOUR_NODES)[0].GetValue(NORMAL);
        const BoundedMatrix<double, 2, 2> stress =
            ComputeStokesCutTestStress(r_reference.Center().Coordinates());
        std::array<double, 2> traction;
        for (IndexType c = 0; c < 2; ++c) {
            traction[c] = stress(c, 0) * normal[0] + stress(c, 1) * normal[1]
                + inner_prod(normal, r_skin_normal) * perturbation[c];
        }

        // u_x, u_y and p of the - side, then of the + side
        Vector expected = ZeroVector(right_hand_side.size());
        IndexType dof = 0;
        for (IndexType side = 0; side < 2; ++side) {
            const double jump_factor = side == 0 ? alpha[0] : -alpha[1];
            for (IndexType part = 0; part < 3; ++part) {
                const auto& r_part = r_geometry.GetGeometryPart(side).GetGeometryPart(part);
                const Matrix& r_N = r_part.ShapeFunctionsValues();
                for (IndexType i = 0; i < r_part.size(); ++i, ++dof) {
                    if (part < 2) {
                        expected[dof] = weight * jump_factor * r_N(0, i) * traction[part];
                    }
                }
            }
        }
        KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, expected, 1.0e-12);
    }
    KRATOS_EXPECT_GT(number_of_cut_facets, 0);
}

// Discrete integration by parts of GSBM on CreateStokesCutModelPart: with
// exact data and the body force f = -\mu \Delta u + \nabla p, the residuals of
// the elements and the Neumann conditions add up to zero in every row
// (momentum, continuity and void). The facet terms are -([[v_\alpha]],
// \sigma n-) where alpha jumps, the quadratures are exact, and the shifted
// traction matches t_N. The elements alone leave a residual.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionNeumannConsistency, KratosIgaFastSuite)
{
    Model model;
    ModelPart& r_iga_model_part = CreateStokesCutModelPart(model, "GsbmStokesConditionNeumann");
    const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
    SetNeumannTestData(model, ZeroVector(3));

    SetStokesCutTestBodyForce(model);

    std::map<const Dof<double>*, double> residual;
    for (auto& r_element : model.GetModelPart("IgaModelPart.Domain").Elements()) {
        AddStokesCutTestResidual(residual, r_element, r_process_info);
    }
    KRATOS_EXPECT_GT(GetLargestStokesCutTestResidual(residual), 1.0e-4);

    for (auto& r_facet : model.GetModelPart("IgaModelPart.Boundary").Conditions()) {
        AddStokesCutTestResidual(residual, r_facet, r_process_info);
    }
    KRATOS_EXPECT_NEAR(GetLargestStokesCutTestResidual(residual), 0.0, 1.0e-12);
}

// Terms of the Dirichlet condition on each cut facet of
// CreateStokesCutModelPart, with g = u + \delta at the closest point. The full
// Taylor series is exact for the exact fields, so S(u) - g = -\delta, and each
// term of the right hand side is known: the consistency w j N_i (\sigma n-)_c,
// with j = \alpha- on the - side and -\alpha+ on the + side, the adjoint
// -\theta w [[\alpha]] V_i . \delta, the q term
// -w [[\alpha]] {{N_p}} n- . \delta on the pressure rows, and the penalty
// pen {{S N_i}} \delta_c. V_i, the viscous traction of the test function, is
// written out by hand.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionDirichletFacetTerms, KratosIgaFastSuite)
{
    struct Variant
    {
        double Theta;
        double Chi;
        bool WithPressure;
    };

    const double mu = StokesCutTestViscosity;
    const double h = 1.0 / 6.0;
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 0.3;
    delta[1] = -0.2;

    const std::array<Variant, 4> variants{{
        {-1.0, 0.0, true}, {1.0, 10.0, true}, {-1.0, 0.0, false}, {1.0, 10.0, false}
    }};
    for (const Variant& r_variant : variants) {
        Model model;
        ModelPart& r_iga_model_part =
            CreateStokesCutModelPart(model, "GsbmStokesConditionDirichlet");
        const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
        Properties& r_properties = r_iga_model_part.GetProperties(1);
        r_properties.SetValue(ADJOINT_FACTOR, r_variant.Theta);
        r_properties.SetValue(PENALTY_FACTOR, r_variant.Chi);
        r_properties.SetValue(ADJOINT_ON_PRESSURE, r_variant.WithPressure);
        SetDirichletTestData(model, delta);

        IndexType number_of_cut_facets = 0;
        for (auto& r_facet : model.GetModelPart("IgaModelPart.Boundary").Conditions()) {
            const auto& r_geometry = r_facet.GetGeometry();
            const std::array<double, 2> alpha{
                r_geometry.GetGeometryPart(0).GetValue(VOLUME_FRACTION),
                r_geometry.GetGeometryPart(1).GetValue(VOLUME_FRACTION)
            };
            if (alpha[0] == alpha[1]) {
                continue;
            }
            ++number_of_cut_facets;
            KRATOS_EXPECT_EQ(r_facet.Check(r_process_info), 0);

            Matrix left_hand_side;
            Vector right_hand_side;
            r_facet.CalculateLocalSystem(left_hand_side, right_hand_side, r_process_info);

            const auto& r_reference = r_geometry.GetGeometryPart(0).GetGeometryPart(0);
            const double weight =
                r_reference.IntegrationPoints()[0].Weight() * r_reference.DeterminantOfJacobian(0);
            const array_1d<double, 3> normal = r_reference.UnitNormal(0);
            const array_1d<double, 3> center = r_reference.Center().Coordinates();
            const array_1d<double, 3> distance =
                r_facet.GetValue(NEIGHBOUR_NODES)[0].Coordinates() - center;
            const BoundedMatrix<double, 2, 2> stress = ComputeStokesCutTestStress(center);

            const double gamma = alpha[1] / (alpha[0] + alpha[1]);
            const double alpha_jump = alpha[0] - alpha[1];
            const double penalty =
                weight * r_variant.Chi * 9.0 * mu / h * std::abs(alpha_jump);

            // u_x, u_y and p of the - side, then of the + side
            Vector expected = ZeroVector(right_hand_side.size());
            IndexType dof = 0;
            for (IndexType side = 0; side < 2; ++side) {
                const auto& r_side = r_geometry.GetGeometryPart(side);
                const double jump_factor = side == 0 ? alpha[0] : -alpha[1];
                const double average_factor = side == 0 ? 1.0 - gamma : gamma;

                for (IndexType part = 0; part < 3; ++part) {
                    const auto& r_part = r_side.GetGeometryPart(part);
                    const Matrix& r_N = r_part.ShapeFunctionsValues();
                    if (part == 2) {
                        for (IndexType i = 0; i < r_part.size(); ++i, ++dof) {
                            if (r_variant.WithPressure) {
                                expected[dof] = -weight * alpha_jump * average_factor * r_N(0, i)
                                    * (normal[0] * delta[0] + normal[1] * delta[1]);
                            }
                        }
                        continue;
                    }

                    const Matrix& r_DN = r_part.ShapeFunctionsLocalGradients()[0];
                    const Vector shifted = ShiftedBoundaryUtilities::ComputeShiftedShapeFunctions(
                        r_part, r_side.GetValue(KNOT_SPAN_INDICES), distance, 3
                    );
                    for (IndexType i = 0; i < r_part.size(); ++i, ++dof) {
                        // (V_i)_m = \mu (\delta_{mc} grad N . n- + n-_c d_m N)
                        const double normal_derivative =
                            r_DN(i, 0) * normal[0] + r_DN(i, 1) * normal[1];
                        double viscous_along_delta = 0.0;
                        for (IndexType m = 0; m < 2; ++m) {
                            viscous_along_delta += average_factor * mu * delta[m]
                                * ((m == part ? normal_derivative : 0.0)
                                    + normal[part] * r_DN(i, m));
                        }

                        expected[dof] = weight * jump_factor * r_N(0, i)
                                * (stress(part, 0) * normal[0] + stress(part, 1) * normal[1])
                            - r_variant.Theta * weight * alpha_jump * viscous_along_delta
                            + penalty * average_factor * shifted[i] * delta[part];
                    }
                }
            }
            KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, expected, 1.0e-10);
        }
        KRATOS_EXPECT_GT(number_of_cut_facets, 0);
    }
}

// Discrete integration by parts of GSBM with the Dirichlet condition: with
// exact data and the body force of the exact fields, the residuals of the
// elements and the conditions add up to zero in every row, with and without
// the q term, as S(u) = g leaves only the consistency term.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionDirichletConsistency, KratosIgaFastSuite)
{
    for (const bool with_pressure : {true, false}) {
        Model model;
        ModelPart& r_iga_model_part =
            CreateStokesCutModelPart(model, "GsbmStokesConditionDirichlet");
        const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
        Properties& r_properties = r_iga_model_part.GetProperties(1);
        r_properties.SetValue(ADJOINT_FACTOR, 1.0);
        r_properties.SetValue(PENALTY_FACTOR, 10.0);
        r_properties.SetValue(ADJOINT_ON_PRESSURE, with_pressure);
        SetDirichletTestData(model, ZeroVector(3));
        SetStokesCutTestBodyForce(model);

        std::map<const Dof<double>*, double> residual;
        for (auto& r_element : model.GetModelPart("IgaModelPart.Domain").Elements()) {
            AddStokesCutTestResidual(residual, r_element, r_process_info);
        }
        KRATOS_EXPECT_GT(GetLargestStokesCutTestResidual(residual), 1.0e-4);

        for (auto& r_facet : model.GetModelPart("IgaModelPart.Boundary").Conditions()) {
            AddStokesCutTestResidual(residual, r_facet, r_process_info);
        }
        KRATOS_EXPECT_NEAR(GetLargestStokesCutTestResidual(residual), 0.0, 1.0e-11);
    }
}

// First solve of the GSBM Stokes problem: element, ghost penalty and Dirichlet
// condition from zero nodal values, assembled into a dense system, recover the
// exact fields in every coefficient, void included, with the q term and the
// standard ghost and without it and the ghost on the divergence. All-Dirichlet
// leaves the constant pressure free, so one pressure coefficient in the void
// is pinned to its exact value.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesConditionDirichletSolve, KratosIgaFastSuite)
{
    for (const bool with_pressure : {true, false}) {
        Model model;
        ModelPart& r_iga_model_part =
            CreateStokesCutModelPart(model, "GsbmStokesConditionDirichlet");
        const ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
        Properties& r_properties = r_iga_model_part.GetProperties(1);
        r_properties.SetValue(ADJOINT_ON_PRESSURE, with_pressure);
        r_properties.SetValue(PRESSURE_GHOST_PENALTY_ON_DIVERGENCE, !with_pressure);
        SetDirichletTestData(model, ZeroVector(3));
        SetStokesCutTestBodyForce(model);

        // The exact coefficients, then zero nodal values
        const auto rows = GetStokesCutTestRows(model, r_process_info);
        Vector exact(rows.size());
        for (const auto& [p_dof, row_index] : rows) {
            exact[row_index] = p_dof->GetSolutionStepValue();
            p_dof->GetSolutionStepValue() = 0.0;
        }

        const double pinned_value = exact[rows.at(GetStokesCutTestPinnedDof(model))];
        SolveStokesCutTestIncrement(model, rows, r_process_info, pinned_value);

        double largest_error = 0.0;
        for (const auto& [p_dof, row_index] : rows) {
            largest_error = std::max(
                largest_error, std::abs(p_dof->GetSolutionStepValue() - exact[row_index])
            );
        }
        KRATOS_EXPECT_NEAR(largest_error, 0.0, 1.0e-8);
    }
}

// Transient solve of the GSBM Stokes problem on CreateStokesCutModelPart: the
// fields (1 + t) u and (1 + t) p, linear in time, with density 2 and the body
// force \rho u + (1 + t) f. From their values at t = 0, a step of BDF1 and two
// of BDF2, exact for fields linear in time, recover every coefficient, void
// included, also with PSPG, whose residual has to include the time term.
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesElementTransientSolve, KratosIgaFastSuite)
{
    const double delta_time = 0.1;
    const double density = 2.0;

    for (const double pspg_factor : {0.0, 0.4}) {
        Model model;
        ModelPart& r_iga_model_part =
            CreateStokesCutModelPart(model, "GsbmStokesConditionDirichlet");
        r_iga_model_part.SetBufferSize(3);
        ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
        r_iga_model_part.GetProperties(1).SetValue(DENSITY, density);
        r_iga_model_part.GetProperties(1).SetValue(PSPG_FACTOR, pspg_factor);

        const auto rows = GetStokesCutTestRows(model, r_process_info);
        Vector initial(rows.size());
        for (const auto& [p_dof, row_index] : rows) {
            initial[row_index] = p_dof->GetSolutionStepValue();
        }

        const std::vector<std::vector<double>> bdf_coefficients{
            {1.0 / delta_time, -1.0 / delta_time},
            {1.5 / delta_time, -2.0 / delta_time, 0.5 / delta_time},
            {1.5 / delta_time, -2.0 / delta_time, 0.5 / delta_time}
        };
        for (IndexType step = 0; step < bdf_coefficients.size(); ++step) {
            const double time = (step + 1) * delta_time;
            const double scale = 1.0 + time;
            r_iga_model_part.CloneTimeStep(time);

            const auto& r_step_coefficients = bdf_coefficients[step];
            Vector coefficients(r_step_coefficients.size());
            std::copy(
                r_step_coefficients.begin(), r_step_coefficients.end(), coefficients.begin()
            );
            r_process_info.SetValue(BDF_COEFFICIENTS, coefficients);

            SetDirichletTestData(model, ZeroVector(3), scale);
            SetStokesCutTestBodyForce(model, scale, density);

            const double pinned_value =
                scale * initial[rows.at(GetStokesCutTestPinnedDof(model))];
            SolveStokesCutTestIncrement(model, rows, r_process_info, pinned_value);

            double largest_error = 0.0;
            for (const auto& [p_dof, row_index] : rows) {
                largest_error = std::max(
                    largest_error,
                    std::abs(p_dof->GetSolutionStepValue() - scale * initial[row_index])
                );
            }
            KRATOS_EXPECT_NEAR(largest_error, 0.0, 1.0e-8);
        }
    }
}

// With BDF_COEFFICIENTS, Check needs DENSITY and a buffer as long as the
// coefficients
KRATOS_TEST_CASE_IN_SUITE(GsbmStokesElementTransientCheck, KratosIgaFastSuite)
{
    Model model;
    ModelPart& r_iga_model_part = CreateStokesCutModelPart(model, "GsbmStokesConditionDirichlet");
    r_iga_model_part.SetBufferSize(3);
    ProcessInfo& r_process_info = r_iga_model_part.GetProcessInfo();
    Properties& r_properties = r_iga_model_part.GetProperties(1);
    auto& r_element = model.GetModelPart("IgaModelPart.Domain").Elements().front();

    r_process_info.SetValue(BDF_COEFFICIENTS, Vector(3, 1.0));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_element.Check(r_process_info), "needs DENSITY");

    r_properties.SetValue(DENSITY, 1.0);
    KRATOS_EXPECT_EQ(r_element.Check(r_process_info), 0);

    r_process_info.SetValue(BDF_COEFFICIENTS, Vector(4, 1.0));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_element.Check(r_process_info), "buffer size is smaller");
}

} // namespace Kratos::Testing
