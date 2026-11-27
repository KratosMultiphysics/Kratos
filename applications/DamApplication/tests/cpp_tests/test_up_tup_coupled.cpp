// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Permanent Dam U-P / T-U-P coupled-workflow contract. The solid displacement
// domain uses StructuralMechanics small-displacement elements while the Dam
// acoustic elements, coupling conditions and coupled schemes stay unchanged.
// Tests verify the coupling sign and the thermal-solid coupling integration.
//
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Project includes
#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "spaces/ublas_space.h"
#include "utilities/math_utils.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_strategies/schemes/dam_UP_scheme.hpp"
#include "custom_strategies/schemes/dam_P_scheme.hpp"
#include "custom_conditions/UP_condition.hpp"
#include "custom_conditions/added_mass_condition.hpp"
#include "custom_elements/wave_equation_element.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data.
constexpr double uptc_test_young_modulus = 2.0e7;
constexpr double uptc_test_poisson_ratio = 0.2;
constexpr double uptc_test_density = 2400.0;
constexpr double uptc_test_thickness = 0.15;
constexpr double uptc_test_reference_temperature = 20.0;
constexpr double test_beta = 0.25;
constexpr double test_gamma = 0.5;

/// Sparse/dense spaces used by the schemes.
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

/// Max absolute entry difference between two matrices.
double uptc_MaxAbsDiff(const Matrix& rA, const Matrix& rB)
{
    double max_diff = 0.0;
    for (std::size_t i = 0; i < rA.size1(); ++i)
        for (std::size_t j = 0; j < rA.size2(); ++j)
            max_diff = std::max(max_diff, std::abs(rA(i, j) - rB(i, j)));
    return max_diff;
}

/// Max absolute entry difference between two vectors.
double uptc_MaxAbsDiff(const Vector& rA, const Vector& rB)
{
    double max_diff = 0.0;
    for (std::size_t i = 0; i < rA.size(); ++i)
        max_diff = std::max(max_diff, std::abs(rA(i) - rB(i)));
    return max_diff;
}

/// Builds a 2D Q4 U-P model: one solid element (registered name) + one
/// UPCondition on the bottom edge, all nodes carrying U and P DOFs.
ModelPart& CreateUPModel(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rSolidElementName,
    const std::string& rLawName,
    ModelPart*& rOutModelPart)
{
    KRATOS_TRY;
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 2;
    r_pi[SPACE_DIMENSION] = 2;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 0.1;
    r_pi[IS_CONVERGED] = true;

    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY, &ACCELERATION,
                    &VOLUME_ACCELERATION, &PRESSURE, &TEMPERATURE, &NODAL_REFERENCE_TEMPERATURE,
                    &INITIAL_STRESS_TENSOR, &REACTION}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }

    // Dt_PRESSURE / Dt2_PRESSURE are DamApplication variables.
    r_model_part.AddNodalSolutionStepVariable(Dt_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(Dt2_PRESSURE);

    const double coords[4][2] = {{0,0},{2.0,0},{2.0,1.0},{0,1.0}};
    for (std::size_t i = 0; i < 4; ++i) {
        Node::Pointer p_node = r_model_part.CreateNewNode(i + 1, coords[i][0], coords[i][1], 0.0);
        p_node->AddDof(DISPLACEMENT_X, REACTION_X);
        p_node->AddDof(DISPLACEMENT_Y, REACTION_Y);
        p_node->AddDof(PRESSURE);
        p_node->FastGetSolutionStepValue(TEMPERATURE) = uptc_test_reference_temperature;
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = uptc_test_reference_temperature;
        Matrix z3(3, 3);
        noalias(z3) = ZeroMatrix(3, 3);
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
        p_node->FastGetSolutionStepValue(PRESSURE) = 1.0e3;
        p_node->FastGetSolutionStepValue(Dt_PRESSURE) = 10.0;
        p_node->FastGetSolutionStepValue(Dt2_PRESSURE) = 100.0;
        auto& a = p_node->FastGetSolutionStepValue(ACCELERATION);
        a[0] = 0.5; a[1] = 0.5; a[2] = 0.0;
        auto& v = p_node->FastGetSolutionStepValue(VELOCITY);
        v[0] = 0.05; v[1] = -0.05; v[2] = 0.0;
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    (*p_props)[YOUNG_MODULUS] = uptc_test_young_modulus;
    (*p_props)[POISSON_RATIO] = uptc_test_poisson_ratio;
    (*p_props)[DENSITY] = uptc_test_density;
    (*p_props)[THICKNESS] = uptc_test_thickness;
    (*p_props)[THERMAL_EXPANSION] = 1.0e-5;
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        p_props->SetValue(CONSTITUTIVE_LAW, ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrain()));
    } else {
        KRATOS_ERROR << "Unknown law " << rLawName << std::endl;
    }

    r_model_part.CreateNewElement(rSolidElementName, 1, {1, 2, 3, 4}, p_props);
    r_model_part.CreateNewCondition("UPCondition2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);

    r_model_part.pGetElement(1)->Initialize(r_pi);
    r_model_part.pGetCondition(1)->Initialize(r_pi);
    rOutModelPart = &r_model_part;
    return r_model_part;
    KRATOS_CATCH("");
}

/// Builds an isolated 2D U-P condition on a two-node segment. The pressure is
/// uniform, the acceleration is zero and the pressure-acceleration coefficient
/// is set explicitly, so the RHS U-block is the pure pressure coupling force.
ModelPart& CreateUPConditionSegment(
    Model& rModel,
    const std::string& rModelPartName,
    const double X0, const double Y0,
    const double X1, const double Y1,
    ModelPart*& rOutModelPart)
{
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 2;
    r_pi[SPACE_DIMENSION] = 2;
    r_pi[DELTA_TIME] = 0.1;
    r_pi[ACCELERATION_PRESSURE_COEFFICIENT] = 400.0;

    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY,
                    &ACCELERATION, &PRESSURE, &REACTION}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }
    r_model_part.AddNodalSolutionStepVariable(Dt_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(Dt2_PRESSURE);

    const double coords[2][2] = {{X0, Y0}, {X1, Y1}};
    for (std::size_t i = 0; i < 2; ++i) {
        Node::Pointer p_node = r_model_part.CreateNewNode(
            i + 1, coords[i][0], coords[i][1], 0.0);
        p_node->AddDof(DISPLACEMENT_X, REACTION_X);
        p_node->AddDof(DISPLACEMENT_Y, REACTION_Y);
        p_node->AddDof(PRESSURE);
        p_node->FastGetSolutionStepValue(PRESSURE) = 1.0e3;
        noalias(p_node->FastGetSolutionStepValue(ACCELERATION)) = ZeroVector(3);
        noalias(p_node->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewCondition(
        "UPCondition2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.pGetCondition(1)->Initialize(r_pi);
    rOutModelPart = &r_model_part;
    return r_model_part;
}

/// Applies a uniform temperature increment (T-U-P thermal input).
void ApplyTemperature(ModelPart& rModelPart, const double rDeltaTemperature)
{
    for (auto& n : rModelPart.Nodes()) {
        n.FastGetSolutionStepValue(TEMPERATURE) = uptc_test_reference_temperature + rDeltaTemperature;
    }
}

} // namespace

//************************************************************************************
// Coupling sign/energy sanity: positive pressure -> solid force along +normal.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(UPCouplingPreservesExpectedSign, KratosDamFastSuite)
{
    // With a positive uniform pressure, the UPCondition must produce a solid
    // traction in the +normal direction, identical for the historical alias and
    // the direct SMA solid.
    //
    // Bottom edge of the Q4 spans (0,0)-(2,0); the element interior is above.
    // The line tangent is t=(1,0). The Dam 2D normal follows the 3D surface
    // convention n = t x e_z = (t_y, -t_x), so n=(0,-1) for this edge.
    // The condition computes
    //   UPMatrix = -Nu^T (n outer Np),  UVector = -UPMatrix * P
    // so UVector = Nu^T n * P * (edge length). With n=(0,-1), P>0 the U-block
    // traction is (0,-P*L/2) per node. The pressure-acceleration term is
    // PVector = -Density * Np * (n . a); with a=(0.5,0.5) it is +500 per node.
    Model model;
    ModelPart* p_hist_alias_mp = nullptr;
    ModelPart* p_direct_sma_mp = nullptr;
    CreateUPModel(model, "SignHistAlias", "SmallDisplacementSolidElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_hist_alias_mp);
    CreateUPModel(model, "SignDirectSma", "SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_direct_sma_mp);

    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme(test_beta, test_gamma, 0.0, 0.0);
    scheme.Initialize(*p_hist_alias_mp);
    scheme.Initialize(*p_direct_sma_mp);
    CompressedMatrix A;
    Vector Dx, b;
    scheme.InitializeSolutionStep(*p_hist_alias_mp, A, Dx, b);
    scheme.InitializeSolutionStep(*p_direct_sma_mp, A, Dx, b);

    // Isolate the coupling condition RHS (pressure -> U traction).
    Vector rhs_hist_alias, rhs_direct_sma;
    Element::EquationIdVectorType eq_hist_alias, eq_direct_sma;
    Matrix lhs_hist_alias, lhs_direct_sma;
    scheme.CalculateSystemContributions(*p_hist_alias_mp->pGetCondition(1), lhs_hist_alias, rhs_hist_alias, eq_hist_alias, p_hist_alias_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_direct_sma_mp->pGetCondition(1), lhs_direct_sma, rhs_direct_sma, eq_direct_sma, p_direct_sma_mp->GetProcessInfo());

    std::cout << "[coupled] alias condition RHS = [" << rhs_hist_alias[0] << ", " << rhs_hist_alias[1] << ", " << rhs_hist_alias[2]
              << ", " << rhs_hist_alias[3] << ", " << rhs_hist_alias[4] << ", " << rhs_hist_alias[5] << "]" << std::endl;

    // The coupling contributions are identical for the historical alias and the
    // direct SMA solid (the condition never touches the solid element).
    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(rhs_hist_alias, rhs_direct_sma), 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(lhs_hist_alias, lhs_direct_sma), 0.0, 1.0e-12);

    // Condition DOF order per node: [Ux, Uy, P].
    // U-block: positive P -> +P*L/2 along n=(0,-1).
    KRATOS_EXPECT_NEAR(rhs_hist_alias[0], 0.0, 1.0e-9);       // node0 Ux
    KRATOS_EXPECT_NEAR(rhs_hist_alias[1], -1000.0, 1.0e-9);   // node0 Uy
    KRATOS_EXPECT_NEAR(rhs_hist_alias[3], 0.0, 1.0e-9);       // node1 Ux
    KRATOS_EXPECT_NEAR(rhs_hist_alias[4], -1000.0, 1.0e-9);   // node1 Uy
    // P-block: -Density * Np * (n . a) integrated over the edge.
    KRATOS_EXPECT_NEAR(rhs_hist_alias[2], 500.0, 1.0e-9);     // node0 P block
    KRATOS_EXPECT_NEAR(rhs_hist_alias[5], 500.0, 1.0e-9);     // node1 P block

    std::cout << "[coupled] positive P gives traction along the corrected "
              << "normal n=(0,-1); alias and direct SMA are identical."
              << std::endl;
}


//************************************************************************************
// 2D Nu matrix: every entry must be explicitly initialized.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(UPCondition2DNuMatrixIsFullyInitialized, KratosDamFastSuite)
{
    Matrix N = ZeroMatrix(1, 2);
    N(0, 0) = 0.25;
    N(0, 1) = 0.75;

    BoundedMatrix<double, 2, 4> nu;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < nu.size1(); ++i) {
        for (std::size_t j = 0; j < nu.size2(); ++j) {
            nu(i, j) = nan;
        }
    }

    PoroConditionUtilities::CalculateNuMatrix(nu, N, 0);

    KRATOS_EXPECT_NEAR(nu(0, 0), 0.25, 0.0);
    KRATOS_EXPECT_NEAR(nu(0, 1), 0.0, 0.0);
    KRATOS_EXPECT_NEAR(nu(0, 2), 0.75, 0.0);
    KRATOS_EXPECT_NEAR(nu(0, 3), 0.0, 0.0);
    KRATOS_EXPECT_NEAR(nu(1, 0), 0.0, 0.0);
    KRATOS_EXPECT_NEAR(nu(1, 1), 0.25, 0.0);
    KRATOS_EXPECT_NEAR(nu(1, 2), 0.0, 0.0);
    KRATOS_EXPECT_NEAR(nu(1, 3), 0.75, 0.0);
}

//************************************************************************************
// 2D normal: horizontal, vertical and inclined boundaries.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(UPCoupling2DNormalDirection, KratosDamFastSuite)
{
    struct Segment2D
    {
        double x0, y0, x1, y1, nx, ny;
    };

    const double inv_sqrt_two = 1.0 / std::sqrt(2.0);
    const std::vector<Segment2D> segments = {
        {0.0, 0.0, 2.0, 0.0, 0.0, -1.0},                  // horizontal
        {0.0, 0.0, 0.0, 2.0, 1.0, 0.0},                   // vertical
        {0.0, 0.0, 1.0, 1.0, inv_sqrt_two, -inv_sqrt_two} // inclined
    };

    for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& segment = segments[i];
        Model model;
        ModelPart* p_mp = nullptr;
        CreateUPConditionSegment(
            model, "NormalSegment" + std::to_string(i),
            segment.x0, segment.y0, segment.x1, segment.y1, p_mp);

        Matrix lhs;
        Vector rhs;
        p_mp->pGetCondition(1)->CalculateLocalSystem(
            lhs, rhs, p_mp->GetProcessInfo());

        const double length = std::hypot(
            segment.x1 - segment.x0, segment.y1 - segment.y0);
        const double expected_force = 1.0e3 * length / 2.0;
        const double tx = (segment.x1 - segment.x0) / length;
        const double ty = (segment.y1 - segment.y0) / length;

        KRATOS_EXPECT_NEAR(rhs[0], expected_force * segment.nx, 1.0e-9);
        KRATOS_EXPECT_NEAR(rhs[1], expected_force * segment.ny, 1.0e-9);
        KRATOS_EXPECT_NEAR(rhs[3], expected_force * segment.nx, 1.0e-9);
        KRATOS_EXPECT_NEAR(rhs[4], expected_force * segment.ny, 1.0e-9);

        // The coupling force must be perpendicular to the boundary tangent.
        KRATOS_EXPECT_NEAR(rhs[0] * tx + rhs[1] * ty, 0.0, 1.0e-9);

        // With zero acceleration the pressure-acceleration block is zero.
        KRATOS_EXPECT_NEAR(rhs[2], 0.0, 1.0e-12);
        KRATOS_EXPECT_NEAR(rhs[5], 0.0, 1.0e-12);
    }

    // Reversing the line node ordering must reverse the normal and the force.
    Model reversed_model;
    ModelPart* p_reversed_mp = nullptr;
    CreateUPConditionSegment(
        reversed_model, "ReversedNormalSegment", 2.0, 0.0, 0.0, 0.0,
        p_reversed_mp);
    Matrix reversed_lhs;
    Vector reversed_rhs;
    p_reversed_mp->pGetCondition(1)->CalculateLocalSystem(
        reversed_lhs, reversed_rhs, p_reversed_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(reversed_rhs[0], 0.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(reversed_rhs[1], 1000.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(reversed_rhs[3], 0.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(reversed_rhs[4], 1000.0, 1.0e-9);
}

//************************************************************************************
// Added mass 2D: uses the same fully initialized Nu matrix deterministically.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(AddedMassCondition2DUsesFullyInitializedNuMatrix, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("AddedMass2D", 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[ACCELERATION_PRESSURE_COEFFICIENT] = 4.0;

    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &ACCELERATION, &ADDED_MASS}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }

    for (std::size_t i = 0; i < 2; ++i) {
        Node::Pointer p_node = r_model_part.CreateNewNode(
            i + 1, 2.0 * static_cast<double>(i), 0.0, 0.0);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->FastGetSolutionStepValue(ADDED_MASS) = 2.0;
        noalias(p_node->FastGetSolutionStepValue(ACCELERATION)) = ZeroVector(3);
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewCondition(
        "AddedMassCondition2D2N", 1,
        std::vector<ModelPart::IndexType>{1, 2}, p_props);

    Matrix lhs_first, lhs_second;
    Vector rhs_first, rhs_second;
    auto p_condition = r_model_part.pGetCondition(1);
    p_condition->CalculateLocalSystem(lhs_first, rhs_first, r_pi);
    p_condition->CalculateLocalSystem(lhs_second, rhs_second, r_pi);

    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(lhs_first, lhs_second), 0.0, 0.0);
    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(rhs_first, rhs_second), 0.0, 0.0);

    // LHS = coefficient * added_mass * integral(Nu^T Nu) over the edge.
    // Line2D2 uses its default single Gauss point at the centre (N0=N1=0.5)
    // and the edge length is 2, so integral(N0*N0)=integral(N0*N1)=0.5.
    const double diagonal = 4.0 * 2.0 * 0.5;
    const double off_diagonal = 4.0 * 2.0 * 0.5;
    KRATOS_EXPECT_NEAR(lhs_first(0, 0), diagonal, 1.0e-12);
    KRATOS_EXPECT_NEAR(lhs_first(1, 1), diagonal, 1.0e-12);
    KRATOS_EXPECT_NEAR(lhs_first(2, 2), diagonal, 1.0e-12);
    KRATOS_EXPECT_NEAR(lhs_first(3, 3), diagonal, 1.0e-12);
    KRATOS_EXPECT_NEAR(lhs_first(0, 2), off_diagonal, 1.0e-12);
    KRATOS_EXPECT_NEAR(lhs_first(1, 3), off_diagonal, 1.0e-12);

    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(rhs_first, ZeroVector(rhs_first.size())), 0.0, 1.0e-12);
}

//************************************************************************************
// T-U-P: thermal solid + coupling.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(TUPCouplingWithThermalSolid, KratosDamFastSuite)
{
    // T-U-P thermal stage: a temperature increment must generate the same solid
    // thermal stress/force with the historical thermo alias as with the direct
    // SMA solid, and the U-P coupling must remain identical.
    Model model;
    ModelPart* p_hist_alias_mp = nullptr;
    ModelPart* p_direct_sma_mp = nullptr;
    CreateUPModel(model, "TUPhistalias", "SmallDisplacementThermoMechanicElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_hist_alias_mp);
    CreateUPModel(model, "TUPdirectsma", "SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_direct_sma_mp);

    // Reset the pressure/acceleration so the comparison isolates thermal force.
    for (auto* p_mp : {p_hist_alias_mp, p_direct_sma_mp}) {
        for (auto& n : p_mp->Nodes()) {
            n.FastGetSolutionStepValue(PRESSURE) = 0.0;
            auto& a = n.FastGetSolutionStepValue(ACCELERATION);
            a[0] = 0.0; a[1] = 0.0; a[2] = 0.0;
            auto& v = n.FastGetSolutionStepValue(VELOCITY);
            v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
            n.FastGetSolutionStepValue(Dt_PRESSURE) = 0.0;
            n.FastGetSolutionStepValue(Dt2_PRESSURE) = 0.0;
        }
    }
    ApplyTemperature(*p_hist_alias_mp, 25.0);
    ApplyTemperature(*p_direct_sma_mp, 25.0);

    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme(test_beta, test_gamma, 0.0, 0.0);
    scheme.Initialize(*p_hist_alias_mp);
    scheme.Initialize(*p_direct_sma_mp);
    CompressedMatrix A;
    Vector Dx, b;
    scheme.InitializeSolutionStep(*p_hist_alias_mp, A, Dx, b);
    scheme.InitializeSolutionStep(*p_direct_sma_mp, A, Dx, b);

    // Solid thermal internal force via the scheme contributions.
    Matrix lhs_hist_alias, lhs_direct_sma;
    Vector rhs_hist_alias, rhs_direct_sma;
    Element::EquationIdVectorType eq_hist_alias, eq_direct_sma;
    scheme.CalculateSystemContributions(*p_hist_alias_mp->pGetElement(1), lhs_hist_alias, rhs_hist_alias, eq_hist_alias, p_hist_alias_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_direct_sma_mp->pGetElement(1), lhs_direct_sma, rhs_direct_sma, eq_direct_sma, p_direct_sma_mp->GetProcessInfo());

    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(lhs_hist_alias, lhs_direct_sma), 0.0, 1.0e-8);
    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(rhs_hist_alias, rhs_direct_sma), 0.0, 1.0e-8);
    // The thermal expansion must produce a non-trivial internal force.
    KRATOS_EXPECT_GT(uptc_MaxAbsDiff(rhs_hist_alias, ZeroVector(rhs_hist_alias.size())), 1.0);

    // Coupling condition remains identical.
    Vector rhs_c_l, rhs_c_s;
    Element::EquationIdVectorType eq_c_l, eq_c_s;
    Matrix lhs_c_l, lhs_c_s;
    scheme.CalculateSystemContributions(*p_hist_alias_mp->pGetCondition(1), lhs_c_l, rhs_c_l, eq_c_l, p_hist_alias_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_direct_sma_mp->pGetCondition(1), lhs_c_s, rhs_c_s, eq_c_s, p_direct_sma_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(uptc_MaxAbsDiff(rhs_c_l, rhs_c_s), 0.0, 1.0e-12);

    std::cout << "[coupled] T-U-P: thermal force == alias thermo force "
              << "(solid RHS diff=" << uptc_MaxAbsDiff(rhs_hist_alias, rhs_direct_sma)
              << ", coupling diff=" << uptc_MaxAbsDiff(rhs_c_l, rhs_c_s) << ")" << std::endl;
}


//************************************************************************************
// 5. Historical-alias / name independence inside U-P.
//************************************************************************************


} // namespace Testing
} // namespace Kratos
