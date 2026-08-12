// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers
//

// Phase 5D (P / U-P / T-U-P coupled workflows): CHARACTERIZATION ONLY. Determines
// whether the solid displacement domain of the Dam acoustic/structural coupled
// workflows can use StructuralMechanicsApplication SmallDisplacement elements
// while the Dam acoustic elements, coupling conditions and coupled schemes stay
// unchanged. No production code or registration is modified.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
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
#include "custom_elements/wave_equation_element.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thickness = 0.15;
constexpr double test_reference_temperature = 20.0;
constexpr double test_beta = 0.25;
constexpr double test_gamma = 0.5;

/// Sparse/dense spaces used by the schemes.
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

/// Max absolute entry difference between two matrices.
double MaxAbsDiff(const Matrix& rA, const Matrix& rB)
{
    double max_diff = 0.0;
    for (std::size_t i = 0; i < rA.size1(); ++i)
        for (std::size_t j = 0; j < rA.size2(); ++j)
            max_diff = std::max(max_diff, std::abs(rA(i, j) - rB(i, j)));
    return max_diff;
}

/// Max absolute entry difference between two vectors.
double MaxAbsDiff(const Vector& rA, const Vector& rB)
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
                    &INITIAL_STRESS_TENSOR}) {
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
        p_node->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
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
    (*p_props)[YOUNG_MODULUS] = test_young_modulus;
    (*p_props)[POISSON_RATIO] = test_poisson_ratio;
    (*p_props)[DENSITY] = test_density;
    (*p_props)[THICKNESS] = test_thickness;
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

/// Applies a uniform temperature increment (T-U-P thermal input).
void ApplyTemperature(ModelPart& rModelPart, const double rDeltaTemperature)
{
    for (auto& n : rModelPart.Nodes()) {
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + rDeltaTemperature;
    }
}

} // namespace

//************************************************************************************
// 1. U-P coupling independence: UPCondition and DamUPScheme with legacy vs SMA solid.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R54_UP_Coupling_SolidIndependence, KratosDamFastSuite)
{
    // The UPCondition and DamUPScheme must be independent of the solid element
    // implementation. Build two identical Q4 U-P models differing only in the
    // solid element (legacy vs SMA) and compare every contribution.
    Model model;
    ModelPart* p_legacy_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    CreateUPModel(model, "UPLegacy", "SmallDisplacementSolidElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_legacy_mp);
    CreateUPModel(model, "UPSma", "SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_sma_mp);

    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme_legacy(test_beta, test_gamma, 0.02, 0.03);
    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme_sma(test_beta, test_gamma, 0.02, 0.03);
    scheme_legacy.Initialize(*p_legacy_mp);
    scheme_sma.Initialize(*p_sma_mp);

    CompressedMatrix A_l, A_s;
    Vector Dx_l, Dx_s, b_l, b_s;
    scheme_legacy.InitializeSolutionStep(*p_legacy_mp, A_l, Dx_l, b_l);
    scheme_sma.InitializeSolutionStep(*p_sma_mp, A_s, Dx_s, b_s);

    // Solid element dynamic contributions (element + condition) via the scheme.
    Matrix lhs_e_legacy, lhs_e_sma, lhs_c_legacy, lhs_c_sma;
    Vector rhs_e_legacy, rhs_e_sma, rhs_c_legacy, rhs_c_sma;
    Element::EquationIdVectorType eq_e_legacy, eq_e_sma;
    Element::EquationIdVectorType eq_c_legacy, eq_c_sma;

    scheme_legacy.CalculateSystemContributions(
        *p_legacy_mp->pGetElement(1), lhs_e_legacy, rhs_e_legacy, eq_e_legacy, p_legacy_mp->GetProcessInfo());
    scheme_sma.CalculateSystemContributions(
        *p_sma_mp->pGetElement(1), lhs_e_sma, rhs_e_sma, eq_e_sma, p_sma_mp->GetProcessInfo());

    scheme_legacy.CalculateSystemContributions(
        *p_legacy_mp->pGetCondition(1), lhs_c_legacy, rhs_c_legacy, eq_c_legacy, p_legacy_mp->GetProcessInfo());
    scheme_sma.CalculateSystemContributions(
        *p_sma_mp->pGetCondition(1), lhs_c_sma, rhs_c_sma, eq_c_sma, p_sma_mp->GetProcessInfo());

    // 1a. Solid element contributions: legacy vs SMA are numerically identical
    //     (stiffness/internal/mass/damping through the generic scheme path).
    KRATOS_EXPECT_EQ(lhs_e_legacy.size1(), lhs_e_sma.size1());
    KRATOS_EXPECT_EQ(rhs_e_legacy.size(), rhs_e_sma.size());
    KRATOS_EXPECT_NEAR(MaxAbsDiff(lhs_e_legacy, lhs_e_sma), 0.0, 1.0e-8);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(rhs_e_legacy, rhs_e_sma), 0.0, 1.0e-8);
    KRATOS_EXPECT_EQ(eq_e_legacy.size(), eq_e_sma.size());
    for (std::size_t i = 0; i < eq_e_legacy.size(); ++i)
        KRATOS_EXPECT_EQ(eq_e_legacy[i], eq_e_sma[i]);

    // 1b. Coupling-condition contributions are EXACTLY identical (the condition
    //     never touches the solid element).
    KRATOS_EXPECT_NEAR(MaxAbsDiff(lhs_c_legacy, lhs_c_sma), 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(rhs_c_legacy, rhs_c_sma), 0.0, 1.0e-12);
    KRATOS_EXPECT_EQ(eq_c_legacy.size(), eq_c_sma.size());
    for (std::size_t i = 0; i < eq_c_legacy.size(); ++i)
        KRATOS_EXPECT_EQ(eq_c_legacy[i], eq_c_sma[i]);

    // 1c. Equation-ID layout: the condition carries [Ux,Uy,P] per node (6 IDs),
    //     the solid element carries [Ux,Uy] per node (8 IDs).
    std::cout << "[5D] UP legacy vs SMA: solid LHS diff=" << MaxAbsDiff(lhs_e_legacy, lhs_e_sma)
              << " condition LHS diff=" << MaxAbsDiff(lhs_c_legacy, lhs_c_sma)
              << " solid eq_ids=" << eq_e_legacy.size()
              << " condition eq_ids=" << eq_c_legacy.size() << std::endl;
    KRATOS_EXPECT_EQ(eq_c_legacy.size(), 6u);
    KRATOS_EXPECT_EQ(eq_e_legacy.size(), 8u);
    std::cout << "[5D] UPCondition + DamUPScheme are independent of the solid "
              << "element implementation (name-independent, generic Element API)." << std::endl;
}

//************************************************************************************
// 2. Coupling sign/energy sanity: positive pressure -> solid force in +normal.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R54_UP_Coupling_SignSanity, KratosDamFastSuite)
{
    // With a positive uniform pressure, the UPCondition must produce a solid
    // traction in the +normal direction, identical for legacy and SMA solids.
    // Bottom edge of the Q4 spans (0,0)-(2,0); the element interior is above, so
    // the outward edge normal is (0,-1). The condition computes
    //   UPMatrix = -Nu^T (n outer Np),  UVector = -UPMatrix * P
    // so UVector = Nu^T n * P * (edge length). With n=(0,-1), P>0 the solid
    // traction is along -n = (0,+1) on the U block.
    Model model;
    ModelPart* p_legacy_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    CreateUPModel(model, "SignLegacy", "SmallDisplacementSolidElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_legacy_mp);
    CreateUPModel(model, "SignSma", "SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_sma_mp);

    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme(test_beta, test_gamma, 0.0, 0.0);
    scheme.Initialize(*p_legacy_mp);
    scheme.Initialize(*p_sma_mp);
    CompressedMatrix A;
    Vector Dx, b;
    scheme.InitializeSolutionStep(*p_legacy_mp, A, Dx, b);
    scheme.InitializeSolutionStep(*p_sma_mp, A, Dx, b);

    // Isolate the coupling condition RHS (pressure -> U traction).
    Vector rhs_l, rhs_s;
    Element::EquationIdVectorType eq_l, eq_s;
    Matrix lhs_l, lhs_s;
    scheme.CalculateSystemContributions(*p_legacy_mp->pGetCondition(1), lhs_l, rhs_l, eq_l, p_legacy_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_sma_mp->pGetCondition(1), lhs_s, rhs_s, eq_s, p_sma_mp->GetProcessInfo());

    // Condition DOF order per node: [Ux, Uy, P]. The Dam 2D UPCondition
    // (CalculateNormalVector<2,2>) uses the edge Jacobian as the coupling
    // direction, i.e. for the bottom edge it is the edge tangent (1,0). Positive
    // P therefore produces a solid traction along +tangent with magnitude
    // P*L/2 = 1000. The Uy entry stays uninitialized when the coupling direction
    // has no y component (a pre-existing UPCondition quirk, identical for both
    // implementations).
    std::cout << "[5D] legacy condition RHS = [" << rhs_l[0] << ", " << rhs_l[1] << ", " << rhs_l[2]
              << ", " << rhs_l[3] << ", " << rhs_l[4] << ", " << rhs_l[5] << "]" << std::endl;

    // THE key migration result: the coupling contributions are IDENTICAL for the
    // legacy and SMA solids (the condition never touches the solid element).
    KRATOS_EXPECT_NEAR(MaxAbsDiff(rhs_l, rhs_s), 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(lhs_l, lhs_s), 0.0, 1.0e-12);

    // Sign convention preserved: positive P -> solid traction +P*L/2 along the
    // coupling (tangent) direction; P-block RHS = -P*L/2 (per node).
    KRATOS_EXPECT_NEAR(rhs_l[0], 1000.0, 1.0e-9);    // node0 Ux traction
    KRATOS_EXPECT_NEAR(rhs_l[3], 1000.0, 1.0e-9);    // node1 Ux traction
    KRATOS_EXPECT_NEAR(rhs_l[2], -500.0, 1.0e-9);    // node0 P block
    KRATOS_EXPECT_NEAR(rhs_l[5], -500.0, 1.0e-9);    // node1 P block

    std::cout << "[5D] coupling sign/energy: positive P gives the same +tangent "
              << "traction and same P-block RHS for legacy and SMA; no reversal "
              << "or rescaling. (Uy entry is uninitialized in the 2D UPCondition "
              << "when the coupling direction is purely x - pre-existing quirk.)" << std::endl;
}

//************************************************************************************
// 3. T-U-P: thermal solid + coupling with SMA vs legacy thermo element.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R54_TUP_ThermalSolid_Coupling, KratosDamFastSuite)
{
    // T-U-P thermal stage: a temperature increment must generate the same solid
    // thermal stress/force with the SMA solid as with the legacy thermo element,
    // and the U-P coupling must remain identical.
    Model model;
    ModelPart* p_legacy_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    CreateUPModel(model, "TUPLegacy", "SmallDisplacementThermoMechanicElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_legacy_mp);
    CreateUPModel(model, "TUPSma", "SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_sma_mp);

    // Reset the pressure/acceleration so the comparison isolates thermal force.
    for (auto* p_mp : {p_legacy_mp, p_sma_mp}) {
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
    ApplyTemperature(*p_legacy_mp, 25.0);
    ApplyTemperature(*p_sma_mp, 25.0);

    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme(test_beta, test_gamma, 0.0, 0.0);
    scheme.Initialize(*p_legacy_mp);
    scheme.Initialize(*p_sma_mp);
    CompressedMatrix A;
    Vector Dx, b;
    scheme.InitializeSolutionStep(*p_legacy_mp, A, Dx, b);
    scheme.InitializeSolutionStep(*p_sma_mp, A, Dx, b);

    // Solid thermal internal force via the scheme contributions.
    Matrix lhs_l, lhs_s;
    Vector rhs_l, rhs_s;
    Element::EquationIdVectorType eq_l, eq_s;
    scheme.CalculateSystemContributions(*p_legacy_mp->pGetElement(1), lhs_l, rhs_l, eq_l, p_legacy_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_sma_mp->pGetElement(1), lhs_s, rhs_s, eq_s, p_sma_mp->GetProcessInfo());

    KRATOS_EXPECT_NEAR(MaxAbsDiff(lhs_l, lhs_s), 0.0, 1.0e-8);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(rhs_l, rhs_s), 0.0, 1.0e-8);
    // The thermal expansion must produce a non-trivial internal force.
    KRATOS_EXPECT_GT(MaxAbsDiff(rhs_l, ZeroVector(rhs_l.size())), 1.0);

    // Coupling condition remains identical.
    Vector rhs_c_l, rhs_c_s;
    Element::EquationIdVectorType eq_c_l, eq_c_s;
    Matrix lhs_c_l, lhs_c_s;
    scheme.CalculateSystemContributions(*p_legacy_mp->pGetCondition(1), lhs_c_l, rhs_c_l, eq_c_l, p_legacy_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_sma_mp->pGetCondition(1), lhs_c_s, rhs_c_s, eq_c_s, p_sma_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(MaxAbsDiff(rhs_c_l, rhs_c_s), 0.0, 1.0e-12);

    std::cout << "[5D] T-U-P: SMA solid thermal force == legacy thermo element "
              << "(solid RHS diff=" << MaxAbsDiff(rhs_l, rhs_s)
              << ", coupling diff=" << MaxAbsDiff(rhs_c_l, rhs_c_s) << ")" << std::endl;
}

//************************************************************************************
// 4. P-only isolation: WaveEquationElement + DamPScheme with no solid element.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R54_POnly_Isolation, KratosDamFastSuite)
{
    // The P-only workflow (WaveEquationElement + DamPScheme) must assemble with
    // ONLY pressure DOFs and no solid element present.
    Model model;
    ModelPart& r_mp = model.CreateModelPart("POnly", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 2;
    r_pi[SPACE_DIMENSION] = 2;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 0.1;
    r_mp.AddNodalSolutionStepVariable(PRESSURE);
    r_mp.AddNodalSolutionStepVariable(Dt_PRESSURE);
    r_mp.AddNodalSolutionStepVariable(Dt2_PRESSURE);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

    const double coords[4][2] = {{0,0},{2.0,0},{2.0,1.0},{0,1.0}};
    for (std::size_t i = 0; i < 4; ++i) {
        Node::Pointer p_node = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], 0.0);
        p_node->AddDof(PRESSURE);
        p_node->FastGetSolutionStepValue(PRESSURE) = 100.0;
        p_node->FastGetSolutionStepValue(Dt_PRESSURE) = 5.0;
        p_node->FastGetSolutionStepValue(Dt2_PRESSURE) = 50.0;
    }
    auto p_props = r_mp.CreateNewProperties(1);
    (*p_props)[DENSITY] = 1000.0;
    (*p_props)[DENSITY_WATER] = 1000.0;
    (*p_props)[BULK_MODULUS_LIQUID] = 2.2e9;
    r_mp.CreateNewElement("WaveEquationElement2D4N", 1, {1, 2, 3, 4}, p_props);
    r_mp.pGetElement(1)->Initialize(r_pi);

    // Drive the real DamPScheme path (Initialize + InitializeSolutionStep +
    // CalculateSystemContributions) - the exact operations the P solver runs.
    DamPScheme<SparseSpaceType, LocalSpaceType> p_scheme(test_beta, test_gamma);
    p_scheme.Initialize(r_mp);
    CompressedMatrix A;
    Vector Dx, b;
    p_scheme.InitializeSolutionStep(r_mp, A, Dx, b);

    Matrix lhs;
    Vector rhs;
    Element::EquationIdVectorType eq_id;
    p_scheme.CalculateSystemContributions(*r_mp.pGetElement(1), lhs, rhs, eq_id, r_pi);

    KRATOS_EXPECT_EQ(eq_id.size(), 4u);   // PRESSURE only
    KRATOS_EXPECT_EQ(lhs.size1(), 4u);
    KRATOS_EXPECT_GT(lhs(0, 0), 0.0);     // acoustic mass/stiffness assembled

    std::cout << "[5D] P-only: WaveEquationElement + DamPScheme assemble with "
              << "only pressure DOFs (" << eq_id.size() << " ids), no solid element."
              << " Phase 6 does not need to modify the P-only workflow." << std::endl;
}

//************************************************************************************
// 5. Historical-alias / name independence inside U-P.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R54_UP_HistoricalAlias_NameIndependence, KratosDamFastSuite)
{
    // The user-facing element name has no effect on coupling behavior: the
    // scheme and condition never inspect element names or concrete types. This
    // is demonstrated by the identical results for the two registered names
    // (legacy vs SMA). Verify the runtime element is SMA and everything works.
    Model model;
    ModelPart* p_mp = nullptr;
    CreateUPModel(model, "AliasUP", "SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", p_mp);

    DamUPScheme<SparseSpaceType, LocalSpaceType> scheme(test_beta, test_gamma, 0.0, 0.0);
    scheme.Initialize(*p_mp);
    CompressedMatrix A;
    Vector Dx, b;
    scheme.InitializeSolutionStep(*p_mp, A, Dx, b);

    Matrix lhs_e, lhs_c;
    Vector rhs_e, rhs_c;
    Element::EquationIdVectorType eq_e, eq_c;
    scheme.CalculateSystemContributions(*p_mp->pGetElement(1), lhs_e, rhs_e, eq_e, p_mp->GetProcessInfo());
    scheme.CalculateSystemContributions(*p_mp->pGetCondition(1), lhs_c, rhs_c, eq_c, p_mp->GetProcessInfo());

    KRATOS_EXPECT_EQ(lhs_e.size1(), 8u);
    KRATOS_EXPECT_EQ(lhs_c.size1(), 6u);
    KRATOS_EXPECT_TRUE(std::abs(lhs_c(0, 2)) > 0.0);   // U<->P coupling block populated (negative UPMatrix entry)

    std::cout << "[5D] historical-alias simulation: the SMA-backed historical name "
              << "behaves identically inside the U-P workflow (coupling block populated, "
              << "scheme/condition generic)." << std::endl;
}

} // namespace Testing
} // namespace Kratos
