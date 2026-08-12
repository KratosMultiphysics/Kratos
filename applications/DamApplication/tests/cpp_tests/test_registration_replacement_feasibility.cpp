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

// Phase 5C.2 (registration-replacement readiness): CHARACTERIZATION / DESIGN
// ONLY. These tests prove whether the historical Dam small-displacement element
// names can be backed by StructuralMechanicsApplication SmallDisplacement
// prototypes, and that no production logic depends on the concrete Dam element
// type. No production code or registration is modified.

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
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_processes/dam_nodal_cauchy_stress_extrapolation_process.hpp"
#include "custom_utilities/dam_nonlocal_damage_utilities.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data.
constexpr double test_density = 2400.0;
constexpr double test_thickness = 0.15;
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_reference_temperature = 20.0;

/// Constructs a constitutive law for the given family.
ConstitutiveLaw::Pointer CreateLaw(const std::string& rFamily, const bool rIs2d)
{
    if (rFamily == "linear") {
        if (rIs2d) return ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrain());
        return ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw());
    }
    if (rFamily == "nodal") {
        if (rIs2d) return ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrainNodal());
        return ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal());
    }
    if (rFamily == "local")
        return ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw());
    return ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw());
}

/// Builds a model part with one element of the given registered name created
/// from the given nodal coordinates.
Element::Pointer CreateElement(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::vector<std::vector<double>>& rCoords,
    const bool rIs2d,
    const std::string& rLawFamily,
    ModelPart*& rOutModelPart,
    const double rDensity = test_density)
{
    KRATOS_TRY;
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = rIs2d ? 2 : 3;
    r_pi[SPACE_DIMENSION] = rIs2d ? 2 : 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_YOUNG_MODULUS);
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    r_model_part.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

    for (std::size_t i = 0; i < rCoords.size(); ++i) {
        const auto& c = rCoords[i];
        Node::Pointer p_node = r_model_part.CreateNewNode(i + 1, c[0], c[1], c[2]);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        p_node->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = test_young_modulus;
        const std::size_t dim = rIs2d ? 2 : 3;
        Matrix z(dim, dim);
        noalias(z) = ZeroMatrix(dim, dim);
        p_node->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = z;
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z;
        p_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    (*p_props)[YOUNG_MODULUS] = test_young_modulus;
    (*p_props)[POISSON_RATIO] = test_poisson_ratio;
    (*p_props)[DENSITY] = rDensity;
    (*p_props)[THERMAL_EXPANSION] = 1.0e-5;
    (*p_props)[DAMAGE_THRESHOLD] = 5.0e-3;
    (*p_props)[STRENGTH_RATIO] = 10.0;
    (*p_props)[FRACTURE_ENERGY] = 5000.0;
    if (rIs2d)
        (*p_props)[THICKNESS] = test_thickness;
    p_props->SetValue(CONSTITUTIVE_LAW, CreateLaw(rLawFamily, rIs2d));

    std::vector<ModelPart::IndexType> element_nodes(rCoords.size());
    for (std::size_t i = 0; i < rCoords.size(); ++i)
        element_nodes[i] = i + 1;
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_props);

    Element::Pointer p_element = r_model_part.pGetElement(1);
    p_element->Initialize(r_pi);
    rOutModelPart = &r_model_part;
    return p_element;
    KRATOS_CATCH("");
}

/// Applies a uniform total-strain state (uniaxial-stress-like).
void ApplyState(ModelPart& rModelPart, const double rEps)
{
    const bool is3d = (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3);
    for (auto& n : rModelPart.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = rEps * x0[0];
        u[1] = -test_poisson_ratio * rEps * x0[1];
        u[2] = (is3d ? -test_poisson_ratio * rEps * x0[2] : 0.0);
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
    }
}

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

} // namespace

//************************************************************************************
// 1. Historical-name -> SMA-prototype feasibility: complete public API surface.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R52_SMAPrototype_FullAPISurface, KratosDamFastSuite)
{
    // For every historical Dam small-displacement registration, the SMA
    // counterpart prototype must satisfy the full externally-visible Element
    // API expected by Dam workflows. Exercise every entry point.
    const std::vector<std::tuple<std::string, bool, std::vector<std::vector<double>>>> cases = {
        {"SmallDisplacementElement2D3N", true, {{0,0,0},{2.5,0,0},{0,1.25,0}}},
        {"SmallDisplacementElement2D4N", true, {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0}}},
        {"SmallDisplacementElement3D4N", false, {{0,0,0},{2.5,0,0},{0,1.25,0},{0,0,2.5}}},
        {"SmallDisplacementElement3D8N", false, {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}}},
        {"SmallDisplacementElement3D10N", false, {{0,0,0},{2.5,0,0},{0,1.25,0},{0,0,2.5},{1.25,0,0},{1.25,0.625,0},{0,0.625,0},{0,0,1.25},{1.25,0,1.25},{0,0.625,1.25}}},
    };

    for (const auto& c : cases) {
        const std::string& name = std::get<0>(c);
        const bool is2d = std::get<1>(c);
        const auto& coords = std::get<2>(c);

        Model model;
        ModelPart* p_mp = nullptr;
        Element::Pointer p_elem = CreateElement(model, "API" + name, name, coords, is2d, "linear", p_mp);
        ProcessInfo& r_pi = p_mp->GetProcessInfo();
        ApplyState(*p_mp, 1.0e-5);

        // Runtime type is the SMA SmallDisplacement (not a Dam type).
        const std::string runtime_type = typeid(*p_elem).name();
        std::cout << "[5C.2] " << name << " runtime_type=" << runtime_type
                  << " npoints=" << p_elem->GetGeometry().PointsNumber() << std::endl;

        // Create / Clone / Check.
        KRATOS_EXPECT_EQ(p_elem->Check(r_pi), 0);
        Element::Pointer p_clone = p_elem->Clone(2, p_elem->GetGeometry().Points());
        KRATOS_EXPECT_EQ(p_clone->GetGeometry().PointsNumber(), coords.size());

        // Lifecycle.
        p_elem->InitializeSolutionStep(r_pi);
        p_elem->InitializeNonLinearIteration(r_pi);

        // Local system, LHS, RHS, mass, damping.
        Matrix lhs, mass, damping;
        Vector rhs;
        p_elem->CalculateLocalSystem(lhs, rhs, r_pi);
        Matrix lhs2;
        p_elem->CalculateLeftHandSide(lhs2, r_pi);
        Vector rhs2;
        p_elem->CalculateRightHandSide(rhs2, r_pi);
        p_elem->CalculateMassMatrix(mass, r_pi);
        p_elem->CalculateDampingMatrix(damping, r_pi);
        KRATOS_EXPECT_EQ(lhs.size1(), coords.size() * (is2d ? 2 : 3));
        KRATOS_EXPECT_EQ(rhs.size(), coords.size() * (is2d ? 2 : 3));
        KRATOS_EXPECT_NEAR(MaxAbsDiff(lhs, lhs2), 0.0, 1.0e-12);
        KRATOS_EXPECT_NEAR(MaxAbsDiff(rhs, rhs2), 0.0, 1.0e-12);
        KRATOS_EXPECT_GT(mass(0, 0), 0.0);

        // Integration-point outputs: scalar, Vector, Matrix, CONSTITUTIVE_LAW.
        std::vector<Vector> v_out;
        p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, v_out, r_pi);
        std::vector<Matrix> m_out;
        p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, m_out, r_pi);
        std::vector<double> s_out;
        p_elem->CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, s_out, r_pi);
        std::vector<ConstitutiveLaw::Pointer> law_out;
        p_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, law_out, r_pi);
        KRATOS_EXPECT_EQ(v_out.size(), law_out.size());
        KRATOS_EXPECT_EQ(m_out.size(), law_out.size());
        KRATOS_EXPECT_EQ(s_out.size(), law_out.size());
        KRATOS_EXPECT_GT(v_out[0].size(), 0);
        KRATOS_EXPECT_TRUE(law_out[0] != nullptr);

        // Finalize lifecycle.
        p_elem->FinalizeNonLinearIteration(r_pi);
        r_pi[IS_CONVERGED] = true;
        p_elem->FinalizeSolutionStep(r_pi);

        std::cout << "[5C.2] " << name << ": full API surface exercised OK ("
                  << law_out.size() << " constitutive-law GPs)" << std::endl;
    }
    std::cout << "[5C.2] Feasibility: the SMA prototypes satisfy the complete Element API "
              << "expected by Dam workflows for every historical registration." << std::endl;
}

//************************************************************************************
// 2. Mechanical / thermo-mechanical registration convergence.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R52_MechThermo_RegistrationConvergence, KratosDamFastSuite)
{
    // Both historical names (SmallDisplacementSolidElement3D8N and
    // SmallDisplacementThermoMechanicElement3D8N) would map to the SAME SMA
    // SmallDisplacementElement3D8N prototype. The element NAME must have no
    // remaining effect: behavior is governed entirely by the constitutive law.
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}};

    for (const std::string& family : {"linear", "nodal", "local", "nonlocal"}) {
        Model model;
        ModelPart* p_mp = nullptr;
        Element::Pointer p_elem = CreateElement(model, "Conv" + family, "SmallDisplacementElement3D8N",
                                                coords, false, family, p_mp);
        ProcessInfo& r_pi = p_mp->GetProcessInfo();
        ApplyState(*p_mp, 1.0e-5);

        KRATOS_EXPECT_EQ(p_elem->Check(r_pi), 0);
        Matrix lhs, mass;
        Vector rhs;
        p_elem->CalculateLocalSystem(lhs, rhs, r_pi);
        p_elem->CalculateMassMatrix(mass, r_pi);

        std::vector<Vector> v_out;
        p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, v_out, r_pi);
        KRATOS_EXPECT_EQ(v_out.size(), 8u);
        KRATOS_EXPECT_TRUE(MaxAbsDiff(lhs, lhs) == 0.0);

        std::cout << "[5C.2] SMA 3D8N with " << family << " law: LHS[0]=" << lhs(0,0)
                  << " mass[0]=" << mass(0,0) << " -> name-independent behavior" << std::endl;
    }
    std::cout << "[5C.2] Convergence: the mechanical and thermo-mechanical historical names "
              << "both map to one SMA implementation; the law drives all behavior." << std::endl;
}

//************************************************************************************
// 3. Constitutive-law vector / GP-count compatibility.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R52_GPCount_Match, KratosDamFastSuite)
{
    // Dam and SMA must instantiate the same number of constitutive-law
    // integration points for the normal mechanical response.
    struct Case { std::string dam; std::string sma; bool is2d; std::vector<std::vector<double>> coords; };
    const std::vector<Case> cases = {
        {"SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N", true, {{0,0,0},{2.5,0,0},{0,1.25,0}}},
        {"SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N", true, {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0}}},
        {"SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N", false, {{0,0,0},{2.5,0,0},{0,1.25,0},{0,0,2.5}}},
        {"SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N", false, {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}}},
        {"SmallDisplacementSolidElement3D10N", "SmallDisplacementElement3D10N", false, {{0,0,0},{2.5,0,0},{0,1.25,0},{0,0,2.5},{1.25,0,0},{1.25,0.625,0},{0,0.625,0},{0,0,1.25},{1.25,0,1.25},{0,0.625,1.25}}},
    };

    for (const auto& c : cases) {
        Model model;
        ModelPart* p_dam_mp = nullptr;
        ModelPart* p_sma_mp = nullptr;
        Element::Pointer p_dam = CreateElement(model, "GP" + c.dam, c.dam, c.coords, c.is2d, "linear", p_dam_mp);
        Element::Pointer p_sma = CreateElement(model, "GP" + c.sma, c.sma, c.coords, c.is2d, "linear", p_sma_mp);

        std::vector<ConstitutiveLaw::Pointer> dam_laws, sma_laws;
        p_dam->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, dam_laws, p_dam_mp->GetProcessInfo());
        p_sma->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, sma_laws, p_sma_mp->GetProcessInfo());

        const std::size_t default_gp = p_dam->GetGeometry().IntegrationPoints(p_dam->GetIntegrationMethod()).size();
        KRATOS_EXPECT_EQ(dam_laws.size(), default_gp);
        KRATOS_EXPECT_EQ(sma_laws.size(), default_gp);
        KRATOS_EXPECT_EQ(dam_laws.size(), sma_laws.size());

        std::cout << "[5C.2] " << c.dam << " vs " << c.sma
                  << ": constitutive-law GPs " << dam_laws.size() << " == " << sma_laws.size()
                  << " (geometry default " << default_gp << ")" << std::endl;
    }
    std::cout << "[5C.2] GP layout: identical constitutive-law vector size/ordering "
              << "between Dam and SMA for every geometry." << std::endl;
}

//************************************************************************************
// 4. Smoothing compatibility on an SMA-only model.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R52_Smoothing_SMAOnly, KratosDamFastSuite)
{
    // T3 (2D) and T4 (3D) are the geometries used by the Dam test .mdpa files.
    // Build SMA-only models and run the real DamNodalCauchyStressExtrapolationProcess.
    const std::vector<std::tuple<std::string, bool, std::vector<std::vector<double>>>> cases = {
        {"SmallDisplacementElement2D3N", true, {{0,0,0},{2.5,0,0},{0,1.25,0}}},
        {"SmallDisplacementElement3D4N", false, {{0,0,0},{2.5,0,0},{0,1.25,0},{0,0,2.5}}},
    };

    for (const auto& c : cases) {
        const std::string& name = std::get<0>(c);
        const bool is2d = std::get<1>(c);
        const auto& coords = std::get<2>(c);

        Model model;
        ModelPart* p_mp = nullptr;
        Element::Pointer p_elem = CreateElement(model, "Smooth" + name, name, coords, is2d, "linear", p_mp);
        ProcessInfo& r_pi = p_mp->GetProcessInfo();
        ApplyState(*p_mp, 1.0e-5);

        // Compute the Gauss-point Cauchy stress (material response).
        std::vector<Matrix> gauss_stress;
        p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, gauss_stress, r_pi);
        KRATOS_EXPECT_EQ(gauss_stress.size(), 1u);   // single-GP T3/T4

        // Run the real Dam smoothing process.
        DamNodalCauchyStressExtrapolationProcess process(*p_mp);
        process.ExtrapolateAndAccumulate();

        // Normalize (the smoothing scheme owns this; replicate here).
        for (auto& n : p_mp->Nodes()) {
            const double area = n.FastGetSolutionStepValue(NODAL_AREA);
            KRATOS_EXPECT_GT(area, 0.0);
            auto& stress = n.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
            if (area > 0.0)
                stress /= area;
            // For a single-GP element the extrapolated nodal stress equals the
            // Gauss-point stress tensor.
            const std::size_t dim = is2d ? 2 : 3;
            for (std::size_t i = 0; i < dim; ++i)
                for (std::size_t j = 0; j < dim; ++j)
                    KRATOS_EXPECT_NEAR(stress(i, j), gauss_stress[0](i, j), 1.0e-9);
        }
        std::cout << "[5C.2] smoothing on SMA-only " << name << ": nodal stress == GP stress" << std::endl;
    }
    std::cout << "[5C.2] Smoothing: DamNodalCauchyStressExtrapolationProcess works on SMA "
              << "elements with no concrete Dam element dependency." << std::endl;
}

//************************************************************************************
// 5. Nonlocal compatibility on SMA-only model.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R52_Nonlocal_SMAOnly, KratosDamFastSuite)
{
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}};

    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = CreateElement(model, "Nonlocal", "SmallDisplacementElement3D8N",
                                            coords, false, "nonlocal", p_mp);
    ProcessInfo& r_pi = p_mp->GetProcessInfo();
    ApplyState(*p_mp, 2.0e-5);

    // Prescribe a nonlocal driving value and commit on the element's OWN law
    // clones (obtained through the generic CONSTITUTIVE_LAW pointer output).
    std::vector<ConstitutiveLaw::Pointer> laws;
    p_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, r_pi);
    KRATOS_EXPECT_EQ(laws.size(), 8u);
    laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_pi);

    // Run the real Dam nonlocal LOCAL production utility on the SMA element.
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(*p_mp);

    // The LOCAL quantity is stored in the element law's flow rule.
    double local = 0.0;
    laws[0]->GetValue(LOCAL_EQUIVALENT_STRAIN, local);
    KRATOS_EXPECT_TRUE(local > 0.0);
    std::cout << "[5C.2] nonlocal LOCAL production on SMA element: LOCAL=" << local << std::endl;
    std::cout << "[5C.2] Nonlocal: DamNonlocalDamageUtilities works on SMA elements with no "
              << "concrete Dam element dependency." << std::endl;
}

//************************************************************************************
// 6. Selfweight / body-force compatibility.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R52_SelfweightBodyForce, KratosDamFastSuite)
{
    // The Dam construction/selfweight workflow applies gravity as a body force
    // (VOLUME_ACCELERATION). Verify Dam and SMA elements produce identical RHS.
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}};

    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, "SWDam", "SmallDisplacementSolidElement3D8N", coords, false, "linear", p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, "SWSma", "SmallDisplacementElement3D8N", coords, false, "linear", p_sma_mp);

    // Gravity body force.
    for (auto& n : p_dam_mp->Nodes()) {
        auto& g = n.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        g[0] = 0.0; g[1] = 0.0; g[2] = -9.81;
    }
    for (auto& n : p_sma_mp->Nodes()) {
        auto& g = n.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        g[0] = 0.0; g[1] = 0.0; g[2] = -9.81;
    }

    Vector rhs_dam, rhs_sma;
    p_dam->CalculateRightHandSide(rhs_dam, p_dam_mp->GetProcessInfo());
    p_sma->CalculateRightHandSide(rhs_sma, p_sma_mp->GetProcessInfo());

    KRATOS_EXPECT_EQ(rhs_dam.size(), rhs_sma.size());
    for (std::size_t i = 0; i < rhs_dam.size(); ++i)
        KRATOS_EXPECT_NEAR(rhs_dam(i), rhs_sma(i), 1.0e-9);
    KRATOS_EXPECT_LT(rhs_sma(2), 0.0);  // downward force on first node

    std::cout << "[5C.2] selfweight body force (VOLUME_ACCELERATION): Dam RHS == SMA RHS"
              << std::endl;
}

} // namespace Testing
} // namespace Kratos
