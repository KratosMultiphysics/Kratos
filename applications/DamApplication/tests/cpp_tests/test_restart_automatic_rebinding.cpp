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

// Phase 6C.1: automatic restart continuation of the Dam thermal damage
// constitutive-law families. After a binary restart the flow/yield/hardening
// transient Properties bindings are reconstructed automatically from the current
// ConstitutiveLaw::Parameters; no manual SetProperties/ReinitializeMaterialProperties
// call is required in any solver, element, scheme, test or user code.

// System includes
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
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
#include "includes/stream_serializer.h"
#include "utilities/math_utils.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_elements/solid_elements/small_displacement.h"

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
constexpr double test_reference_temperature = 20.0;

/// Relative comparison tolerance for the restart-vs-reference match.
constexpr double restart_tolerance = 1.0e-8;

/// A serializable holder for a single element.
struct ElementHolder
{
    std::vector<Element::Pointer> elements;
    void save(Serializer& rSerializer) const { rSerializer.save("Elements", elements); }
    void load(Serializer& rSerializer) { rSerializer.load("Elements", elements); }
};

/// A serializable holder for a single constitutive law.
struct LawHolder
{
    ConstitutiveLaw::Pointer p_law;
    void save(Serializer& rSerializer) const { rSerializer.save("Law", p_law); }
    void load(Serializer& rSerializer) { rSerializer.load("Law", p_law); }
};

/// Builds a single-element model part (SMA SmallDisplacement) with the given
/// registered element name, 3D or 2D, and the given damage law.
Element::Pointer BuildDamageModel(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    ConstitutiveLaw::Pointer pLaw,
    const bool rIs3d,
    ModelPart*& rOutModelPart)
{
    KRATOS_TRY;
    ModelPart& r_mp = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = rIs3d ? 3 : 2;
    r_pi[SPACE_DIMENSION] = rIs3d ? 3 : 2;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY, &ACCELERATION,
                    &VOLUME_ACCELERATION, &TEMPERATURE, &NODAL_REFERENCE_TEMPERATURE,
                    &INITIAL_STRESS_TENSOR}) {
        r_mp.AddNodalSolutionStepVariable(*v);
    }

    const std::size_t number_of_nodes = rIs3d ? 8 : 4;
    const double c3[8][3] = {{0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,2},{2,0,2},{2,1,2},{0,1,2}};
    const double c2[4][3] = {{0,0,0},{2,0,0},{2,1,0},{0,1,0}};
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const double* c = rIs3d ? c3[i] : c2[i];
        Node::Pointer p = r_mp.CreateNewNode(i + 1, c[0], c[1], c[2]);
        p->AddDof(DISPLACEMENT_X); p->AddDof(DISPLACEMENT_Y); p->AddDof(DISPLACEMENT_Z);
        p->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        p->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
        p->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }
    auto p_props = r_mp.CreateNewProperties(1);
    (*p_props)[YOUNG_MODULUS] = test_young_modulus;
    (*p_props)[POISSON_RATIO] = test_poisson_ratio;
    (*p_props)[DENSITY] = test_density;
    (*p_props)[THERMAL_EXPANSION] = 1.0e-5;
    (*p_props)[DAMAGE_THRESHOLD] = 5.0e-3;
    (*p_props)[STRENGTH_RATIO] = 10.0;
    (*p_props)[FRACTURE_ENERGY] = 5000.0;
    p_props->SetValue(CONSTITUTIVE_LAW, pLaw);

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i)
        element_nodes[i] = i + 1;
    r_mp.CreateNewElement(rElementName, 1, element_nodes, p_props);
    Element::Pointer p_element = r_mp.pGetElement(1);
    p_element->Initialize(r_pi);
    rOutModelPart = &r_mp;
    return p_element;
    KRATOS_CATCH("");
}

/// Applies a uniform axial strain + temperature increment to the element and
/// records (cauchy_xx, damage, lhs_00) for the CURRENT committed state.
void ApplyAndRecord(Element& rElement, ModelPart& rMp, const double rEps, const double rDeltaT,
                    double& rCauchyXX, double& rDamage, double& rLhs00)
{
    ProcessInfo& r_pi = rMp.GetProcessInfo();
    const bool is3d = (rMp.GetProcessInfo()[DOMAIN_SIZE] == 3);
    for (auto& n : rElement.GetGeometry()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = rEps * x0[0];
        u[1] = -test_poisson_ratio * rEps * x0[1];
        u[2] = (is3d ? -test_poisson_ratio * rEps * x0[2] : 0.0);
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + rDeltaT;
    }

    // Trial response for the CURRENT committed state.
    std::vector<Vector> stress;
    rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress, r_pi);
    rCauchyXX = stress[0][0];

    std::vector<ConstitutiveLaw::Pointer> laws;
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, r_pi);
    double d = 0.0;
    laws[0]->GetValue(DAMAGE_VARIABLE, d);
    rDamage = d;

    Matrix lhs;
    rElement.CalculateLeftHandSide(lhs, r_pi);
    rLhs00 = lhs(0, 0);

    // Commit the step.
    r_pi[IS_CONVERGED] = true;
    rElement.FinalizeSolutionStep(r_pi);
}

/// Serializes an element to a binary archive string.
std::string SerializeElement(Element& rElement)
{
    StreamSerializer serializer;
    serializer.save("Holder", ElementHolder{{&rElement}});
    return serializer.GetStringRepresentation();
}

/// Loads an element from a binary archive string.
Element::Pointer LoadElement(const std::string& rArchive, ProcessInfo& rPi)
{
    StreamSerializer loader(rArchive);
    loader.SetLoadState();
    ElementHolder loaded;
    loader.load("Holder", loaded);
    Element::Pointer p = loaded.elements[0];
    for (auto& n : p->GetGeometry()) {
        n.AddDof(DISPLACEMENT_X); n.AddDof(DISPLACEMENT_Y); n.AddDof(DISPLACEMENT_Z);
    }
    rPi[DOMAIN_SIZE] = p->GetGeometry().WorkingSpaceDimension();
    rPi[SPACE_DIMENSION] = p->GetGeometry().WorkingSpaceDimension();
    rPi[DELTA_TIME] = 1.0;
    rPi[IS_CONVERGED] = true;
    return p;
}

/// Returns a 3D uniaxial-stress total-strain Voigt vector for the axial strain.
Vector UniaxialTotalStrainVector(const ProcessInfo& rPi, const double rEps)
{
    Vector e(6);
    e[0] = rEps;
    e[1] = -test_poisson_ratio * rEps;
    e[2] = -test_poisson_ratio * rEps;
    e[3] = 0.0; e[4] = 0.0; e[5] = 0.0;
    return e;
}

/// Applies a displacement field (no temperature change) to an element.
void ApplyDisplacement(Element& rElement, ModelPart& rMp, const double rEps)
{
    const bool is3d = (rMp.GetProcessInfo()[DOMAIN_SIZE] == 3);
    for (auto& n : rElement.GetGeometry()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = rEps * x0[0];
        u[1] = -test_poisson_ratio * rEps * x0[1];
        u[2] = (is3d ? -test_poisson_ratio * rEps * x0[2] : 0.0);
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
    }
}

/// Creates a wrapper model part owning the loaded element's geometry nodes and
/// carrying the given ProcessInfo (so solution-step data is accessible).
ModelPart& WrapLoadedElement(Model& rModel, Element& rElement, const ProcessInfo& rPi)
{
    ModelPart& r_wrap = rModel.CreateModelPart("Wrap", 2);
    for (std::size_t i = 0; i < rElement.GetGeometry().size(); ++i) {
        r_wrap.AddNode(rElement.GetGeometry().pGetPoint(i));
    }
    r_wrap.GetProcessInfo() = rPi;
    return r_wrap;
}

/// Reads the committed damage of an element.
double ElementDamage(Element& rElement, const ProcessInfo& rPi)
{
    std::vector<ConstitutiveLaw::Pointer> laws;
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, rPi);
    double d = 0.0;
    laws[0]->GetValue(DAMAGE_VARIABLE, d);
    return d;
}

} // namespace

//************************************************************************************
// 1. Local 3D damage restart continuation (reference vs restart branches).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_LocalDamage_3D_RestartContinuation, KratosDamFastSuite)
{
    // Steps: elastic, damage, unload, reload-below-max, reload-beyond-max,
    // further damage with temperature.
    const std::vector<std::pair<double, double>> steps = {
        {1.0e-6, 20.0}, {2.0e-5, 20.0}, {1.0e-5, 20.0}, {1.5e-5, 20.0},
        {3.0e-5, 20.0}, {3.5e-5, 50.0}};
    const std::size_t restart_after = 3;   // steps 1..3 committed before restart

    // Reference branch (uninterrupted).
    Model ref_model;
    ModelPart* p_ref_mp = nullptr;
    Element::Pointer p_ref = BuildDamageModel(ref_model, "Ref", "SmallDisplacementSolidElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_ref_mp);
    std::vector<double> ref_cauchy, ref_damage, ref_lhs;
    for (const auto& s : steps) {
        double c, d, l;
        ApplyAndRecord(*p_ref, *p_ref_mp, s.first, s.second, c, d, l);
        ref_cauchy.push_back(c); ref_damage.push_back(d); ref_lhs.push_back(l);
    }

    // Restart branch (serialize after step 3, destroy, reload, continue).
    Model rst_model;
    ModelPart* p_rst_mp = nullptr;
    Element::Pointer p_rst = BuildDamageModel(rst_model, "Rst", "SmallDisplacementSolidElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_rst_mp);
    for (std::size_t i = 0; i < restart_after; ++i) {
        double c, d, l;
        ApplyAndRecord(*p_rst, *p_rst_mp, steps[i].first, steps[i].second, c, d, l);
    }
    const std::string archive = SerializeElement(*p_rst);
    ProcessInfo rst_pi;
    Element::Pointer p_rst_loaded = LoadElement(archive, rst_pi);

    // NO manual material repair here.
    std::vector<double> rst_cauchy, rst_damage, rst_lhs;
    for (std::size_t i = restart_after; i < steps.size(); ++i) {
        double c, d, l;
        // Rebuild a model-part wrapper for the loaded element's geometry nodes so
        // the solution-step data is accessible.
        // (The loaded element owns its nodes; ApplyAndRecord uses the element's
        // geometry directly, so we pass a lightweight model part holding them.)
        Model wrapper_model;
        ModelPart& r_wrap = WrapLoadedElement(wrapper_model, *p_rst_loaded, rst_pi);
        ApplyAndRecord(*p_rst_loaded, r_wrap, steps[i].first, steps[i].second, c, d, l);
        rst_cauchy.push_back(c); rst_damage.push_back(d); rst_lhs.push_back(l);
    }

    // Compare restarted evolution against the uninterrupted reference.
    for (std::size_t i = 0; i < rst_cauchy.size(); ++i) {
        const std::size_t ref_idx = restart_after + i;
        KRATOS_EXPECT_NEAR(rst_cauchy[i], ref_cauchy[ref_idx],
                           restart_tolerance * std::max(1.0, std::abs(ref_cauchy[ref_idx])));
        KRATOS_EXPECT_NEAR(rst_damage[i], ref_damage[ref_idx], 1.0e-12);
        KRATOS_EXPECT_NEAR(rst_lhs[i], ref_lhs[ref_idx],
                           restart_tolerance * std::max(1.0, std::abs(ref_lhs[ref_idx])));
    }
    std::cout << "[6C.1] local 3D restart: continuation matches reference ("
              << rst_cauchy.size() << " steps, damage after restart = "
              << rst_damage.back() << ")" << std::endl;
}

//************************************************************************************
// 2. Local 2D plane-strain restart continuation.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_LocalDamage_2DPlaneStrain_RestartContinuation, KratosDamFastSuite)
{
    const std::vector<std::pair<double, double>> steps = {
        {1.0e-6, 20.0}, {2.0e-5, 20.0}, {1.0e-5, 20.0}, {3.0e-5, 50.0}};
    const std::size_t restart_after = 2;

    Model ref_model;
    ModelPart* p_ref_mp = nullptr;
    Element::Pointer p_ref = BuildDamageModel(ref_model, "Ref2D", "SmallDisplacementSolidElement2D4N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamagePlaneStrain2DLaw()), false, p_ref_mp);
    std::vector<double> ref_cauchy, ref_damage;
    for (const auto& s : steps) {
        double c, d, l;
        ApplyAndRecord(*p_ref, *p_ref_mp, s.first, s.second, c, d, l);
        ref_cauchy.push_back(c); ref_damage.push_back(d);
    }

    Model rst_model;
    ModelPart* p_rst_mp = nullptr;
    Element::Pointer p_rst = BuildDamageModel(rst_model, "Rst2D", "SmallDisplacementSolidElement2D4N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamagePlaneStrain2DLaw()), false, p_rst_mp);
    for (std::size_t i = 0; i < restart_after; ++i) {
        double c, d, l;
        ApplyAndRecord(*p_rst, *p_rst_mp, steps[i].first, steps[i].second, c, d, l);
    }
    const std::string archive = SerializeElement(*p_rst);
    ProcessInfo rst_pi;
    Element::Pointer p_rst_loaded = LoadElement(archive, rst_pi);
    for (std::size_t i = restart_after; i < steps.size(); ++i) {
        double c, d, l;
        Model wrapper_model;
        ModelPart& r_wrap = WrapLoadedElement(wrapper_model, *p_rst_loaded, rst_pi);
        ApplyAndRecord(*p_rst_loaded, r_wrap, steps[i].first, steps[i].second, c, d, l);
        KRATOS_EXPECT_NEAR(c, ref_cauchy[restart_after + (i - restart_after)],
                           restart_tolerance * std::max(1.0, std::abs(c)));
        KRATOS_EXPECT_NEAR(d, ref_damage[restart_after + (i - restart_after)], 1.0e-12);
    }
    std::cout << "[6C.1] local 2D plane-strain restart: continuation matches reference" << std::endl;
}

//************************************************************************************
// 3. Direct law-only binary serialization restart.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_DirectLaw_Serialization_Restart, KratosDamFastSuite)
{
    // Law-only round trip: commit damage, serialize, load, first response with
    // valid Parameters, continuation. No manual rebind.
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = BuildDamageModel(model, "LawOnly", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_mp);
    ProcessInfo& r_pi = p_mp->GetProcessInfo();
    double c, d, l;
    ApplyAndRecord(*p_elem, *p_mp, 1.0e-6, 20.0, c, d, l);   // elastic
    ApplyAndRecord(*p_elem, *p_mp, 2.0e-5, 20.0, c, d, l);   // damages
    const double committed = ElementDamage(*p_elem, r_pi);
    KRATOS_EXPECT_GT(committed, 0.0);

    std::vector<ConstitutiveLaw::Pointer> laws;
    p_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, r_pi);

    LawHolder holder;
    holder.p_law = laws[0];
    StreamSerializer serializer;
    serializer.save("Holder", holder);
    const std::string archive = serializer.GetStringRepresentation();

    LawHolder loaded;
    StreamSerializer loader(archive);
    loader.SetLoadState();
    loader.load("Holder", loaded);

    // Build valid ConstitutiveLaw::Parameters for a first response.
    const auto& r_geom = p_elem->GetGeometry();
    Vector N(r_geom.size());
    r_geom.ShapeFunctionsValues(N, r_geom.IntegrationPoints(p_elem->GetIntegrationMethod())[0].Coordinates());
    ConstitutiveLaw::Parameters params(r_geom, *p_elem->pGetProperties(), r_pi);
    params.SetShapeFunctionsValues(N);
    // A tensile-dominant total strain so the thermal-corrected mechanical
    // strain produces positive stress in the x direction.
    Vector total_strain(6);
    noalias(total_strain) = ZeroVector(6);
    total_strain[0] = 5.0e-4;
    params.SetStrainVector(total_strain);
    Vector stress(6);
    Matrix D(6, 6);
    params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    params.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    params.SetStressVector(stress);
    params.SetConstitutiveMatrix(D);

    // FIRST response after load - NO manual material repair.
    loaded.p_law->CalculateMaterialResponseCauchy(params);
    double d_after = 0.0;
    loaded.p_law->GetValue(DAMAGE_VARIABLE, d_after);
    KRATOS_EXPECT_GT(stress[0], 0.0);
    KRATOS_EXPECT_GT(d_after, 0.0);

    std::cout << "[6C.1] direct law restart: first response works, damage=" << d_after << std::endl;
}

//************************************************************************************
// 4. Nonlocal Simo-Ju restart continuation.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_Nonlocal_SimoJu_RestartContinuation, KratosDamFastSuite)
{
    const std::vector<std::pair<double, double>> steps = {
        {1.0e-6, 20.0}, {2.0e-5, 20.0}, {1.0e-5, 20.0}, {3.0e-5, 50.0}};
    const std::size_t restart_after = 2;

    auto run_branch = [&](Model& rModel, const std::string& rName, const std::size_t rStart) {
        ModelPart* p_mp = nullptr;
        Element::Pointer p = BuildDamageModel(rModel, rName, "SmallDisplacementElement3D8N",
            ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw()), true, p_mp);
        std::vector<double> cauchy, damage;
        for (std::size_t i = rStart; i < steps.size(); ++i) {
            // Prescribe the nonlocal driving through the normal workflow.
            std::vector<ConstitutiveLaw::Pointer> laws;
            p->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, p_mp->GetProcessInfo());
            laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.0e-2 * (1.0 + i), p_mp->GetProcessInfo());
            double c, d, l;
            ApplyAndRecord(*p, *p_mp, steps[i].first, steps[i].second, c, d, l);
            cauchy.push_back(c); damage.push_back(d);
        }
        return std::pair<std::vector<double>, std::vector<double>>{cauchy, damage};
    };

    Model ref_model;
    auto ref = run_branch(ref_model, "RefNLSJ", 0);

    Model rst_model;
    ModelPart* p_rst_mp = nullptr;
    Element::Pointer p_rst = BuildDamageModel(rst_model, "RstNLSJ", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw()), true, p_rst_mp);
    for (std::size_t i = 0; i < restart_after; ++i) {
        std::vector<ConstitutiveLaw::Pointer> laws;
        p_rst->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, p_rst_mp->GetProcessInfo());
        laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.0e-2 * (1.0 + i), p_rst_mp->GetProcessInfo());
        double c, d, l;
        ApplyAndRecord(*p_rst, *p_rst_mp, steps[i].first, steps[i].second, c, d, l);
    }
    const std::string archive = SerializeElement(*p_rst);
    ProcessInfo rst_pi;
    Element::Pointer p_rst_loaded = LoadElement(archive, rst_pi);
    for (std::size_t i = restart_after; i < steps.size(); ++i) {
        std::vector<ConstitutiveLaw::Pointer> laws;
        p_rst_loaded->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, rst_pi);
        laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.0e-2 * (1.0 + i), rst_pi);
        double c, d, l;
        Model wrapper_model;
        ModelPart& r_wrap = WrapLoadedElement(wrapper_model, *p_rst_loaded, rst_pi);
        ApplyAndRecord(*p_rst_loaded, r_wrap, steps[i].first, steps[i].second, c, d, l);
        const std::size_t ref_idx = restart_after + (i - restart_after);
        KRATOS_EXPECT_NEAR(c, ref.first[ref_idx], restart_tolerance * std::max(1.0, std::abs(c)));
        KRATOS_EXPECT_NEAR(d, ref.second[ref_idx], 1.0e-12);
    }
    std::cout << "[6C.1] nonlocal Simo-Ju restart: continuation matches reference" << std::endl;
}

//************************************************************************************
// 5. Nonlocal Modified-Mises restart continuation (same rebinding path).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_Nonlocal_ModifiedMises_RestartContinuation, KratosDamFastSuite)
{
    // Modified-Mises uses the same flow/yield/hardening transient Properties
    // dependency and the same automatic rebinding (it shares the nonlocal
    // ReinitializeMaterialProperties path); a focused representative test is
    // sufficient (justified inheritance coverage).
    const std::vector<std::pair<double, double>> steps = {
        {1.0e-6, 20.0}, {2.0e-5, 20.0}, {3.0e-5, 50.0}};
    const std::size_t restart_after = 1;

    Model ref_model;
    ModelPart* p_ref_mp = nullptr;
    Element::Pointer p_ref = BuildDamageModel(ref_model, "RefNLMM", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamage3DLaw()), true, p_ref_mp);
    std::vector<double> ref_cauchy, ref_damage;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        std::vector<ConstitutiveLaw::Pointer> laws;
        p_ref->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, p_ref_mp->GetProcessInfo());
        laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.0e-2 * (1.0 + i), p_ref_mp->GetProcessInfo());
        double c, d, l;
        ApplyAndRecord(*p_ref, *p_ref_mp, steps[i].first, steps[i].second, c, d, l);
        ref_cauchy.push_back(c); ref_damage.push_back(d);
    }

    Model rst_model;
    ModelPart* p_rst_mp = nullptr;
    Element::Pointer p_rst = BuildDamageModel(rst_model, "RstNLMM", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamage3DLaw()), true, p_rst_mp);
    for (std::size_t i = 0; i < restart_after; ++i) {
        std::vector<ConstitutiveLaw::Pointer> laws;
        p_rst->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, p_rst_mp->GetProcessInfo());
        laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.0e-2 * (1.0 + i), p_rst_mp->GetProcessInfo());
        double c, d, l;
        ApplyAndRecord(*p_rst, *p_rst_mp, steps[i].first, steps[i].second, c, d, l);
    }
    const std::string archive = SerializeElement(*p_rst);
    ProcessInfo rst_pi;
    Element::Pointer p_rst_loaded = LoadElement(archive, rst_pi);
    for (std::size_t i = restart_after; i < steps.size(); ++i) {
        std::vector<ConstitutiveLaw::Pointer> laws;
        p_rst_loaded->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, rst_pi);
        laws[0]->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.0e-2 * (1.0 + i), rst_pi);
        double c, d, l;
        Model wrapper_model;
        ModelPart& r_wrap = WrapLoadedElement(wrapper_model, *p_rst_loaded, rst_pi);
        ApplyAndRecord(*p_rst_loaded, r_wrap, steps[i].first, steps[i].second, c, d, l);
        const std::size_t ref_idx = restart_after + (i - restart_after);
        KRATOS_EXPECT_NEAR(c, ref_cauchy[ref_idx], restart_tolerance * std::max(1.0, std::abs(c)));
        KRATOS_EXPECT_NEAR(d, ref_damage[ref_idx], 1.0e-12);
    }
    std::cout << "[6C.1] nonlocal Modified-Mises restart: continuation matches reference" << std::endl;
}

//************************************************************************************
// 6. Multi-Properties isolation after restart.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_MultiProperties_Isolation, KratosDamFastSuite)
{
    // Two laws with materially different Properties must be rebound to their OWN
    // current Properties after restart (no cross-binding).
    Model model_a, model_b;
    ModelPart* p_mp_a = nullptr;
    ModelPart* p_mp_b = nullptr;
    Element::Pointer p_a = BuildDamageModel(model_a, "PropA", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_mp_a);
    Element::Pointer p_b = BuildDamageModel(model_b, "PropB", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_mp_b);

    // Change material parameters on B.
    (*p_b->pGetProperties())[YOUNG_MODULUS] = 4.0e7;
    (*p_b->pGetProperties())[DAMAGE_THRESHOLD] = 3.0e-3;

    double c_a, d_a, l_a;
    ApplyAndRecord(*p_a, *p_mp_a, 2.0e-5, 20.0, c_a, d_a, l_a);
    double c_b, d_b, l_b;
    ApplyAndRecord(*p_b, *p_mp_b, 2.0e-5, 20.0, c_b, d_b, l_b);
    // Materially different damage/response.
    KRATOS_EXPECT_NE(d_a, d_b);

    // Restart both.
    ProcessInfo pi_a, pi_b;
    Element::Pointer p_a_loaded = LoadElement(SerializeElement(*p_a), pi_a);
    Element::Pointer p_b_loaded = LoadElement(SerializeElement(*p_b), pi_b);

    // Continue both; each must use ITS OWN Properties.
    double c_a2, d_a2, l_a2;
    {
        Model w; ModelPart& rw = WrapLoadedElement(w, *p_a_loaded, pi_a);
        ApplyAndRecord(*p_a_loaded, rw, 3.0e-5, 20.0, c_a2, d_a2, l_a2);
    }
    double c_b2, d_b2, l_b2;
    {
        Model w; ModelPart& rw = WrapLoadedElement(w, *p_b_loaded, pi_b);
        ApplyAndRecord(*p_b_loaded, rw, 3.0e-5, 20.0, c_b2, d_b2, l_b2);
    }

    // Reference: uninterrupted A and B through the same extra step.
    Model ref_a, ref_b;
    ModelPart* p_ra_mp = nullptr; ModelPart* p_rb_mp = nullptr;
    Element::Pointer p_ra = BuildDamageModel(ref_a, "RA", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_ra_mp);
    Element::Pointer p_rb = BuildDamageModel(ref_b, "RB", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_rb_mp);
    (*p_rb->pGetProperties())[YOUNG_MODULUS] = 4.0e7;
    (*p_rb->pGetProperties())[DAMAGE_THRESHOLD] = 3.0e-3;
    double ra_c, ra_d, ra_l, rb_c, rb_d, rb_l;
    ApplyAndRecord(*p_ra, *p_ra_mp, 2.0e-5, 20.0, ra_c, ra_d, ra_l);
    ApplyAndRecord(*p_rb, *p_rb_mp, 2.0e-5, 20.0, rb_c, rb_d, rb_l);
    ApplyAndRecord(*p_ra, *p_ra_mp, 3.0e-5, 20.0, ra_c, ra_d, ra_l);
    ApplyAndRecord(*p_rb, *p_rb_mp, 3.0e-5, 20.0, rb_c, rb_d, rb_l);

    KRATOS_EXPECT_NEAR(c_a2, ra_c, restart_tolerance * std::max(1.0, std::abs(ra_c)));
    KRATOS_EXPECT_NEAR(d_a2, ra_d, 1.0e-12);
    KRATOS_EXPECT_NEAR(c_b2, rb_c, restart_tolerance * std::max(1.0, std::abs(rb_c)));
    KRATOS_EXPECT_NEAR(d_b2, rb_d, 1.0e-12);
    KRATOS_EXPECT_NE(d_a2, d_b2);   // the two materials remain distinct after restart
    std::cout << "[6C.1] multi-Properties isolation: each restarted law rebound to "
              << "its OWN current Properties (dA=" << d_a2 << ", dB=" << d_b2 << ")" << std::endl;
}

//************************************************************************************
// 7. Specialized outputs right after restart (no manual repair).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_SpecializedOutputs_AfterRestart, KratosDamFastSuite)
{
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = BuildDamageModel(model, "SpecOut", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_mp);
    double c, d, l;
    ApplyAndRecord(*p_elem, *p_mp, 1.0e-6, 40.0, c, d, l);
    ApplyAndRecord(*p_elem, *p_mp, 2.0e-5, 40.0, c, d, l);
    const double committed_damage = ElementDamage(*p_elem, p_mp->GetProcessInfo());
    KRATOS_EXPECT_GT(committed_damage, 0.0);

    ProcessInfo rst_pi;
    Element::Pointer p_loaded = LoadElement(SerializeElement(*p_elem), rst_pi);

    // Immediately after deserialization request the specialized outputs.
    std::vector<Vector> th_strain, th_stress, mech_stress;
    p_loaded->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, th_strain, rst_pi);
    p_loaded->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, th_stress, rst_pi);
    p_loaded->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, mech_stress, rst_pi);
    std::vector<Matrix> th_tensor;
    p_loaded->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, th_tensor, rst_pi);
    KRATOS_EXPECT_EQ(th_strain[0].size(), 6u);
    KRATOS_EXPECT_GT(th_strain[0][0], 0.0);   // thermal strain from the serialized temperature field
    // Decomposition: total == mech - therm.
    std::vector<Vector> cauchy;
    p_loaded->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy, rst_pi);
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_NEAR(cauchy[0](i), mech_stress[0](i) - th_stress[0](i), 1.0e-9);

    // No state commit / no damage change from the outputs.
    const double damage_after = ElementDamage(*p_loaded, rst_pi);
    KRATOS_EXPECT_NEAR(damage_after, committed_damage, 1.0e-12);
    std::cout << "[6C.1] specialized outputs after restart work without manual repair" << std::endl;
}

//************************************************************************************
// 8. Clone after restart + repeated-call determinism.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_CloneAndDeterminism_AfterRestart, KratosDamFastSuite)
{
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = BuildDamageModel(model, "CloneRst", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_mp);
    double c, d, l;
    ApplyAndRecord(*p_elem, *p_mp, 1.0e-6, 20.0, c, d, l);
    ApplyAndRecord(*p_elem, *p_mp, 2.0e-5, 20.0, c, d, l);
    const double committed_damage = ElementDamage(*p_elem, p_mp->GetProcessInfo());
    KRATOS_EXPECT_GT(committed_damage, 0.0);

    ProcessInfo rst_pi;
    Element::Pointer p_loaded = LoadElement(SerializeElement(*p_elem), rst_pi);

    // Clone the loaded law.
    std::vector<ConstitutiveLaw::Pointer> laws;
    p_loaded->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, rst_pi);
    ConstitutiveLaw::Pointer p_clone = laws[0]->Clone();

    // Repeated responses on the original (3x) are identical and read-only.
    std::vector<double> c0;
    for (int rep = 0; rep < 3; ++rep) {
        std::vector<Vector> stress;
        p_loaded->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress, rst_pi);
        c0.push_back(stress[0][0]);
    }
    KRATOS_EXPECT_NEAR(c0[0], c0[1], 1.0e-12);
    KRATOS_EXPECT_NEAR(c0[1], c0[2], 1.0e-12);
    KRATOS_EXPECT_NEAR(ElementDamage(*p_loaded, rst_pi), committed_damage, 1.0e-12);

    // Clone is deterministic and shares no corrupt transient state.
    double d_clone = 0.0;
    p_clone->GetValue(DAMAGE_VARIABLE, d_clone);
    KRATOS_EXPECT_NEAR(d_clone, committed_damage, 1.0e-12);
    std::cout << "[6C.1] clone-after-restart + repeated-call determinism OK" << std::endl;
}

//************************************************************************************
// 9. Trial / rollback / commit after restart.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6C1_TrialRollbackCommit_AfterRestart, KratosDamFastSuite)
{
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = BuildDamageModel(model, "Rollback", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_mp);
    double c, d, l;
    ApplyAndRecord(*p_elem, *p_mp, 1.0e-6, 20.0, c, d, l);
    ApplyAndRecord(*p_elem, *p_mp, 2.0e-5, 20.0, c, d, l);
    const double committed_before = ElementDamage(*p_elem, p_mp->GetProcessInfo());

    ProcessInfo rst_pi;
    Element::Pointer p_loaded = LoadElement(SerializeElement(*p_elem), rst_pi);

    // Trial response (loading beyond the committed max).
    Model w; ModelPart& rw = WrapLoadedElement(w, *p_loaded, rst_pi);
    ApplyDisplacement(*p_loaded, rw, 3.0e-5);
    std::vector<Vector> stress;
    p_loaded->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress, rst_pi);

    // Non-converged finalize -> rollback, committed damage unchanged.
    rst_pi[IS_CONVERGED] = false;
    p_loaded->FinalizeSolutionStep(rst_pi);
    KRATOS_EXPECT_NEAR(ElementDamage(*p_loaded, rst_pi), committed_before, 1.0e-12);

    // Converged finalize -> commit advances damage exactly once.
    rst_pi[IS_CONVERGED] = true;
    p_loaded->FinalizeSolutionStep(rst_pi);
    const double committed_after = ElementDamage(*p_loaded, rst_pi);
    KRATOS_EXPECT_GT(committed_after, committed_before);

    // Reference rollback/commit sequence (uninterrupted).
    Model ref_model;
    ModelPart* p_ref_mp = nullptr;
    Element::Pointer p_ref = BuildDamageModel(ref_model, "RefRB", "SmallDisplacementElement3D8N",
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), true, p_ref_mp);
    double rc, rd, rl;
    ApplyAndRecord(*p_ref, *p_ref_mp, 1.0e-6, 20.0, rc, rd, rl);
    ApplyAndRecord(*p_ref, *p_ref_mp, 2.0e-5, 20.0, rc, rd, rl);
    KRATOS_EXPECT_NEAR(ElementDamage(*p_ref, p_ref_mp->GetProcessInfo()), committed_before, 1.0e-12);
    ApplyDisplacement(*p_ref, *p_ref_mp, 3.0e-5);
    p_ref_mp->GetProcessInfo()[IS_CONVERGED] = false;
    p_ref->FinalizeSolutionStep(p_ref_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(ElementDamage(*p_ref, p_ref_mp->GetProcessInfo()), committed_before, 1.0e-12);
    p_ref_mp->GetProcessInfo()[IS_CONVERGED] = true;
    p_ref->FinalizeSolutionStep(p_ref_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(ElementDamage(*p_ref, p_ref_mp->GetProcessInfo()), committed_after, 1.0e-12);
    std::cout << "[6C.1] trial/rollback/commit after restart matches reference" << std::endl;
}

} // namespace Testing
} // namespace Kratos
