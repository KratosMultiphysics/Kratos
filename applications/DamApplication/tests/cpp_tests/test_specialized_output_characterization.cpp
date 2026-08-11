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

// Phase 5B.1: characterize the specialized thermo-mechanical integration-point
// output API (THERMAL_STRAIN/STRESS_VECTOR/TENSOR, MECHANICAL_STRESS_VECTOR/TENSOR)
// still implemented in SmallDisplacementThermoMechanicElement. Characterization
// only: no production code is modified.

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
#include "geometries/hexahedra_3d_8.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"

// StructuralMechanicsApplication small-displacement element
#include "custom_elements/solid_elements/small_displacement.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerance (characterization: relative agreement with an absolute
/// floor for near-zero components).
constexpr double characterization_relative_tolerance = 1.0e-6;
constexpr double characterization_absolute_floor = 1.0e-10;

/// Comparison helper: relative tolerance with absolute floor.
bool NearComponent(const double rValue, const double rReference)
{
    return std::abs(rValue - rReference) <=
           std::max(characterization_relative_tolerance * std::abs(rReference),
                    characterization_absolute_floor);
}

/// Material data.
constexpr double test_poisson_ratio = 0.2;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;

/// Test-only element subclasses exposing the constitutive-law vector.
class TestThermoMechanicElement : public SmallDisplacementThermoMechanicElement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TestThermoMechanicElement);
    using BaseType = SmallDisplacementThermoMechanicElement;
    TestThermoMechanicElement(IndexType NewId, GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties) {}
    ConstitutiveLaw& GetConstitutiveLaw(std::size_t i) { return *mConstitutiveLawVector[i]; }
    ConstitutiveLaw::Pointer GetConstitutiveLawPointer(std::size_t i) { return mConstitutiveLawVector[i]; }
};

class TestSmallDisplacementElement : public SmallDisplacement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TestSmallDisplacementElement);
    using BaseType = SmallDisplacement;
    TestSmallDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry,
                                 PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties) {}
    ConstitutiveLaw& GetConstitutiveLaw(std::size_t i) { return *mConstitutiveLawVector[i]; }
    ConstitutiveLaw::Pointer GetConstitutiveLawPointer(std::size_t i) { return mConstitutiveLawVector[i]; }
};

/// Builds a single H8 model part with the given law and element type.
template<class TTestElement>
ModelPart& CreateOutputModelPart(
    Model& rModel,
    const std::string& rName,
    TTestElement*& rOutElement,
    ConstitutiveLaw::Pointer pLaw,
    const bool rNodalProperties)
{
    KRATOS_TRY;
    ModelPart& r_model_part = rModel.CreateModelPart(rName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY, &ACCELERATION,
                    &VOLUME_ACCELERATION, &TEMPERATURE, &NODAL_REFERENCE_TEMPERATURE,
                    &NODAL_YOUNG_MODULUS, &NODAL_CAUCHY_STRESS_TENSOR, &NODAL_AREA,
                    &INITIAL_STRESS_TENSOR}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }
    const double coords[8][3] = {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}};
    for (std::size_t i = 0; i < 8; ++i) {
        Node::Pointer n = r_model_part.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        n->AddDof(DISPLACEMENT_X); n->AddDof(DISPLACEMENT_Y); n->AddDof(DISPLACEMENT_Z);
        n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        n->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = 2.0e7;
        Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
        n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }
    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = 2.0e7;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = 5.0e-3;
    (*p_prop)[STRENGTH_RATIO] = 10.0;
    (*p_prop)[FRACTURE_ENERGY] = 5000.0;
    p_prop->SetValue(CONSTITUTIVE_LAW, pLaw);
    Geometry<Node>::PointsArrayType pts;
    for (std::size_t i = 0; i < 8; ++i) pts.push_back(r_model_part.pGetNode(i + 1));
    auto p_elem = Kratos::make_intrusive<TTestElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pts)), p_prop);
    r_model_part.AddElement(p_elem);
    rOutElement = p_elem.get();
    return r_model_part;
    KRATOS_CATCH("");
}

/// Applies a uniaxial-STRESS state with uniform temperature change.
void ApplyState(ModelPart& rModelPart, const double rEpsilonX, const double rDeltaTemperature)
{
    for (auto& n : rModelPart.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = rEpsilonX * x0[0];
        u[1] = -test_poisson_ratio * rEpsilonX * x0[1];
        u[2] = -test_poisson_ratio * rEpsilonX * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + rDeltaTemperature;
    }
}

/// Computes the elastic C matrix (3D).
Matrix C3D(const double E, const double nu)
{
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    Matrix C(6, 6);
    noalias(C) = ZeroMatrix(6, 6);
    C(0,0)=lambda+2*mu; C(0,1)=lambda;       C(0,2)=lambda;
    C(1,0)=lambda;       C(1,1)=lambda+2*mu; C(1,2)=lambda;
    C(2,0)=lambda;       C(2,1)=lambda;       C(2,2)=lambda+2*mu;
    C(3,3)=mu; C(4,4)=mu; C(5,5)=mu;
    return C;
}

/// Returns the 3D total-strain Voigt vector for the applied uniaxial-STRESS state.
Vector TotalStrain3D(const double eps)
{
    Vector e(6);
    e[0] = eps; e[1] = -test_poisson_ratio * eps; e[2] = -test_poisson_ratio * eps;
    e[3] = 0.0; e[4] = 0.0; e[5] = 0.0;
    return e;
}

} // namespace

//************************************************************************************
// 1. Non-nodal linear reference: what does the legacy element output return?
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SpecializedOutputs_LinearReference, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestThermoMechanicElement>(
        model, "RefLinear", p_elem, ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw()), false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    const double eps = 2.0e-6, dT = 40.0;
    ApplyState(r_mp, eps, dT);

    const double E = 2.0e7;
    const Matrix C = C3D(E, test_poisson_ratio);
    const Vector e_total = TotalStrain3D(eps);
    const double alpha_dT = test_thermal_expansion * dT;
    Vector e_th(6);
    e_th[0] = alpha_dT; e_th[1] = alpha_dT; e_th[2] = alpha_dT; e_th[3] = 0; e_th[4] = 0; e_th[5] = 0;
    const Vector mech_stress = prod(C, e_total);
    const Vector therm_stress = prod(C, e_th);
    const Vector total_stress = prod(C, e_total - e_th);

    std::vector<Vector> out;
    p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_TRUE(NearComponent(out[0](i), total_stress[i]));

    p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_TRUE(NearComponent(out[0](i), mech_stress[i]));

    p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_TRUE(NearComponent(out[0](i), therm_stress[i]));

    // THERMAL_STRAIN_VECTOR: legacy BUG - returns the TOTAL strain (the law
    // response is invoked without COMPUTE_STRESS / COMPUTE_CONSTITUTIVE_TENSOR /
    // VOLUMETRIC_TENSOR_ONLY, so it does nothing and the strain buffer keeps the
    // element-provided total strain). The intended value is epsilon_th.
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, out, r_mp.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_NEAR(out[0](i), e_total[i], 1.0e-12);  // returns total strain (documented bug)
    std::cout << "[5B.1] linear reference: total==mech-therm="
              << (std::abs((mech_stress[0]-therm_stress[0])-total_stress[0]) < 1.0e-8)
              << " THERMAL_STRAIN returns total strain (bug)" << std::endl;
}

//************************************************************************************
// 2. Nodal linear family: same legacy semantics, but with nodal-E interpolation.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SpecializedOutputs_NodalLinear, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestThermoMechanicElement>(
        model, "NodalOutput", p_elem, ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal()), true);
    p_elem->Initialize(r_mp.GetProcessInfo());
    const double eps = 2.0e-6, dT = 40.0;
    ApplyState(r_mp, eps, dT);
    // Non-uniform nodal E.
    for (auto& n : r_mp.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        n.FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = 2.0e7 + 1.0e6 * (x0[0] + x0[1]);
    }

    std::vector<Vector> out;
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(out[0].size(), 6);
    p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(out[0].size(), 6);
    // THERMAL_STRAIN_VECTOR returns the total strain (same legacy bug path).
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(out[0].size(), 6);
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_NEAR(out[0](i), TotalStrain3D(eps)[i], 1.0e-12);
    std::cout << "[5B.1] nodal linear: legacy outputs follow the linear path; THERMAL_STRAIN bug present" << std::endl;
}

//************************************************************************************
// 3. Local-damage: determine the damage-scaling semantics.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SpecializedOutputs_LocalDamage, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestThermoMechanicElement>(
        model, "LocalOut", p_elem, ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    const double eps = 2.0e-5, dT = 40.0;  // damaged state
    ApplyState(r_mp, eps, dT);

    // Commit a damaged state (converged).
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    p_elem->FinalizeSolutionStep(r_mp.GetProcessInfo());

    double d = 0.0;
    p_elem->GetConstitutiveLaw(0).GetValue(DAMAGE_VARIABLE, d);
    KRATOS_EXPECT_TRUE(d > 0.0);
    std::cout << "[5B.1] local damage d=" << d << std::endl;

    const Matrix C = C3D(2.0e7, test_poisson_ratio);
    const Vector e_total = TotalStrain3D(eps);
    const double alpha_dT = test_thermal_expansion * dT;
    Vector e_th(6);
    e_th[0]=alpha_dT; e_th[1]=alpha_dT; e_th[2]=alpha_dT; e_th[3]=0; e_th[4]=0; e_th[5]=0;

    // Semantic 1: undamaged decomposition (mech=C*e_total, therm=C*e_th).
    Vector mech_undam = prod(C, e_total);
    Vector therm_undam = prod(C, e_th);
    // Semantic 2: damage-scaled ((1-d)*C*e_total etc).
    Vector mech_damaged = (1.0 - d) * mech_undam;
    Vector therm_damaged = (1.0 - d) * therm_undam;
    // Total damaged stress.
    std::vector<Vector> out;
    p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    const Vector total_observed = out[0];
    const Vector total_damaged = (1.0 - d) * prod(C, e_total - e_th);
    KRATOS_EXPECT_NEAR(total_observed[0], total_damaged[0],
                       characterization_relative_tolerance * std::abs(total_observed[0]));

    p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    const Vector mech_observed = out[0];
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
    const Vector therm_observed = out[0];

    // Determine which semantic: compare mech_observed against mech_undam vs mech_damaged.
    const bool mech_is_damaged = std::abs(mech_observed[0] - mech_damaged[0]) <
                                 std::abs(mech_observed[0] - mech_undam[0]);
    const bool therm_is_damaged = std::abs(therm_observed[0] - therm_damaged[0]) <
                                  std::abs(therm_observed[0] - therm_undam[0]);
    std::cout << "[5B.1] local damage: mech_is_damage_scaled=" << mech_is_damaged
              << " therm_is_damage_scaled=" << therm_is_damaged
              << " total==mech-therm="
              << (std::abs((mech_observed[0]-therm_observed[0])-total_observed[0]) <
                  characterization_relative_tolerance * std::abs(total_observed[0])) << std::endl;
    KRATOS_EXPECT_TRUE(mech_is_damaged);
    KRATOS_EXPECT_TRUE(therm_is_damaged);
    KRATOS_EXPECT_TRUE(std::abs((mech_observed[0]-therm_observed[0])-total_observed[0]) <
                       characterization_relative_tolerance * std::abs(total_observed[0]));
}

//************************************************************************************
// 4. Nonlocal Simo-Ju and Modified-Mises: damage-scaled semantics + side effects.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SpecializedOutputs_NonlocalDamage, KratosDamFastSuite)
{
    for (bool modified_mises : {false, true}) {
        Model model;
        auto p_law = modified_mises
            ? ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamage3DLaw())
            : ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw());
        TestThermoMechanicElement* p_elem = nullptr;
        ModelPart& r_mp = CreateOutputModelPart<TestThermoMechanicElement>(
            model, modified_mises ? "ModOut" : "SJOut", p_elem, p_law, false);
        p_elem->Initialize(r_mp.GetProcessInfo());
        const double eps = 2.0e-5, dT = 40.0;
        ApplyState(r_mp, eps, dT);

        // Commit a damaged state driven by a prescribed nonlocal strain.
        r_mp.GetProcessInfo()[IS_CONVERGED] = true;
        p_elem->GetConstitutiveLaw(0).SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_mp.GetProcessInfo());
        p_elem->FinalizeSolutionStep(r_mp.GetProcessInfo());

        double d = 0.0;
        p_elem->GetConstitutiveLaw(0).GetValue(DAMAGE_VARIABLE, d);
        KRATOS_EXPECT_TRUE(d > 0.0);

        // Side-effect audit: record committed state before outputs.
        double threshold_before = 0.0, local_before = 0.0;
        p_elem->GetConstitutiveLaw(0).GetValue(STATE_VARIABLE, threshold_before);
        p_elem->GetConstitutiveLaw(0).GetValue(LOCAL_EQUIVALENT_STRAIN, local_before);

        std::vector<Vector> out;
        std::vector<Matrix> m_out;
        for (std::size_t repeat = 0; repeat < 3; ++repeat) {
            p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
            p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
            p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, out, r_mp.GetProcessInfo());
            p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, m_out, r_mp.GetProcessInfo());
        }
        double threshold_after = 0.0, local_after = 0.0;
        p_elem->GetConstitutiveLaw(0).GetValue(STATE_VARIABLE, threshold_after);
        p_elem->GetConstitutiveLaw(0).GetValue(LOCAL_EQUIVALENT_STRAIN, local_after);
        KRATOS_EXPECT_NEAR(threshold_after, threshold_before, 1.0e-15);
        KRATOS_EXPECT_NEAR(local_after, local_before, 1.0e-15);

        // Damage-scaled semantics.
        const Matrix C = C3D(2.0e7, test_poisson_ratio);
        const Vector e_total = TotalStrain3D(eps);
        const double alpha_dT = test_thermal_expansion * dT;
        Vector e_th(6);
        e_th[0]=alpha_dT; e_th[1]=alpha_dT; e_th[2]=alpha_dT; e_th[3]=0; e_th[4]=0; e_th[5]=0;
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, out, r_mp.GetProcessInfo());
        const Vector mech_observed = out[0];
        const Vector mech_damaged = (1.0 - d) * prod(C, e_total);
        KRATOS_EXPECT_NEAR(mech_observed[0], mech_damaged[0],
                           characterization_relative_tolerance * std::abs(mech_observed[0]));
        std::cout << "[5B.1] nonlocal " << (modified_mises ? "ModMises" : "SimoJu")
                  << ": d=" << d << " mech damage-scaled, no committed-state side effect" << std::endl;
    }
}

//************************************************************************************
// 5. SMA behavior for the specialized outputs (current laws, pre-Phase-5B.2).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SpecializedOutputs_SMABehavior, KratosDamFastSuite)
{
    // For the SMA element, the specialized outputs go through
    // CalculateOnConstitutiveLaw -> CalculateValue. Only the non-nodal linear law
    // implements CalculateValue for these variables; the others return the base
    // default (zero).
    for (bool linear : {true, false}) {
        Model model;
        auto p_law = linear
            ? ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw())
            : ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw());
        TestSmallDisplacementElement* p_elem = nullptr;
        ModelPart& r_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
            model, linear ? "SMAOutLinear" : "SMAOutDamage", p_elem, p_law, false);
        p_elem->Initialize(r_mp.GetProcessInfo());
        ApplyState(r_mp, 2.0e-6, 40.0);

        std::vector<Vector> out;
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, out, r_mp.GetProcessInfo());
        KRATOS_EXPECT_EQ(out.size(), 8);
        if (linear) {
            // The non-nodal linear law implements CalculateValue: correct eps_th.
            KRATOS_EXPECT_NEAR(out[0][0], test_thermal_expansion * 40.0, 1.0e-15);
        } else {
            // The damage law does not implement CalculateValue for this variable:
            // returns zero (base GetValue default).
            KRATOS_EXPECT_NEAR(out[0][0], 0.0, 1.0e-15);
        }
        std::cout << "[5B.1] SMA specialized output " << (linear ? "linear=correct" : "damage=zero") << std::endl;
    }
}

//************************************************************************************
// 6. Vector/tensor dimensions and tensor conversion.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SpecializedOutputs_Dimensions, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestThermoMechanicElement>(
        model, "DimOut", p_elem, ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw()), false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    ApplyState(r_mp, 2.0e-6, 40.0);

    std::vector<Vector> v_out;
    std::vector<Matrix> m_out;
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, v_out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(v_out[0].size(), 6);   // 3D Voigt size 6
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, v_out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(v_out[0].size(), 6);
    p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, v_out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(v_out[0].size(), 6);
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, m_out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(m_out[0].size1(), 3);
    KRATOS_EXPECT_EQ(m_out[0].size2(), 3);
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, m_out, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(m_out[0].size1(), 3);
    KRATOS_EXPECT_EQ(m_out[0].size2(), 3);
    std::cout << "[5B.1] dimensions: 3D Voigt 6, tensors 3x3" << std::endl;
}

} // namespace Testing
} // namespace Kratos
