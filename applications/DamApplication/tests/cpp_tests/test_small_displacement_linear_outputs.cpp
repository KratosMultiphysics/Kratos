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
class TestThermoMechanicElement : public SmallDisplacement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TestThermoMechanicElement);
    using BaseType = SmallDisplacement;
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

    // THERMAL_STRAIN_VECTOR: intentional bug fix (Phase 5B.2). The legacy
    // element used to return the TOTAL strain (the law response was invoked
    // without COMPUTE_STRESS / COMPUTE_CONSTITUTIVE_TENSOR / VOLUMETRIC_TENSOR_ONLY,
    // so the strain buffer kept the element-provided total strain). The
    // constitutive-law output now returns the actual thermal strain epsilon_th.
    p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, out, r_mp.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_EXPECT_TRUE(NearComponent(out[0](i), e_th[i]));  // returns epsilon_th (bug fixed)
    std::cout << "[5B.1] linear reference: total==mech-therm="
              << (std::abs((mech_stress[0]-therm_stress[0])-total_stress[0]) < 1.0e-8)
              << " THERMAL_STRAIN returns epsilon_th (bug fixed)" << std::endl;
}


//************************************************************************************
// 2. Nodal linear family: same legacy semantics, but with nodal-E interpolation.
//************************************************************************************


//************************************************************************************
// 3. Local-damage: determine the damage-scaling semantics.
//************************************************************************************


//************************************************************************************
// 4. Nonlocal Simo-Ju and Modified-Mises: damage-scaled semantics + side effects.
//************************************************************************************


//************************************************************************************
// 5. SMA behavior for the specialized outputs (current laws, pre-Phase-5B.2).
//************************************************************************************


//************************************************************************************
// 6. Vector/tensor dimensions and tensor conversion.
//************************************************************************************


} // namespace Testing
} // namespace Kratos
