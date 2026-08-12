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

// Phase 5B.2: the specialized thermo-mechanical Vector/Matrix integration-point
// outputs (THERMAL_STRAIN/STRESS_VECTOR/TENSOR, MECHANICAL_STRESS_VECTOR/TENSOR)
// are owned by the Dam constitutive laws through the parameter-aware
// CalculateValue path and are available through the generic element
// integration-point interface (StructuralMechanics SmallDisplacement and the
// legacy Dam SmallDisplacementThermoMechanicElement).

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
#include "geometries/quadrilateral_2d_4.h"
#include "utilities/math_utils.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress_nodal.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_stress_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"

// StructuralMechanicsApplication small-displacement element
#include "custom_elements/solid_elements/small_displacement.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerance (relative agreement with an absolute floor).
constexpr double characterization_relative_tolerance = 1.0e-6;
constexpr double characterization_absolute_floor = 1.0e-10;

/// Material data.
constexpr double test_poisson_ratio = 0.2;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_damage_threshold = 5.0e-3;

/// Comparison helper: relative tolerance with absolute floor.
bool Near(const double rValue, const double rReference)
{
    return std::abs(rValue - rReference) <=
           std::max(characterization_relative_tolerance * std::abs(rReference),
                    characterization_absolute_floor);
}

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
};

/// Test-only nonlocal law subclass exposing mNonlocalEquivalentStrain so the
/// tests can verify that the specialized outputs do not change it.
class TestNonlocalLaw : public ThermalSimoJuNonlocalDamage3DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TestNonlocalLaw);
    using BaseType = ThermalSimoJuNonlocalDamage3DLaw;
    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new TestNonlocalLaw(*this));
    }
    double NonlocalStrain() const { return mNonlocalEquivalentStrain; }
};

/// Geometry kind used by the tests.
enum class GeometryKind { Hexa3D, Quadrilateral2D };

/// Builds a single-element model part with the given law and element type.
template<class TElement>
ModelPart& CreateOutputModelPart(
    Model& rModel,
    const std::string& rName,
    TElement*& rOutElement,
    ConstitutiveLaw::Pointer pLaw,
    const GeometryKind rGeometryKind,
    const bool rNodalProperties)
{
    KRATOS_TRY;
    const bool is_3d = (rGeometryKind == GeometryKind::Hexa3D);
    ModelPart& r_model_part = rModel.CreateModelPart(rName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = is_3d ? 3 : 2;
    r_pi[SPACE_DIMENSION] = is_3d ? 3 : 2;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY, &ACCELERATION,
                    &VOLUME_ACCELERATION, &TEMPERATURE, &NODAL_REFERENCE_TEMPERATURE,
                    &NODAL_YOUNG_MODULUS, &NODAL_CAUCHY_STRESS_TENSOR, &NODAL_AREA,
                    &INITIAL_STRESS_TENSOR}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }

    const double coords_3d[8][3] = {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}};
    const double coords_2d[4][3] = {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0}};
    const std::size_t number_of_nodes = is_3d ? 8 : 4;
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const double* c = is_3d ? coords_3d[i] : coords_2d[i];
        Node::Pointer n = r_model_part.CreateNewNode(i + 1, c[0], c[1], c[2]);
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
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = 10.0;
    (*p_prop)[FRACTURE_ENERGY] = 5000.0;
    p_prop->SetValue(CONSTITUTIVE_LAW, pLaw);

    Geometry<Node>::PointsArrayType pts;
    for (std::size_t i = 0; i < number_of_nodes; ++i) pts.push_back(r_model_part.pGetNode(i + 1));
    Geometry<Node>::Pointer p_geometry;
    if (is_3d)
        p_geometry = Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pts));
    else
        p_geometry = Geometry<Node>::Pointer(new Quadrilateral2D4<Node>(pts));

    auto p_elem = Kratos::make_intrusive<TElement>(1, p_geometry, p_prop);
    r_model_part.AddElement(p_elem);
    rOutElement = p_elem.get();
    return r_model_part;
    KRATOS_CATCH("");
}

/// Interpolates a nodal scalar at a Gauss point with the shape functions.
double InterpolateScalarAtGP(
    const Geometry<Node>& rGeometry,
    const Variable<double>& rVariable,
    const std::size_t rPointNumber,
    const GeometryData::IntegrationMethod rMethod)
{
    const auto& r_integration_points = rGeometry.IntegrationPoints(rMethod);
    Vector N(rGeometry.PointsNumber());
    rGeometry.ShapeFunctionsValues(N, r_integration_points[rPointNumber].Coordinates());
    double value = 0.0;
    for (std::size_t i = 0; i < rGeometry.PointsNumber(); ++i) {
        value += N[i] * rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
    return value;
}

/// Applies a uniform total-strain Voigt state (shear components are already the
/// engineering shear gamma) and a uniform temperature increment.
void ApplyState(ModelPart& rModelPart, const Vector& rEpsVoigt, const double rDeltaTemperature)
{
    for (auto& n : rModelPart.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = rEpsVoigt[0] * x0[0] + 0.5 * rEpsVoigt[3] * x0[1] + 0.5 * rEpsVoigt[5] * x0[2];
        u[1] = 0.5 * rEpsVoigt[3] * x0[0] + rEpsVoigt[1] * x0[1] + 0.5 * rEpsVoigt[4] * x0[2];
        u[2] = 0.5 * rEpsVoigt[5] * x0[0] + 0.5 * rEpsVoigt[4] * x0[1] + rEpsVoigt[2] * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + rDeltaTemperature;
    }
}

/// Returns the 3D isotropic constitutive matrix for the given Young modulus.
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

/// Returns the plane-strain constitutive matrix.
Matrix CPlaneStrain(const double E, const double nu)
{
    Matrix C(3, 3);
    noalias(C) = ZeroMatrix(3, 3);
    const double c11 = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double c12 = c11 * nu / (1.0 - nu);
    const double c33 = c11 * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
    C(0,0)=c11; C(1,1)=c11; C(2,2)=c33;
    C(0,1)=c12; C(1,0)=c12;
    return C;
}

/// Returns the plane-stress constitutive matrix.
Matrix CPlaneStress(const double E, const double nu)
{
    Matrix C(3, 3);
    noalias(C) = ZeroMatrix(3, 3);
    const double c11 = E / (1.0 - nu * nu);
    const double c12 = c11 * nu;
    const double c33 = c11 * (1.0 - nu) * 0.5;
    C(0,0)=c11; C(1,1)=c11; C(2,2)=c33;
    C(0,1)=c12; C(1,0)=c12;
    return C;
}

/// The 3D (6-component) or 2D (3-component) total strain for a uniaxial-stress
/// state of the given axial strain.
Vector UniaxialTotalStrain(const bool rIs3d, const double rEps)
{
    if (rIs3d) {
        Vector e(6);
        e[0] = rEps; e[1] = -test_poisson_ratio * rEps; e[2] = -test_poisson_ratio * rEps;
        e[3] = 0.0; e[4] = 0.0; e[5] = 0.0;
        return e;
    }
    Vector e(3);
    e[0] = rEps; e[1] = -test_poisson_ratio * rEps; e[2] = 0.0;
    return e;
}

/// Analytical thermal strain base (before the plane-strain (1+nu) factor).
Vector ThermalStrainBase(const GeometryKind rGeo, const double rDeltaTemperature)
{
    const double alpha_dT = test_thermal_expansion * rDeltaTemperature;
    if (rGeo == GeometryKind::Hexa3D) {
        Vector e(6);
        e[0] = alpha_dT; e[1] = alpha_dT; e[2] = alpha_dT;
        e[3] = 0.0; e[4] = 0.0; e[5] = 0.0;
        return e;
    }
    Vector e(3);
    e[0] = alpha_dT; e[1] = alpha_dT; e[2] = 0.0;
    return e;
}

/// Reads the element-computed total strain (Voigt) at the first integration point.
Vector ElementTotalStrain(Element& rElement, const ProcessInfo& rPi)
{
    std::vector<Vector> out;
    rElement.CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, out, rPi);
    return out[0];
}

/// Reads one element Vector output at the first integration point.
Vector ElementVectorOutput(Element& rElement, const Variable<Vector>& rVariable, const ProcessInfo& rPi)
{
    std::vector<Vector> out;
    rElement.CalculateOnIntegrationPoints(rVariable, out, rPi);
    return out[0];
}

/// Reads one element Matrix output at the first integration point.
Matrix ElementMatrixOutput(Element& rElement, const Variable<Matrix>& rVariable, const ProcessInfo& rPi)
{
    std::vector<Matrix> out;
    rElement.CalculateOnIntegrationPoints(rVariable, out, rPi);
    return out[0];
}

/// The current damage variable d as exposed by the law.
double CurrentDamage(ConstitutiveLaw& rLaw)
{
    double d = 0.0;
    rLaw.GetValue(DAMAGE_VARIABLE, d);
    return d;
}

/// Verifies the three vector outputs against analytical values built from the
/// given constitutive matrix, total strain, thermal strain and damage factor.
void VerifyVectorOutputs(
    const std::string& rLabel,
    Element& rElement,
    const ProcessInfo& rPi,
    const Matrix& rC,
    const Vector& rEpsilonTh,
    const double rDamageFactor)
{
    const Vector r_epsilon_total = ElementTotalStrain(rElement, rPi);
    const Vector thermal_strain = ElementVectorOutput(rElement, THERMAL_STRAIN_VECTOR, rPi);
    const Vector thermal_stress = ElementVectorOutput(rElement, THERMAL_STRESS_VECTOR, rPi);
    const Vector mechanical_stress = ElementVectorOutput(rElement, MECHANICAL_STRESS_VECTOR, rPi);
    const Vector cauchy = ElementVectorOutput(rElement, CAUCHY_STRESS_VECTOR, rPi);

    const Vector expected_thermal_stress = rDamageFactor * prod(rC, rEpsilonTh);
    const Vector expected_mechanical_stress = rDamageFactor * prod(rC, r_epsilon_total);

    // THERMAL_STRAIN_VECTOR == epsilon_th (intentional bug fix; the legacy
    // element used to return the total strain).
    for (std::size_t i = 0; i < rEpsilonTh.size(); ++i)
        KRATOS_EXPECT_TRUE(Near(thermal_strain(i), rEpsilonTh(i)));
    // THERMAL_STRESS_VECTOR / MECHANICAL_STRESS_VECTOR damage-scaled.
    for (std::size_t i = 0; i < rEpsilonTh.size(); ++i) {
        KRATOS_EXPECT_TRUE(Near(thermal_stress(i), expected_thermal_stress(i)));
        KRATOS_EXPECT_TRUE(Near(mechanical_stress(i), expected_mechanical_stress(i)));
    }
    // Decomposition: total == mechanical - thermal.
    for (std::size_t i = 0; i < rEpsilonTh.size(); ++i)
        KRATOS_EXPECT_TRUE(Near(cauchy(i), mechanical_stress(i) - thermal_stress(i)));
    std::cout << "[5B.2] " << rLabel << ": e_th/mech/therm/decomp OK" << std::endl;
}

/// Verifies vector/tensor consistency (including the shear convention).
void VerifyTensorConsistency(
    const std::string& rLabel,
    Element& rElement,
    const ProcessInfo& rPi)
{
    const Vector thermal_strain = ElementVectorOutput(rElement, THERMAL_STRAIN_VECTOR, rPi);
    const Vector thermal_stress = ElementVectorOutput(rElement, THERMAL_STRESS_VECTOR, rPi);
    const Vector mechanical_stress = ElementVectorOutput(rElement, MECHANICAL_STRESS_VECTOR, rPi);
    const Matrix thermal_strain_tensor = ElementMatrixOutput(rElement, THERMAL_STRAIN_TENSOR, rPi);
    const Matrix thermal_stress_tensor = ElementMatrixOutput(rElement, THERMAL_STRESS_TENSOR, rPi);
    const Matrix mechanical_stress_tensor = ElementMatrixOutput(rElement, MECHANICAL_STRESS_TENSOR, rPi);

    const bool is_3d = (thermal_strain.size() == 6);
    const std::size_t dimension = is_3d ? 3 : 2;
    const std::size_t shear_index = is_3d ? 3 : 2;

    const Matrix expected_strain_tensor = MathUtils<double>::StrainVectorToTensor(thermal_strain);
    const Matrix expected_thermal_tensor = MathUtils<double>::StressVectorToTensor(thermal_stress);
    const Matrix expected_mechanical_tensor = MathUtils<double>::StressVectorToTensor(mechanical_stress);

    for (std::size_t i = 0; i < dimension; ++i) {
        for (std::size_t j = 0; j < dimension; ++j) {
            KRATOS_EXPECT_TRUE(Near(thermal_strain_tensor(i, j), expected_strain_tensor(i, j)));
            KRATOS_EXPECT_TRUE(Near(thermal_stress_tensor(i, j), expected_thermal_tensor(i, j)));
            KRATOS_EXPECT_TRUE(Near(mechanical_stress_tensor(i, j), expected_mechanical_tensor(i, j)));
        }
    }

    // Shear convention: strain tensor xy == 0.5 * Voigt shear; stress tensor
    // xy == Voigt shear.
    KRATOS_EXPECT_TRUE(Near(thermal_strain_tensor(0, 1), 0.5 * thermal_strain(shear_index)));
    KRATOS_EXPECT_TRUE(Near(thermal_stress_tensor(0, 1), thermal_stress(shear_index)));
    KRATOS_EXPECT_TRUE(Near(mechanical_stress_tensor(0, 1), mechanical_stress(shear_index)));
    std::cout << "[5B.2] " << rLabel << ": tensor consistency + shear convention OK" << std::endl;
}

} // namespace

//************************************************************************************
// 1. Has() stays false: the specialized outputs keep the parameter-aware dispatch.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(B52_SpecializedOutputs_HasRemainsFalse, KratosDamFastSuite)
{
    const std::vector<ConstitutiveLaw::Pointer> laws = {
        ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw()),
        ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal()),
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()),
        ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw()),
        ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamage3DLaw())};
    for (auto& p_law : laws) {
        for (auto* p_variable : {&THERMAL_STRAIN_VECTOR, &THERMAL_STRESS_VECTOR, &MECHANICAL_STRESS_VECTOR}) {
            KRATOS_EXPECT_FALSE(p_law->Has(*p_variable));
        }
        for (auto* p_variable : {&THERMAL_STRAIN_TENSOR, &THERMAL_STRESS_TENSOR, &MECHANICAL_STRESS_TENSOR}) {
            KRATOS_EXPECT_FALSE(p_law->Has(*p_variable));
        }
    }
    std::cout << "[5B.2] Has() false for all specialized outputs on all families" << std::endl;
}

//************************************************************************************
// 2. Nodal-linear acceptance: non-uniform E / temperature / reference temperature.
//************************************************************************************

namespace
{
void VerifyNodalLinear(
    const std::string& rLabel,
    const GeometryKind rGeo,
    ConstitutiveLaw::Pointer pLaw,
    const double rPlaneFactor)   // (1+nu) for plane strain, 1.0 otherwise
{
    Model model;
    TestSmallDisplacementElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
        model, rLabel, p_elem, pLaw, rGeo, true);
    p_elem->Initialize(r_mp.GetProcessInfo());

    // Apply the displacement state first (uniform strain).
    const bool is_3d = (rGeo == GeometryKind::Hexa3D);
    const double eps = 2.0e-6;
    ApplyState(r_mp, UniaxialTotalStrain(is_3d, eps), 0.0);

    // Non-uniform nodal state (E, temperature, reference temperature).
    for (auto& n : r_mp.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        n.FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = 2.0e7 + 1.0e6 * (x0[0] + x0[1]);
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 10.0 + 2.0 * x0[0];
        n.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature - 5.0 + x0[1];
    }

    const GeometryData::IntegrationMethod integration_method = p_elem->GetIntegrationMethod();
    const Geometry<Node>& r_geometry = p_elem->GetGeometry();
    const double e_interp = InterpolateScalarAtGP(r_geometry, NODAL_YOUNG_MODULUS, 0, integration_method);
    const double t_interp = InterpolateScalarAtGP(r_geometry, TEMPERATURE, 0, integration_method);
    const double tref_interp = InterpolateScalarAtGP(r_geometry, NODAL_REFERENCE_TEMPERATURE, 0, integration_method);
    const double delta_temperature = t_interp - tref_interp;

    const Matrix C = (rGeo == GeometryKind::Hexa3D) ? C3D(e_interp, test_poisson_ratio)
                     : (rPlaneFactor > 1.0) ? CPlaneStrain(e_interp, test_poisson_ratio)
                                            : CPlaneStress(e_interp, test_poisson_ratio);

    Vector e_th = ThermalStrainBase(rGeo, delta_temperature);
    if (rPlaneFactor > 1.0) { e_th[0] *= rPlaneFactor; e_th[1] *= rPlaneFactor; }

    const ProcessInfo& r_pi = r_mp.GetProcessInfo();
    const Vector r_epsilon_total = ElementTotalStrain(*p_elem, r_pi);
    const Vector thermal_strain = ElementVectorOutput(*p_elem, THERMAL_STRAIN_VECTOR, r_pi);
    const Vector thermal_stress = ElementVectorOutput(*p_elem, THERMAL_STRESS_VECTOR, r_pi);
    const Vector mechanical_stress = ElementVectorOutput(*p_elem, MECHANICAL_STRESS_VECTOR, r_pi);
    const Vector cauchy = ElementVectorOutput(*p_elem, CAUCHY_STRESS_VECTOR, r_pi);

    // THERMAL_STRAIN_VECTOR == epsilon_th (intentional bug fix).
    for (std::size_t i = 0; i < e_th.size(); ++i)
        KRATOS_EXPECT_TRUE(Near(thermal_strain(i), e_th(i)));
    // Stresses.
    const Vector expected_thermal_stress = prod(C, e_th);
    const Vector expected_mechanical_stress = prod(C, r_epsilon_total);
    for (std::size_t i = 0; i < e_th.size(); ++i) {
        KRATOS_EXPECT_TRUE(Near(thermal_stress(i), expected_thermal_stress(i)));
        KRATOS_EXPECT_TRUE(Near(mechanical_stress(i), expected_mechanical_stress(i)));
        KRATOS_EXPECT_TRUE(Near(cauchy(i), mechanical_stress(i) - thermal_stress(i)));
    }
    // Tensor consistency.
    VerifyTensorConsistency(rLabel, *p_elem, r_pi);
    std::cout << "[5B.2] " << rLabel << ": nodal linear acceptance OK (E_interp="
              << e_interp << ")" << std::endl;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(B52_NodalLinear_NonUniform3D, KratosDamFastSuite)
{
    VerifyNodalLinear("NodalLinear3D",
                      GeometryKind::Hexa3D,
                      ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal()),
                      1.0);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NodalLinear_NonUniformPlaneStrain, KratosDamFastSuite)
{
    VerifyNodalLinear("NodalLinearPlaneStrain",
                      GeometryKind::Quadrilateral2D,
                      ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrainNodal()),
                      1.0 + test_poisson_ratio);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NodalLinear_NonUniformPlaneStress, KratosDamFastSuite)
{
    VerifyNodalLinear("NodalLinearPlaneStress",
                      GeometryKind::Quadrilateral2D,
                      ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStressNodal()),
                      1.0);
}

//************************************************************************************
// 3. Local-damage acceptance over a full load/unload/reload sequence.
//************************************************************************************

namespace
{
template<class TElement>
void VerifyLocalDamageSequence(
    const std::string& rLabel,
    const GeometryKind rGeo,
    ConstitutiveLaw::Pointer pLaw,
    const double rPlaneFactor)
{
    Model model;
    TElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TElement>(model, rLabel, p_elem, pLaw, rGeo, false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    ProcessInfo& r_pi = r_mp.GetProcessInfo();

    const bool is_3d = (rGeo == GeometryKind::Hexa3D);
    const Matrix C = (rGeo == GeometryKind::Hexa3D) ? C3D(2.0e7, test_poisson_ratio)
                     : (rPlaneFactor > 1.0) ? CPlaneStrain(2.0e7, test_poisson_ratio)
                                            : CPlaneStress(2.0e7, test_poisson_ratio);

    // State sequences: {eps_x, dT}.
    const std::vector<std::pair<double, double>> states = {
        {0.5e-6, 40.0},   // elastic (below threshold)
        {2.0e-5, 40.0},   // damaged
        {1.0e-5, 40.0},   // unloading
        {1.5e-5, 40.0},   // reload below max
        {3.0e-5, 40.0},   // reload beyond max
        {2.0e-5, 80.0}    // combined thermo-mechanical damaged
    };

    for (std::size_t s = 0; s < states.size(); ++s) {
        const double eps = states[s].first;
        const double dT = states[s].second;
        ApplyState(r_mp, UniaxialTotalStrain(is_3d, eps), dT);
        r_pi[IS_CONVERGED] = true;
        p_elem->FinalizeSolutionStep(r_pi);

        // Committed damage factor.
        const double d = CurrentDamage(p_elem->GetConstitutiveLaw(0));
        const double damage_factor = 1.0 - d;

        Vector e_th_cur = ThermalStrainBase(rGeo, dT);
        if (rPlaneFactor > 1.0) { e_th_cur[0] *= rPlaneFactor; e_th_cur[1] *= rPlaneFactor; }

        VerifyVectorOutputs(rLabel + "_state" + std::to_string(s),
                            *p_elem, r_pi, C, e_th_cur, damage_factor);
        VerifyTensorConsistency(rLabel + "_state" + std::to_string(s), *p_elem, r_pi);

        // Committed state unchanged by the output calls.
        const double d_after = CurrentDamage(p_elem->GetConstitutiveLaw(0));
        KRATOS_EXPECT_NEAR(d_after, d, 1.0e-15);
    }
    std::cout << "[5B.2] " << rLabel << ": local-damage sequence OK" << std::endl;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(B52_LocalDamage_StateSequence3D, KratosDamFastSuite)
{
    VerifyLocalDamageSequence<TestSmallDisplacementElement>(
        "LocalSeq3D", GeometryKind::Hexa3D,
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(B52_LocalDamage_StateSequencePlaneStrain, KratosDamFastSuite)
{
    VerifyLocalDamageSequence<TestSmallDisplacementElement>(
        "LocalSeqPlaneStrain", GeometryKind::Quadrilateral2D,
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamagePlaneStrain2DLaw()), 1.0 + test_poisson_ratio);
}

KRATOS_TEST_CASE_IN_SUITE(B52_LocalDamage_StateSequencePlaneStress, KratosDamFastSuite)
{
    VerifyLocalDamageSequence<TestSmallDisplacementElement>(
        "LocalSeqPlaneStress", GeometryKind::Quadrilateral2D,
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamagePlaneStress2DLaw()), 1.0);
}

//************************************************************************************
// 4. Nonlocal-damage acceptance (Simo-Ju and Modified-Mises), 3D and 2D.
//************************************************************************************

namespace
{
template<class TElement>
void VerifyNonlocalDamage(
    const std::string& rLabel,
    const GeometryKind rGeo,
    ConstitutiveLaw::Pointer pLaw,
    const double rPlaneFactor,
    const double rNonlocalDriving)
{
    Model model;
    TElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TElement>(model, rLabel, p_elem, pLaw, rGeo, false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    ProcessInfo& r_pi = r_mp.GetProcessInfo();

    const bool is_3d = (rGeo == GeometryKind::Hexa3D);
    const Matrix C = (rGeo == GeometryKind::Hexa3D) ? C3D(2.0e7, test_poisson_ratio)
                     : (rPlaneFactor > 1.0) ? CPlaneStrain(2.0e7, test_poisson_ratio)
                                            : CPlaneStress(2.0e7, test_poisson_ratio);

    const double eps = 2.0e-5, dT = 40.0;
    ApplyState(r_mp, UniaxialTotalStrain(is_3d, eps), dT);

    // Prescribe nonlocal driving and commit.
    r_pi[IS_CONVERGED] = true;
    ConstitutiveLaw& r_law = p_elem->GetConstitutiveLaw(0);
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, rNonlocalDriving, r_pi);
    p_elem->FinalizeSolutionStep(r_pi);

    const double d = CurrentDamage(r_law);
    KRATOS_EXPECT_TRUE(d > 0.0);
    const double damage_factor = 1.0 - d;

    Vector e_th_cur = ThermalStrainBase(rGeo, dT);
    if (rPlaneFactor > 1.0) { e_th_cur[0] *= rPlaneFactor; e_th_cur[1] *= rPlaneFactor; }

    // Side-effect audit: record committed/LOCAL before outputs.
    double local_before = 0.0, threshold_before = 0.0;
    r_law.GetValue(LOCAL_EQUIVALENT_STRAIN, local_before);
    r_law.GetValue(STATE_VARIABLE, threshold_before);

    // Request all six outputs three times (deterministic).
    for (std::size_t repeat = 0; repeat < 3; ++repeat) {
        std::vector<Vector> v_out;
        std::vector<Matrix> m_out;
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, m_out, r_pi);
    }

    double local_after = 0.0, threshold_after = 0.0;
    r_law.GetValue(LOCAL_EQUIVALENT_STRAIN, local_after);
    r_law.GetValue(STATE_VARIABLE, threshold_after);
    KRATOS_EXPECT_NEAR(local_after, local_before, 1.0e-15);
    KRATOS_EXPECT_NEAR(threshold_after, threshold_before, 1.0e-15);
    KRATOS_EXPECT_NEAR(CurrentDamage(r_law), d, 1.0e-15);

    // Damage-scaled outputs + decomposition.
    VerifyVectorOutputs(rLabel, *p_elem, r_pi, C, e_th_cur, damage_factor);
    VerifyTensorConsistency(rLabel, *p_elem, r_pi);
    std::cout << "[5B.2] " << rLabel << ": d=" << d << " no local/committed side effect"
              << std::endl;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(B52_NonlocalDamage_SimoJu3D, KratosDamFastSuite)
{
    VerifyNonlocalDamage<TestSmallDisplacementElement>(
        "NonlocalSJ3D", GeometryKind::Hexa3D,
        ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw()), 1.0, 1.2e-2);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NonlocalDamage_ModifiedMises3D, KratosDamFastSuite)
{
    VerifyNonlocalDamage<TestSmallDisplacementElement>(
        "NonlocalMM3D", GeometryKind::Hexa3D,
        ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamage3DLaw()), 1.0, 1.2e-2);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NonlocalDamage_SimoJuPlaneStrain, KratosDamFastSuite)
{
    VerifyNonlocalDamage<TestSmallDisplacementElement>(
        "NonlocalSJPlaneStrain", GeometryKind::Quadrilateral2D,
        ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamagePlaneStrain2DLaw()),
        1.0 + test_poisson_ratio, 1.2e-2);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NonlocalDamage_SimoJuPlaneStress, KratosDamFastSuite)
{
    VerifyNonlocalDamage<TestSmallDisplacementElement>(
        "NonlocalSJPlaneStress", GeometryKind::Quadrilateral2D,
        ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamagePlaneStress2DLaw()),
        1.0, 1.2e-2);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NonlocalDamage_ModifiedMisesPlaneStrain, KratosDamFastSuite)
{
    VerifyNonlocalDamage<TestSmallDisplacementElement>(
        "NonlocalMMPlaneStrain", GeometryKind::Quadrilateral2D,
        ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw()),
        1.0 + test_poisson_ratio, 1.2e-2);
}

KRATOS_TEST_CASE_IN_SUITE(B52_NonlocalDamage_ModifiedMisesPlaneStress, KratosDamFastSuite)
{
    VerifyNonlocalDamage<TestSmallDisplacementElement>(
        "NonlocalMMPlaneStress", GeometryKind::Quadrilateral2D,
        ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw()),
        1.0, 1.2e-2);
}

//************************************************************************************
// 5. Side-effect regression: all six outputs requested 3x at a fixed state.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(B52_SideEffect_LocalDamage3D, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
        model, "SideEffectLocal3D", p_elem,
        ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), GeometryKind::Hexa3D, false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    ApplyState(r_mp, UniaxialTotalStrain(true, 2.0e-5), 40.0);
    r_pi[IS_CONVERGED] = true;
    p_elem->FinalizeSolutionStep(r_pi);

    ConstitutiveLaw& r_law = p_elem->GetConstitutiveLaw(0);
    const double d_before = CurrentDamage(r_law);
    double state_before = 0.0;
    r_law.GetValue(STATE_VARIABLE, state_before);

    Matrix lhs_before, lhs_after;
    Vector rhs_before, rhs_after;
    p_elem->CalculateLocalSystem(lhs_before, rhs_before, r_pi);
    const Vector total_strain_before = ElementTotalStrain(*p_elem, r_pi);
    const Vector cauchy_before = ElementVectorOutput(*p_elem, CAUCHY_STRESS_VECTOR, r_pi);
    const Vector pk2_before = ElementVectorOutput(*p_elem, PK2_STRESS_VECTOR, r_pi);

    for (std::size_t repeat = 0; repeat < 3; ++repeat) {
        std::vector<Vector> v_out;
        std::vector<Matrix> m_out;
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, m_out, r_pi);
    }

    KRATOS_EXPECT_NEAR(CurrentDamage(r_law), d_before, 1.0e-15);
    double state_after = 0.0;
    r_law.GetValue(STATE_VARIABLE, state_after);
    KRATOS_EXPECT_NEAR(state_after, state_before, 1.0e-15);
    p_elem->CalculateLocalSystem(lhs_after, rhs_after, r_pi);
    const Vector total_strain_after = ElementTotalStrain(*p_elem, r_pi);
    const Vector cauchy_after = ElementVectorOutput(*p_elem, CAUCHY_STRESS_VECTOR, r_pi);
    const Vector pk2_after = ElementVectorOutput(*p_elem, PK2_STRESS_VECTOR, r_pi);
    KRATOS_EXPECT_EQ(total_strain_after.size(), total_strain_before.size());
    for (std::size_t i = 0; i < total_strain_before.size(); ++i)
        KRATOS_EXPECT_NEAR(total_strain_after(i), total_strain_before(i), 1.0e-15);
    for (std::size_t i = 0; i < cauchy_before.size(); ++i) {
        KRATOS_EXPECT_NEAR(cauchy_after(i), cauchy_before(i), 1.0e-12);
        KRATOS_EXPECT_NEAR(pk2_after(i), pk2_before(i), 1.0e-12);
    }
    KRATOS_EXPECT_EQ(rhs_after.size(), rhs_before.size());
    for (std::size_t i = 0; i < rhs_before.size(); ++i)
        KRATOS_EXPECT_NEAR(rhs_after(i), rhs_before(i), 1.0e-12);
    std::cout << "[5B.2] SideEffectLocal3D: no side effect on committed/strain/Cauchy/PK2/LHS/RHS"
              << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(B52_SideEffect_NonlocalDamage3D, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_elem = nullptr;
    ModelPart& r_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
        model, "SideEffectNonlocal3D", p_elem,
        ConstitutiveLaw::Pointer(new TestNonlocalLaw()), GeometryKind::Hexa3D, false);
    p_elem->Initialize(r_mp.GetProcessInfo());
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    ApplyState(r_mp, UniaxialTotalStrain(true, 2.0e-5), 40.0);
    r_pi[IS_CONVERGED] = true;
    TestNonlocalLaw& r_law = dynamic_cast<TestNonlocalLaw&>(p_elem->GetConstitutiveLaw(0));
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_pi);
    p_elem->FinalizeSolutionStep(r_pi);

    const double d_before = CurrentDamage(r_law);
    const double nonlocal_before = r_law.NonlocalStrain();
    double local_before = 0.0, state_before = 0.0;
    r_law.GetValue(LOCAL_EQUIVALENT_STRAIN, local_before);
    r_law.GetValue(STATE_VARIABLE, state_before);

    Matrix lhs_before, lhs_after;
    Vector rhs_before, rhs_after;
    p_elem->CalculateLocalSystem(lhs_before, rhs_before, r_pi);
    const Vector cauchy_before = ElementVectorOutput(*p_elem, CAUCHY_STRESS_VECTOR, r_pi);

    for (std::size_t repeat = 0; repeat < 3; ++repeat) {
        std::vector<Vector> v_out;
        std::vector<Matrix> m_out;
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, m_out, r_pi);
    }

    KRATOS_EXPECT_NEAR(CurrentDamage(r_law), d_before, 1.0e-15);
    KRATOS_EXPECT_NEAR(r_law.NonlocalStrain(), nonlocal_before, 1.0e-15);
    double local_after = 0.0, state_after = 0.0;
    r_law.GetValue(LOCAL_EQUIVALENT_STRAIN, local_after);
    r_law.GetValue(STATE_VARIABLE, state_after);
    KRATOS_EXPECT_NEAR(local_after, local_before, 1.0e-15);
    KRATOS_EXPECT_NEAR(state_after, state_before, 1.0e-15);
    p_elem->CalculateLocalSystem(lhs_after, rhs_after, r_pi);
    const Vector cauchy_after = ElementVectorOutput(*p_elem, CAUCHY_STRESS_VECTOR, r_pi);
    KRATOS_EXPECT_EQ(rhs_after.size(), rhs_before.size());
    for (std::size_t i = 0; i < rhs_before.size(); ++i)
        KRATOS_EXPECT_NEAR(rhs_after(i), rhs_before(i), 1.0e-12);
    for (std::size_t i = 0; i < cauchy_before.size(); ++i)
        KRATOS_EXPECT_NEAR(cauchy_after(i), cauchy_before(i), 1.0e-12);
    std::cout << "[5B.2] SideEffectNonlocal3D: no side effect on NONLOCAL/LOCAL/committed/LHS/RHS"
              << std::endl;
}

//************************************************************************************
// 6. Shear convention with a non-zero shear state.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(B52_ShearConvention, KratosDamFastSuite)
{
    for (bool is_3d : {true, false}) {
        Model model;
        TestSmallDisplacementElement* p_elem = nullptr;
        const GeometryKind geo = is_3d ? GeometryKind::Hexa3D : GeometryKind::Quadrilateral2D;
        ConstitutiveLaw::Pointer p_law = is_3d
            ? ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw())
            : ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStressNodal());
        ModelPart& r_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
            model, is_3d ? "Shear3D" : "Shear2D", p_elem, p_law, geo, false);
        p_elem->Initialize(r_mp.GetProcessInfo());

        // Non-zero shear: e_xy = 0.5*gamma.
        const double gamma = 1.0e-5;
        Vector e_total = UniaxialTotalStrain(is_3d, 1.0e-6);
        const std::size_t shear_index = is_3d ? 3 : 2;
        e_total[shear_index] = gamma;
        ApplyState(r_mp, e_total, 25.0);

        const ProcessInfo& r_pi = r_mp.GetProcessInfo();
        const Vector thermal_strain = ElementVectorOutput(*p_elem, THERMAL_STRAIN_VECTOR, r_pi);
        const Vector thermal_stress = ElementVectorOutput(*p_elem, THERMAL_STRESS_VECTOR, r_pi);
        const Vector mechanical_stress = ElementVectorOutput(*p_elem, MECHANICAL_STRESS_VECTOR, r_pi);
        const Matrix thermal_strain_tensor = ElementMatrixOutput(*p_elem, THERMAL_STRAIN_TENSOR, r_pi);
        const Matrix thermal_stress_tensor = ElementMatrixOutput(*p_elem, THERMAL_STRESS_TENSOR, r_pi);
        const Matrix mechanical_stress_tensor = ElementMatrixOutput(*p_elem, MECHANICAL_STRESS_TENSOR, r_pi);

        // Strain tensor: e_xy = gamma_xy / 2.
        KRATOS_EXPECT_TRUE(Near(thermal_strain_tensor(0, 1), 0.5 * thermal_strain(shear_index)));
        // Stress tensor: sigma_xy = tau_xy (no 1/2).
        KRATOS_EXPECT_TRUE(Near(thermal_stress_tensor(0, 1), thermal_stress(shear_index)));
        KRATOS_EXPECT_TRUE(Near(mechanical_stress_tensor(0, 1), mechanical_stress(shear_index)));
        // Exact Vector/Tensor consistency.
        KRATOS_EXPECT_TRUE(Near(thermal_strain_tensor(0, 1),
                                0.5 * ElementVectorOutput(*p_elem, THERMAL_STRAIN_VECTOR, r_pi)(shear_index)));
        std::cout << "[5B.2] shear " << (is_3d ? "3D" : "2D")
                  << ": strain xy = gamma/2, stress xy = tau" << std::endl;
    }
}

//************************************************************************************
// 7. Dimensions: 3D vector 6 / tensor 3x3; 2D vector 3 / tensor 2x2.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(B52_Dimensions, KratosDamFastSuite)
{
    for (bool is_3d : {true, false}) {
        Model model;
        TestSmallDisplacementElement* p_elem = nullptr;
        const GeometryKind geo = is_3d ? GeometryKind::Hexa3D : GeometryKind::Quadrilateral2D;
        ConstitutiveLaw::Pointer p_law;
        if (is_3d)
            p_law = ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw());
        else
            p_law = ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStressNodal());
        ModelPart& r_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
            model, is_3d ? "Dim3D" : "Dim2D", p_elem, p_law, geo, false);
        p_elem->Initialize(r_mp.GetProcessInfo());
        ApplyState(r_mp, UniaxialTotalStrain(is_3d, 1.0e-6), 25.0);

        const ProcessInfo& r_pi = r_mp.GetProcessInfo();
        std::vector<Vector> v_out;
        std::vector<Matrix> m_out;
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, v_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, m_out, r_pi);
        p_elem->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, m_out, r_pi);

        const std::size_t expected_voigt = is_3d ? 6 : 3;
        const std::size_t expected_dim = is_3d ? 3 : 2;
        KRATOS_EXPECT_EQ(v_out[0].size(), expected_voigt);
        KRATOS_EXPECT_EQ(m_out[0].size1(), expected_dim);
        KRATOS_EXPECT_EQ(m_out[0].size2(), expected_dim);
        std::cout << "[5B.2] dimensions " << (is_3d ? "3D" : "2D")
                  << ": voigt " << expected_voigt << ", tensor " << expected_dim << "x" << expected_dim << std::endl;
    }
}

//************************************************************************************
// 8. Legacy-versus-SMA comparison over all families.
//************************************************************************************

namespace
{
void VerifyLegacyVersusSMA(
    const std::string& rLabel,
    const GeometryKind rGeo,
    ConstitutiveLaw::Pointer pLaw,
    const bool rNonlocal)
{
    ConstitutiveLaw::Pointer p_legacy_law = pLaw->Clone();
    Model model;
    TestThermoMechanicElement* p_legacy = nullptr;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_legacy_mp = CreateOutputModelPart<TestThermoMechanicElement>(
        model, rLabel + "Legacy", p_legacy, p_legacy_law, rGeo, false);
    ModelPart& r_candidate_mp = CreateOutputModelPart<TestSmallDisplacementElement>(
        model, rLabel + "SMA", p_candidate, pLaw, rGeo, false);
    p_legacy->Initialize(r_legacy_mp.GetProcessInfo());
    p_candidate->Initialize(r_candidate_mp.GetProcessInfo());

    const bool is_3d = (rGeo == GeometryKind::Hexa3D);
    const double eps = 1.5e-5, dT = 40.0;
    ApplyState(r_legacy_mp, UniaxialTotalStrain(is_3d, eps), dT);
    ApplyState(r_candidate_mp, UniaxialTotalStrain(is_3d, eps), dT);

    ProcessInfo& r_legacy_pi = r_legacy_mp.GetProcessInfo();
    ProcessInfo& r_candidate_pi = r_candidate_mp.GetProcessInfo();
    r_legacy_pi[IS_CONVERGED] = true;
    r_candidate_pi[IS_CONVERGED] = true;
    if (rNonlocal) {
        p_legacy->GetConstitutiveLaw(0).SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_legacy_pi);
        p_candidate->GetConstitutiveLaw(0).SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_candidate_pi);
    }
    p_legacy->FinalizeSolutionStep(r_legacy_pi);
    p_candidate->FinalizeSolutionStep(r_candidate_pi);

    // Stress vector outputs: candidate == legacy.
    for (auto* p_variable : {&THERMAL_STRESS_VECTOR, &MECHANICAL_STRESS_VECTOR}) {
        std::vector<Vector> legacy_out, candidate_out;
        p_legacy->CalculateOnIntegrationPoints(*p_variable, legacy_out, r_legacy_pi);
        p_candidate->CalculateOnIntegrationPoints(*p_variable, candidate_out, r_candidate_pi);
        KRATOS_EXPECT_EQ(candidate_out.size(), legacy_out.size());
        for (std::size_t i = 0; i < legacy_out.size(); ++i) {
            KRATOS_EXPECT_EQ(candidate_out[i].size(), legacy_out[i].size());
            for (std::size_t j = 0; j < legacy_out[i].size(); ++j)
                KRATOS_EXPECT_TRUE(Near(candidate_out[i](j), legacy_out[i](j)));
        }
    }
    // Stress tensor outputs: candidate == legacy.
    for (auto* p_variable : {&THERMAL_STRESS_TENSOR, &MECHANICAL_STRESS_TENSOR}) {
        std::vector<Matrix> legacy_out, candidate_out;
        p_legacy->CalculateOnIntegrationPoints(*p_variable, legacy_out, r_legacy_pi);
        p_candidate->CalculateOnIntegrationPoints(*p_variable, candidate_out, r_candidate_pi);
        KRATOS_EXPECT_EQ(candidate_out.size(), legacy_out.size());
        for (std::size_t i = 0; i < legacy_out.size(); ++i) {
            KRATOS_EXPECT_EQ(candidate_out[i].size1(), legacy_out[i].size1());
            for (std::size_t j = 0; j < legacy_out[i].size1(); ++j)
                for (std::size_t k = 0; k < legacy_out[i].size2(); ++k)
                    KRATOS_EXPECT_TRUE(Near(candidate_out[i](j, k), legacy_out[i](j, k)));
        }
    }

    // THERMAL_STRAIN: both now return the analytical thermal strain (intentional
    // bug fix on both elements; the legacy element no longer returns the total
    // strain).
    const Vector legacy_thermal_strain = ElementVectorOutput(*p_legacy, THERMAL_STRAIN_VECTOR, r_legacy_pi);
    const Vector candidate_thermal_strain = ElementVectorOutput(*p_candidate, THERMAL_STRAIN_VECTOR, r_candidate_pi);
    KRATOS_EXPECT_EQ(legacy_thermal_strain.size(), candidate_thermal_strain.size());
    for (std::size_t i = 0; i < legacy_thermal_strain.size(); ++i)
        KRATOS_EXPECT_NEAR(legacy_thermal_strain(i), candidate_thermal_strain(i), 1.0e-12);
    KRATOS_EXPECT_TRUE(Near(legacy_thermal_strain[0], test_thermal_expansion * dT));

    std::cout << "[5B.2] " << rLabel << ": legacy vs SMA match for stress outputs; "
              << "THERMAL_STRAIN == analytical on both (bug fixed)" << std::endl;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(B52_LegacyVsSMA_Linear, KratosDamFastSuite)
{
    VerifyLegacyVersusSMA("VsLinear", GeometryKind::Hexa3D,
                          ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw()), false);
}

KRATOS_TEST_CASE_IN_SUITE(B52_LegacyVsSMA_NodalLinear, KratosDamFastSuite)
{
    VerifyLegacyVersusSMA("VsNodalLinear", GeometryKind::Hexa3D,
                          ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal()), false);
}

KRATOS_TEST_CASE_IN_SUITE(B52_LegacyVsSMA_LocalDamage, KratosDamFastSuite)
{
    VerifyLegacyVersusSMA("VsLocalDamage", GeometryKind::Hexa3D,
                          ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw()), false);
}

KRATOS_TEST_CASE_IN_SUITE(B52_LegacyVsSMA_NonlocalDamage, KratosDamFastSuite)
{
    VerifyLegacyVersusSMA("VsNonlocalDamage", GeometryKind::Hexa3D,
                          ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw()), true);
}

} // namespace Testing
} // namespace Kratos
