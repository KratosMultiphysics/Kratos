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

// Phase 4D.1: constitutive-law side of the thermal NONLOCAL damage
// compatibility with StructuralMechanicsApplication::SmallDisplacement, and
// the generic StructuralMechanics mechanism for computing
// LOCAL_EQUIVALENT_STRAIN through
//   Element::CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, ...).
//
// Production scope: DamApplication only. Core, StructuralMechanics,
// Poromechanics and GeoMechanics are not modified. The nonlinear-iteration
// scheme orchestration is NOT integrated (Phase 4D.2).

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
#include "includes/stream_serializer.h"
#include "includes/variables.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"
#include "custom_utilities/dam_nonlocal_damage_utilities.hpp"

// StructuralMechanicsApplication small-displacement element
#include "custom_elements/solid_elements/small_displacement.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances (same philosophy as the previous characterization).
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;

/// Material data shared by all model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_damage_threshold = 5.0e-3;
constexpr double test_strength_ratio = 10.0;
constexpr double test_fracture_energy = 5000.0;

/// Test-only diagnostic subclass exposing the nonlocal internal state.
class DiagnosticSimoJuNonlocalDamage3DLaw : public ThermalSimoJuNonlocalDamage3DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DiagnosticSimoJuNonlocalDamage3DLaw);

    DiagnosticSimoJuNonlocalDamage3DLaw() : ThermalSimoJuNonlocalDamage3DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new DiagnosticSimoJuNonlocalDamage3DLaw(*this));
    }

    double GetLocalEquivalentStrain() const
    {
        return mpFlowRule->GetThermalVariables().DeltaPlasticDissipation;
    }

    double GetNonlocalEquivalentStrain() const
    {
        return mNonlocalEquivalentStrain;
    }

    double GetThresholdVariable() const
    {
        return mpFlowRule->GetInternalVariables().EquivalentPlasticStrain;
    }

    double GetDamageVariable() const
    {
        return mpFlowRule->GetInternalVariables().DeltaPlasticStrain;
    }
};

/// Test-only element subclass exposing the constitutive-law vector.
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

/// Creates a model part with a single SMA element of the given type and the
/// Simo-Ju nonlocal law.
template<class TLaw>
ModelPart& CreateSmaModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::size_t rDimension,
    TestSmallDisplacementElement*& rOutElement)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[IS_CONVERGED] = true;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    r_model_part.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

    const Element& r_prototype = KratosComponents<Element>::Get(rElementName);
    const auto& r_geometry = r_prototype.GetGeometry();
    Matrix local_coordinates;
    r_geometry.PointsLocalCoordinates(local_coordinates);

    const double geometry_scale = 2.5;
    array_1d<double, 3> geometry_offset;
    geometry_offset[0] = 0.75;
    geometry_offset[1] = 1.25;
    geometry_offset[2] = (rDimension == 3) ? 0.5 : 0.0;

    const std::size_t number_of_nodes = r_geometry.PointsNumber();
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const double x = geometry_scale * local_coordinates(i, 0) + geometry_offset[0];
        const double y = geometry_scale * local_coordinates(i, 1) + geometry_offset[1];
        const double z = (rDimension == 3) ? geometry_scale * local_coordinates(i, 2) + geometry_offset[2] : 0.0;
        r_model_part.CreateNewNode(i + 1, x, y, z);
    }

    Geometry<Node>::PointsArrayType points;
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        points.push_back(r_model_part.pGetNode(i + 1));
    }
    Geometry<Node>::Pointer p_geometry = r_geometry.Create(points);

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<TLaw>()->Clone());

    auto p_element = Kratos::make_intrusive<TestSmallDisplacementElement>(1, p_geometry, p_prop);
    r_model_part.AddElement(p_element);
    rOutElement = p_element.get();

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix zero_initial_stress(rDimension, rDimension);
        noalias(zero_initial_stress) = ZeroMatrix(rDimension, rDimension);
        r_node.FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = zero_initial_stress;
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Applies a non-degenerate uniaxial-STRESS field.
void ApplyUniaxialState(ModelPart& rModelPart, const double rEpsilonX)
{
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = rEpsilonX * r_x0[0];
        r_displacement[1] = -test_poisson_ratio * rEpsilonX * r_x0[1];
        r_displacement[2] = -test_poisson_ratio * rEpsilonX * r_x0[2];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
    }
}

void ApplyTemperatureChange(ModelPart& rModelPart, const double rDeltaTemperature)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + rDeltaTemperature;
    }
}

/// Applies a free-thermal-expansion field.
void ApplyFreeThermalExpansion(ModelPart& rModelPart, const double rDeltaTemperature)
{
    const double thermal_strain = test_thermal_expansion * rDeltaTemperature;
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = thermal_strain * r_x0[0];
        r_displacement[1] = thermal_strain * r_x0[1];
        r_displacement[2] = thermal_strain * r_x0[2];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
        r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + rDeltaTemperature;
    }
}

} // namespace

//************************************************************************************
// 1. SMA scalar dispatch: Has(LOCAL_EQUIVALENT_STRAIN) is false and the generic
//    CalculateOnIntegrationPoints -> CalculateOnConstitutiveLaw -> CalculateValue
//    path recomputes LOCAL from the current kinematics.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_ScalarDispatch, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "Dispatch4D1", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));
    // The dispatch contract: Has() must remain false so that the SMA element
    // takes the CalculateOnConstitutiveLaw path (not the stored GetValue path).
    KRATOS_EXPECT_FALSE(r_law.Has(LOCAL_EQUIVALENT_STRAIN));

    // Call ONLY the generic integration-point interface. No manual strain or
    // ConstitutiveLaw::Parameters is constructed.
    ApplyUniaxialState(r_mp, 2.0e-6);
    std::vector<double> local_values;
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(local_values.size(), 8);
    const double expected = std::sqrt(test_young_modulus) * 2.0e-6;
    for (double v : local_values) {
        KRATOS_EXPECT_NEAR(v, expected, 1.0e-9);
    }
    std::cout << "[4D.1] SMA dispatch: LOCAL recomputed = " << local_values[0] << std::endl;
}

//************************************************************************************
// 2. 3D SMA LOCAL acceptance (central Phase-4D.1 test).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_LocalAcceptance3D, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "Accept3D", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    ApplyUniaxialState(r_mp, 2.0e-6);

    // 1. Call ONLY CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN).
    std::vector<double> local_values;
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
    KRATOS_EXPECT_EQ(local_values.size(), 8);

    const double expected = std::sqrt(test_young_modulus) * 2.0e-6;
    for (std::size_t gp = 0; gp < local_values.size(); ++gp) {
        // 1. returned LOCAL is correct.
        KRATOS_EXPECT_NEAR(local_values[gp], expected, 1.0e-9);
        // 2. the same value is returned by law->GetValue(LOCAL_EQUIVALENT_STRAIN).
        const auto& r_diag = dynamic_cast<const DiagnosticSimoJuNonlocalDamage3DLaw&>(
            *p_element->GetConstitutiveLawPointer(gp));
        KRATOS_EXPECT_NEAR(r_diag.GetLocalEquivalentStrain(), local_values[gp], 1.0e-12);
        // 3. NONLOCAL unchanged (0, no averaging has run).
        KRATOS_EXPECT_NEAR(r_diag.GetNonlocalEquivalentStrain(), 0.0, 1.0e-15);
        // 4. damage/history not committed.
        KRATOS_EXPECT_NEAR(r_diag.GetThresholdVariable(), test_damage_threshold, 1.0e-15);
        KRATOS_EXPECT_NEAR(r_diag.GetDamageVariable(), 0.0, 1.0e-15);
    }
    std::cout << "[4D.1] 3D SMA LOCAL acceptance: " << local_values[0] << std::endl;
}

//************************************************************************************
// 3. State-change acceptance: LOCAL tracks the current displacement.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_StateChange, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "StateChange", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    ApplyUniaxialState(r_mp, 2.0e-6);
    std::vector<double> local_old;
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_old, r_mp.GetProcessInfo());

    // Change the displacement WITHOUT reinitializing the element.
    ApplyUniaxialState(r_mp, 4.0e-6);
    std::vector<double> local_new;
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_new, r_mp.GetProcessInfo());

    const double expected_old = std::sqrt(test_young_modulus) * 2.0e-6;
    const double expected_new = std::sqrt(test_young_modulus) * 4.0e-6;
    KRATOS_EXPECT_NEAR(local_old[0], expected_old, 1.0e-9);
    KRATOS_EXPECT_NEAR(local_new[0], expected_new, 1.0e-9);
    KRATOS_EXPECT_TRUE(std::abs(local_new[0] - local_old[0]) > 1.0e-12);  // LOCAL_old != LOCAL_new
    std::cout << "[4D.1] state change: LOCAL " << local_old[0] << " -> " << local_new[0] << std::endl;
}

//************************************************************************************
// 4. Thermal LOCAL acceptance through the SMA integration-point call.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_ThermalLocal, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "ThermalLocal4D1", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    std::vector<double> local_values;

    // Free thermal expansion: LOCAL must be effectively zero.
    ApplyFreeThermalExpansion(r_mp, 50.0);
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
    for (double v : local_values) {
        KRATOS_EXPECT_NEAR(v, 0.0, 1.0e-12);
    }

    // Restrained uniform heating (zero displacement, DeltaT = 50): compressive
    // mechanical strain drives a non-zero LOCAL.
    ApplyUniaxialState(r_mp, 0.0);
    ApplyTemperatureChange(r_mp, 50.0);
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
    for (double v : local_values) {
        KRATOS_EXPECT_TRUE(v > 0.0);
    }

    // Combined mechanical + thermal.
    ApplyUniaxialState(r_mp, 2.0e-6);
    ApplyTemperatureChange(r_mp, 50.0);
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
    for (double v : local_values) {
        KRATOS_EXPECT_TRUE(v > 0.0);
    }

    // Zero-temperature mechanical loading: LOCAL = sqrt(E)*eps.
    ApplyUniaxialState(r_mp, 2.0e-6);
    ApplyTemperatureChange(r_mp, 0.0);
    r_element.CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
    KRATOS_EXPECT_NEAR(local_values[0], std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
    std::cout << "[4D.1] thermal LOCAL OK" << std::endl;
}

//************************************************************************************
// 5. Legacy-vs-new LOCAL: legacy element hook vs SMA CalculateOnIntegrationPoints.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_LegacyVsNewLocal, KratosDamFastSuite)
{
    Model model;
    // Legacy model part.
    ModelPart& r_legacy_mp = model.CreateModelPart("LegacyNL", 2);
    r_legacy_mp.GetProcessInfo()[DOMAIN_SIZE] = 3;
    r_legacy_mp.GetProcessInfo()[SPACE_DIMENSION] = 3;
    r_legacy_mp.GetProcessInfo()[IS_RESTARTED] = false;
    r_legacy_mp.GetProcessInfo()[DELTA_TIME] = 1.0;
    r_legacy_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_legacy_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_legacy_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_legacy_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_legacy_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_legacy_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_legacy_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_legacy_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_legacy_mp.AddNodalSolutionStepVariable(NODAL_AREA);
    r_legacy_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
    const Element& r_proto = KratosComponents<Element>::Get("SmallDisplacementThermoMechanicElement3D8N");
    const auto& r_proto_geom = r_proto.GetGeometry();
    Matrix lc;
    r_proto_geom.PointsLocalCoordinates(lc);
    const double geometry_scale = 2.5;
    array_1d<double, 3> offset;
    offset[0] = 0.75; offset[1] = 1.25; offset[2] = 0.5;
    for (std::size_t i = 0; i < r_proto_geom.PointsNumber(); ++i) {
        r_legacy_mp.CreateNewNode(i + 1,
            geometry_scale * lc(i, 0) + offset[0],
            geometry_scale * lc(i, 1) + offset[1],
            geometry_scale * lc(i, 2) + offset[2]);
    }
    Geometry<Node>::PointsArrayType pts;
    for (std::size_t i = 0; i < r_proto_geom.PointsNumber(); ++i) pts.push_back(r_legacy_mp.pGetNode(i + 1));
    Geometry<Node>::Pointer p_geom = r_proto_geom.Create(pts);
    auto p_prop = r_legacy_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuNonlocalDamage3DLaw().Clone());
    auto p_legacy = Kratos::make_intrusive<TestThermoMechanicElement>(1, p_geom, p_prop);
    r_legacy_mp.AddElement(p_legacy);
    p_legacy->Initialize(r_legacy_mp.GetProcessInfo());
    for (auto& r_node : r_legacy_mp.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix z3(3, 3);
        noalias(z3) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }

    // Candidate SMA model part (identical geometry).
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_candidate_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "CandidateNL", "SmallDisplacementElement3D8N", 3, p_candidate);
    p_candidate->Initialize(r_candidate_mp.GetProcessInfo());

    // Apply a non-degenerate state.
    ApplyUniaxialState(r_legacy_mp, 2.0e-6);
    ApplyUniaxialState(r_candidate_mp, 2.0e-6);

    // Legacy LOCAL via the single production source (Dam LOCAL-update utility).
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_legacy_mp);
    const auto& r_legacy_diag = dynamic_cast<const DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_legacy->GetConstitutiveLaw(0));
    const double legacy_local = r_legacy_diag.GetLocalEquivalentStrain();

    // Candidate LOCAL via the generic SMA integration-point call.
    std::vector<double> candidate_local;
    p_candidate->CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, candidate_local, r_candidate_mp.GetProcessInfo());

    // GP-by-GP comparison (scale-aware for the tiny algebraic roundoff).
    for (std::size_t gp = 0; gp < candidate_local.size(); ++gp) {
        // Scale-aware: the tiny H8 element-strain roundoff (legacy vs SMA
        // kinematic paths) is amplified by sqrt(E) in the damage-driving
        // quantity; this is algebraic roundoff, not a law mismatch.
        KRATOS_EXPECT_NEAR(candidate_local[gp], legacy_local,
                           std::max(1.0e-10, 1.0e-8 * std::abs(legacy_local)));
    }
    std::cout << "[4D.1] legacy LOCAL=" << legacy_local
              << " candidate GP0=" << candidate_local[0] << std::endl;
}

//************************************************************************************
// 6. PK2/Kirchhoff/Cauchy equivalence at law level with a prescribed nonlocal
//    strain.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_MeasureEquivalence, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "MeasureEq", "SmallDisplacementElement3D8N", 3, p_element);
    p_element->Initialize(r_mp.GetProcessInfo());
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));
    r_law.InitializeMaterial(p_element->GetProperties(), p_element->GetGeometry(), Vector());

    const Vector shape_values = row(p_element->GetGeometry().ShapeFunctionsValues(
        p_element->GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        p_element->GetGeometry().ShapeFunctionsLocalGradients(p_element->GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-6;
    strain[1] = 0.0; strain[2] = 0.0; strain[3] = 0.0; strain[4] = 0.0; strain[5] = 0.0;

    r_mp.GetProcessInfo()[IS_CONVERGED] = true;

    auto run_measure = [&](ConstitutiveLaw& law, const ConstitutiveLaw::StressMeasure measure,
                           const double nonlocal, Vector& stress_out, Matrix& tangent_out) {
        Vector sw = strain;
        stress_out.resize(6, false);
        tangent_out.resize(6, 6, false);
        ConstitutiveLaw::Parameters values(
            p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
        values.SetShapeFunctionsValues(shape_values);
        values.SetShapeFunctionsDerivatives(shape_derivatives);
        values.SetDeformationGradientF(identity);
        values.SetDeterminantF(1.0);
        values.SetStrainVector(sw);
        values.SetStressVector(stress_out);
        values.SetConstitutiveMatrix(tangent_out);
        Flags& options = values.GetOptions();
        options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, nonlocal, r_mp.GetProcessInfo());
        law.CalculateMaterialResponse(values, measure);
    };

    for (double nv : {2.0e-3, 6.0e-3, 1.0e-2, 1.5e-2}) {
        Vector s_c, s_p, s_k;
        Matrix t_c, t_p, t_k;
        run_measure(r_law, ConstitutiveLaw::StressMeasure_Cauchy, nv, s_c, t_c);
        run_measure(r_law, ConstitutiveLaw::StressMeasure_PK2, nv, s_p, t_p);
        run_measure(r_law, ConstitutiveLaw::StressMeasure_Kirchhoff, nv, s_k, t_k);
        for (std::size_t i = 0; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(s_p[i], s_c[i], 1.0e-9);
            KRATOS_EXPECT_NEAR(s_k[i], s_c[i], 1.0e-9);
        }
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                KRATOS_EXPECT_NEAR(t_p(i, j), t_c(i, j), 1.0e-6);
                KRATOS_EXPECT_NEAR(t_k(i, j), t_c(i, j), 1.0e-6);
            }
        }
        // Commit through the Cauchy finalize, then repeat with the other
        // measures to check the unload/reload behaviour equivalence.
        Vector sw = strain;
        Vector stress(6);
        Matrix tangent(6, 6);
        ConstitutiveLaw::Parameters values(
            p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
        values.SetShapeFunctionsValues(shape_values);
        values.SetShapeFunctionsDerivatives(shape_derivatives);
        values.SetDeformationGradientF(identity);
        values.SetDeterminantF(1.0);
        values.SetStrainVector(sw);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(tangent);
        Flags& options = values.GetOptions();
        options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_law.FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
    }
    std::cout << "[4D.1] measure equivalence OK" << std::endl;
}

//************************************************************************************
// 7. Lifecycle: SMA InitializeSolutionStep / FinalizeSolutionStep run and the
//    finalize honours IS_CONVERGED through all three measures.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_LifecycleRollback, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "Life4D1", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());
    r_element.InitializeSolutionStep(r_mp.GetProcessInfo());
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));

    const Vector shape_values = row(p_element->GetGeometry().ShapeFunctionsValues(
        p_element->GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        p_element->GetGeometry().ShapeFunctionsLocalGradients(p_element->GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-5;
    Vector stress(6);
    Matrix tangent(6, 6);
    Vector sw = strain;
    ConstitutiveLaw::Parameters values(
        p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(sw);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    Flags& options = values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    // Committed nonlocal damaged state.
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    r_law.FinalizeMaterialResponseCauchy(values);
    const double committed_threshold = r_law.GetThresholdVariable();

    // Rejected PK2 finalize with a larger nonlocal: restore.
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 2.0e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponsePK2(values);
    r_mp.GetProcessInfo()[IS_CONVERGED] = false;
    r_law.FinalizeMaterialResponsePK2(values);
    KRATOS_EXPECT_NEAR(r_law.GetThresholdVariable(), committed_threshold, 1.0e-12);

    // Converged Kirchhoff finalize: commit.
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_law.FinalizeMaterialResponseKirchhoff(values);
    KRATOS_EXPECT_TRUE(r_law.GetThresholdVariable() > committed_threshold);

    // SMA FinalizeSolutionStep runs (stateful, RequiresFinalizeMaterialResponse).
    r_element.FinalizeSolutionStep(r_mp.GetProcessInfo());
    std::cout << "[4D.1] lifecycle/rollback OK, committed=" << committed_threshold << std::endl;
}

//************************************************************************************
// 8. Input-strain non-mutation: normal responses and LOCAL CalculateValue leave
//    the total strain vector bit-identical.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_InputStrainNotMutated, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "NonMut4D1", "SmallDisplacementElement3D8N", 3, p_element);
    p_element->Initialize(r_mp.GetProcessInfo());
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));
    r_law.InitializeMaterial(p_element->GetProperties(), p_element->GetGeometry(), Vector());
    ApplyTemperatureChange(r_mp, 30.0);

    Vector total_strain(6);
    total_strain[0] = 2.0e-5; total_strain[1] = 1.0e-6; total_strain[2] = 0.0;
    total_strain[3] = 3.0e-6; total_strain[4] = 0.0; total_strain[5] = 0.0;
    const Vector original = total_strain;

    const Vector shape_values = row(p_element->GetGeometry().ShapeFunctionsValues(
        p_element->GetIntegrationMethod()), 0);
    Matrix identity = IdentityMatrix(3, 3);

    auto run = [&](const std::string& label, auto fn) {
        Vector sw = total_strain;
        Vector stress(6);
        Matrix tangent(6, 6);
        ConstitutiveLaw::Parameters values(
            p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
        values.SetShapeFunctionsValues(shape_values);
        values.SetDeformationGradientF(identity);
        values.SetDeterminantF(1.0);
        values.SetStrainVector(sw);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(tangent);
        Flags& options = values.GetOptions();
        options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        fn(r_law, values);
        for (std::size_t i = 0; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(sw[i], total_strain[i], 1.0e-15) << label;
        }
    };

    run("Cauchy", [](ConstitutiveLaw& law, ConstitutiveLaw::Parameters& v) { law.CalculateMaterialResponseCauchy(v); });
    run("PK2", [](ConstitutiveLaw& law, ConstitutiveLaw::Parameters& v) { law.CalculateMaterialResponsePK2(v); });
    run("LocalCalcValue", [](ConstitutiveLaw& law, ConstitutiveLaw::Parameters& v) {
        double value = 0.0;
        law.CalculateValue(v, LOCAL_EQUIVALENT_STRAIN, value);
    });
    for (std::size_t i = 0; i < 6; ++i) {
        KRATOS_EXPECT_NEAR(total_strain[i], original[i], 1.0e-15);
    }
}

//************************************************************************************
// 9. 2D inherited family: SMA 2D elements + Simo-Ju and Modified-Mises laws.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_2DFamily, KratosDamFastSuite)
{
    struct Case { const char* law_name; const char* element; bool modified_mises; bool plane_stress; };
    const Case cases[] = {
        {"SimoJu_PS", "SmallDisplacementElement2D4N", false, false},
        {"SimoJu_PStress", "SmallDisplacementElement2D4N", false, true},
        {"ModMises_PS", "SmallDisplacementElement2D4N", true, false},
        {"ModMises_PStress", "SmallDisplacementElement2D4N", true, true},
    };

    for (const auto& c : cases) {
        Model model;
        TestSmallDisplacementElement* p_element = nullptr;
        // Always use the correct 2D nonlocal law for the element.
        ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
            model, std::string("NL2D_") + c.law_name, c.element, 2, p_element);
        ConstitutiveLaw::Pointer p_law;
        if (c.modified_mises) {
            p_law = c.plane_stress
                ? ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw())
                : ConstitutiveLaw::Pointer(new ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw());
        } else if (c.plane_stress) {
            p_law = ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamagePlaneStress2DLaw());
        } else {
            p_law = ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamagePlaneStrain2DLaw());
        }
        p_element->GetProperties().SetValue(CONSTITUTIVE_LAW, p_law->Clone());
        p_element->Initialize(r_mp.GetProcessInfo());

        ApplyUniaxialState(r_mp, 2.0e-6);
        std::vector<double> local_values;
        p_element->CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());
        KRATOS_EXPECT_EQ(local_values.size(), 4);
        for (double v : local_values) {
            KRATOS_EXPECT_TRUE(v > 0.0);
        }
        std::cout << "[4D.1] 2D " << c.law_name << " LOCAL=" << local_values[0] << std::endl;
    }
}

//************************************************************************************
// 10. Serialization / restart at a damaged nonlocal state.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_SerializationRestart, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "SerRestart", "SmallDisplacementElement3D8N", 3, p_element);
    p_element->Initialize(r_mp.GetProcessInfo());

    auto p_law = Kratos::make_shared<ThermalSimoJuNonlocalDamage3DLaw>();
    p_law->InitializeMaterial(p_element->GetProperties(), p_element->GetGeometry(), Vector());

    const Vector shape_values = row(p_element->GetGeometry().ShapeFunctionsValues(
        p_element->GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        p_element->GetGeometry().ShapeFunctionsLocalGradients(p_element->GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-5;
    Vector stress(6);
    Matrix tangent(6, 6);
    Vector sw = strain;
    ConstitutiveLaw::Parameters values(
        p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(sw);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    Flags& options = values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;

    // Damaged nonlocal state with a stored LOCAL quantity (via CalculateValue).
    double local_value = 0.0;
    p_law->CalculateValue(values, LOCAL_EQUIVALENT_STRAIN, local_value);
    p_law->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_mp.GetProcessInfo());
    p_law->CalculateMaterialResponseCauchy(values);
    p_law->FinalizeMaterialResponseCauchy(values);

    double committed_threshold = 0.0;
    p_law->GetValue(STATE_VARIABLE, committed_threshold);
    KRATOS_EXPECT_TRUE(committed_threshold > test_damage_threshold);

    // Reference stress before restart.
    Vector stress_pre(6);
    Matrix tangent_pre(6, 6);
    Vector sw2 = strain;
    ConstitutiveLaw::Parameters values_pre(
        p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
    values_pre.SetShapeFunctionsValues(shape_values);
    values_pre.SetShapeFunctionsDerivatives(shape_derivatives);
    values_pre.SetDeformationGradientF(identity);
    values_pre.SetDeterminantF(1.0);
    values_pre.SetStrainVector(sw2);
    values_pre.SetStressVector(stress_pre);
    values_pre.SetConstitutiveMatrix(tangent_pre);
    Flags& opts_pre = values_pre.GetOptions();
    opts_pre.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    opts_pre.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    opts_pre.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    p_law->CalculateMaterialResponseCauchy(values_pre);
    const Vector reference_stress = stress_pre;

    // Serialize / deserialize.
    ConstitutiveLaw::Pointer p_to_serialize = p_law;
    StreamSerializer serializer;
    serializer.save("ThermalSimoJuNonlocalDamage3DLaw", p_to_serialize);
    ConstitutiveLaw::Pointer p_loaded;
    serializer.load("ThermalSimoJuNonlocalDamage3DLaw", p_loaded);
    KRATOS_EXPECT_TRUE(p_loaded != nullptr);

    // Restart continuation: re-establish the transient material Properties on
    // the restored law (not serialized) without touching the committed state.
    auto* p_loaded_dam = dynamic_cast<ThermalNonlocalDamage3DLaw*>(p_loaded.get());
    KRATOS_EXPECT_TRUE(p_loaded_dam != nullptr);
    p_loaded_dam->ReinitializeMaterialProperties(p_element->GetProperties());

    double loaded_threshold = 0.0;
    p_loaded->GetValue(STATE_VARIABLE, loaded_threshold);
    KRATOS_EXPECT_NEAR(loaded_threshold, committed_threshold, 1.0e-12);

    // The restored NONLOCAL strain is preserved: the same response reproduces
    // the same stress/damage as before restart.
    Vector stress_after(6);
    Matrix tangent_after(6, 6);
    Vector sw3 = strain;
    ConstitutiveLaw::Parameters values_after(
        p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
    values_after.SetShapeFunctionsValues(shape_values);
    values_after.SetShapeFunctionsDerivatives(shape_derivatives);
    values_after.SetDeformationGradientF(identity);
    values_after.SetDeterminantF(1.0);
    values_after.SetStrainVector(sw3);
    values_after.SetStressVector(stress_after);
    values_after.SetConstitutiveMatrix(tangent_after);
    Flags& opts_after = values_after.GetOptions();
    opts_after.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    opts_after.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    opts_after.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    p_loaded->CalculateMaterialResponseCauchy(values_after);
    for (std::size_t i = 0; i < 6; ++i) {
        KRATOS_EXPECT_NEAR(stress_after[i], reference_stress[i], 1.0e-9);
    }
    std::cout << "[4D.1] serialization: threshold " << committed_threshold
              << " -> " << loaded_threshold << ", stress reproduced" << std::endl;
}

//************************************************************************************
// 11. Clone: preserves threshold, LOCAL and NONLOCAL; independent evolution.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamageSMA_CloneState, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateSmaModelPart<DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "Clone4D1", "SmallDisplacementElement3D8N", 3, p_element);
    p_element->Initialize(r_mp.GetProcessInfo());
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));

    ApplyUniaxialState(r_mp, 2.0e-6);
    std::vector<double> local_values;
    p_element->CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, local_values, r_mp.GetProcessInfo());

    const Vector shape_values = row(p_element->GetGeometry().ShapeFunctionsValues(
        p_element->GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        p_element->GetGeometry().ShapeFunctionsLocalGradients(p_element->GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-5;
    Vector stress(6);
    Matrix tangent(6, 6);
    Vector sw = strain;
    ConstitutiveLaw::Parameters values(
        p_element->GetGeometry(), p_element->GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(sw);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    Flags& options = values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    r_law.FinalizeMaterialResponseCauchy(values);

    const double committed_threshold = r_law.GetThresholdVariable();
    const double committed_nonlocal = r_law.GetNonlocalEquivalentStrain();
    const double stored_local = r_law.GetLocalEquivalentStrain();

    auto p_clone = Kratos::make_shared<DiagnosticSimoJuNonlocalDamage3DLaw>(r_law);
    KRATOS_EXPECT_NEAR(p_clone->GetThresholdVariable(), committed_threshold, 1.0e-15);
    KRATOS_EXPECT_NEAR(p_clone->GetNonlocalEquivalentStrain(), committed_nonlocal, 1.0e-15);
    KRATOS_EXPECT_NEAR(p_clone->GetLocalEquivalentStrain(), stored_local, 1.0e-15);

    // Independent evolution of the original.
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 2.0e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    r_law.FinalizeMaterialResponseCauchy(values);
    KRATOS_EXPECT_NEAR(p_clone->GetThresholdVariable(), committed_threshold, 1.0e-15);
    std::cout << "[4D.1] clone: original threshold=" << r_law.GetThresholdVariable()
              << ", clone threshold=" << p_clone->GetThresholdVariable() << std::endl;
}

} // namespace Testing
} // namespace Kratos
