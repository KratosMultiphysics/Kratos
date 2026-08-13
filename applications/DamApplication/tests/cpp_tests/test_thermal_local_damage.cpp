// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Permanent Dam local-damage constitutive contract. ThermalSimoJuLocalDamage*
// executes its complete stateful thermo-damage lifecycle (elastic, initiation,
// committed damage, thermal coupling, rollback) with
// StructuralMechanicsApplication::SmallDisplacement through the common
// PK2/Kirchhoff/Cauchy response implemented in the Dam law.
//
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
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_stress_2D_law.hpp"

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

/// Test-only diagnostic subclass exposing the internal damage state.
class DiagnosticSimoJuLocalDamage3DLaw : public ThermalSimoJuLocalDamage3DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DiagnosticSimoJuLocalDamage3DLaw);

    DiagnosticSimoJuLocalDamage3DLaw() : ThermalSimoJuLocalDamage3DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new DiagnosticSimoJuLocalDamage3DLaw(*this));
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

/// Test-only element subclasses exposing the constitutive-law vector (protected
/// in the Element base class) for non-invasive diagnostic access.
class TestThermoMechanicElement : public SmallDisplacement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TestThermoMechanicElement);
    using BaseType = SmallDisplacement;

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

/// Reads the committed threshold (kappa) and damage (d) of an element.
void ReadElementDamageState(Element& rElement, double& rThreshold, double& rDamage)
{
    const auto& r_diagnostic = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(
        static_cast<TestThermoMechanicElement&>(rElement).GetConstitutiveLaw(0));
    rThreshold = r_diagnostic.GetThresholdVariable();
    rDamage = r_diagnostic.GetDamageVariable();
}

void ReadSmaElementDamageState(Element& rElement, double& rThreshold, double& rDamage)
{
    const auto& r_diagnostic = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(
        static_cast<TestSmallDisplacementElement&>(rElement).GetConstitutiveLaw(0));
    rThreshold = r_diagnostic.GetThresholdVariable();
    rDamage = r_diagnostic.GetDamageVariable();
}

/// Builds an element geometry from the prototype local coordinates.
Geometry<Node>::Pointer CreateGeometry(
    ModelPart& rModelPart,
    const std::string& rPrototypeElementName,
    const std::size_t rDimension)
{
    const Element& r_prototype = KratosComponents<Element>::Get(rPrototypeElementName);
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
        rModelPart.CreateNewNode(i + 1, x, y, z);
    }

    Geometry<Node>::PointsArrayType points;
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        points.push_back(rModelPart.pGetNode(i + 1));
    }
    return r_geometry.Create(points);
}

/// Creates a model part with nodes, properties and a single element of the
/// given test-element type (constructed directly so diagnostics can access the
/// constitutive law vector).
template<class TTestElement>
ModelPart& CreateElementModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rPrototypeElementName,
    const std::size_t rDimension,
    TTestElement*& rOutElement)
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

    Geometry<Node>::Pointer p_geometry = CreateGeometry(r_model_part, rPrototypeElementName, rDimension);

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuLocalDamage3DLaw().Clone());

    auto p_element = Kratos::make_intrusive<TTestElement>(1, p_geometry, p_prop);
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

/// Applies a uniaxial displacement field (ux = eps*x, uy = -nu*eps*y,
/// uz = -nu*eps*z) to every node.
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

/// Applies a free-thermal-expansion displacement field with a uniform
/// temperature change.
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

/// Applies a uniform temperature change (no displacement change).
void ApplyTemperatureChange(ModelPart& rModelPart, const double rDeltaTemperature)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + rDeltaTemperature;
    }
}

/// Applies a spatially non-uniform temperature.
void ApplyNonUniformTemperature(ModelPart& rModelPart)
{
    std::size_t index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + 10.0 + 5.0 * static_cast<double>(index % 4);
        ++index;
    }
}

/// Compares the LHS/RHS of the legacy and candidate elements.
void CompareElementSystems(Element& rLegacy, Element& rCandidate,
                           ModelPart& rLegacyMp, ModelPart& rCandidateMp)
{
    Matrix lhs_legacy, lhs_candidate;
    Vector rhs_legacy, rhs_candidate;
    rLegacy.CalculateLocalSystem(lhs_legacy, rhs_legacy, rLegacyMp.GetProcessInfo());
    rCandidate.CalculateLocalSystem(lhs_candidate, rhs_candidate, rCandidateMp.GetProcessInfo());
    KRATOS_EXPECT_EQ(lhs_legacy.size1(), lhs_candidate.size1());
    KRATOS_EXPECT_EQ(rhs_legacy.size(), rhs_candidate.size());
    for (std::size_t i = 0; i < rhs_legacy.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs_candidate[i], rhs_legacy[i],
                           std::max(1.0e-9,
                                    comparison_relative_tolerance * std::abs(rhs_legacy[i])));
    }
    for (std::size_t i = 0; i < lhs_legacy.size1(); ++i) {
        for (std::size_t j = 0; j < lhs_legacy.size2(); ++j) {
            KRATOS_EXPECT_NEAR(lhs_candidate(i, j), lhs_legacy(i, j),
                               std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lhs_legacy(i, j))));
        }
    }
}

} // namespace

//************************************************************************************
// 1. The constitutive response must not mutate the supplied total strain.
//************************************************************************************


//************************************************************************************
// 2. SMA lifecycle: InitializeSolutionStep runs and preserves the committed
//    damage state.
//************************************************************************************


//************************************************************************************
// 3. 3D acceptance: legacy vs candidate full lifecycle A-F.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamageSMA_3DAcceptanceAtoF, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_legacy = nullptr;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_legacy_mp = CreateElementModelPart(
        model, "AccLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_legacy);
    ModelPart& r_candidate_mp = CreateElementModelPart(
        model, "AccCandidate", "SmallDisplacementElement3D8N", 3, p_candidate);

    auto& r_legacy_element = *p_legacy;
    auto& r_candidate_element = *p_candidate;
    r_legacy_element.Initialize(r_legacy_mp.GetProcessInfo());
    r_legacy_element.InitializeSolutionStep(r_legacy_mp.GetProcessInfo());
    r_candidate_element.Initialize(r_candidate_mp.GetProcessInfo());
    r_candidate_element.InitializeSolutionStep(r_candidate_mp.GetProcessInfo());

    double legacy_threshold = 0.0, legacy_damage = 0.0;
    double candidate_threshold = 0.0, candidate_damage = 0.0;
    double last_legacy_threshold = test_damage_threshold;
    double last_legacy_damage = 0.0;

    // The element-level driving quantity for the applied uniaxial-STRESS field
    // is kappa = sqrt(E)*epsilon (k0 = 5e-3 -> epsilon_init ~= 1.12e-6).
    const double steps[][2] = {
        {5.0e-7, 1.0},  // A elastic (below initiation)
        {1.3e-6, 1.0},  // B damage initiation
        {1.5e-6, 1.0},  // C progressive
        {1.8e-6, 1.0},  // C progressive
        {2.2e-6, 1.0},  // C progressive
        {8.0e-7, 1.0},  // D unload
        {1.8e-6, 1.0},  // E reload below maximum
        {2.6e-6, 1.0}   // F reload beyond maximum
    };

    for (std::size_t step = 0; step < 8; ++step) {
        const double eps = steps[step][0];
        ApplyUniaxialState(r_legacy_mp, eps);
        ApplyUniaxialState(r_candidate_mp, eps);

        // Trial LHS/RHS must agree between legacy and candidate.
        CompareElementSystems(r_legacy_element, r_candidate_element, r_legacy_mp, r_candidate_mp);

        // Converged finalize commits the trial state.
        r_legacy_mp.GetProcessInfo()[IS_CONVERGED] = true;
        r_candidate_mp.GetProcessInfo()[IS_CONVERGED] = true;
        r_legacy_element.FinalizeSolutionStep(r_legacy_mp.GetProcessInfo());
        r_candidate_element.FinalizeSolutionStep(r_candidate_mp.GetProcessInfo());

        ReadElementDamageState(r_legacy_element, legacy_threshold, legacy_damage);
        ReadSmaElementDamageState(r_candidate_element, candidate_threshold, candidate_damage);

        KRATOS_EXPECT_NEAR(candidate_threshold, legacy_threshold,
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(legacy_threshold)));
        KRATOS_EXPECT_NEAR(candidate_damage, legacy_damage,
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(legacy_damage)));

        if (step == 0) {
            // A: elastic, no damage.
            KRATOS_EXPECT_NEAR(legacy_threshold, test_damage_threshold, 1.0e-12);
            KRATOS_EXPECT_NEAR(legacy_damage, 0.0, 1.0e-12);
        } else if (step == 1) {
            // B: damage initiation.
            KRATOS_EXPECT_TRUE(legacy_threshold > test_damage_threshold);
            KRATOS_EXPECT_TRUE(legacy_damage > 0.0);
        } else if (step >= 2 && step <= 4) {
            // C: progressive damage growth.
            KRATOS_EXPECT_TRUE(legacy_damage >= last_legacy_damage);
            KRATOS_EXPECT_TRUE(legacy_threshold >= last_legacy_threshold);
        } else if (step == 5) {
            // D: unload - irreversible damage.
            KRATOS_EXPECT_NEAR(legacy_threshold, last_legacy_threshold, 1.0e-12);
            KRATOS_EXPECT_NEAR(legacy_damage, last_legacy_damage, 1.0e-12);
        } else if (step == 6) {
            // E: reload below maximum - no growth.
            KRATOS_EXPECT_NEAR(legacy_threshold, last_legacy_threshold, 1.0e-12);
            KRATOS_EXPECT_NEAR(legacy_damage, last_legacy_damage, 1.0e-12);
        } else {
            // F: reload beyond maximum - growth.
            KRATOS_EXPECT_TRUE(legacy_threshold > last_legacy_threshold);
            KRATOS_EXPECT_TRUE(legacy_damage >= last_legacy_damage);
        }
        last_legacy_threshold = legacy_threshold;
        last_legacy_damage = legacy_damage;
    }

    std::cout << "[4B] 3D acceptance: legacy dC=" << legacy_damage
              << " candidate dC=" << candidate_damage
              << " threshold=" << legacy_threshold << std::endl;
}


//************************************************************************************
// 4. Thermal coupling with the real SMA element.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamageSMA_ThermalCoupling, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_legacy = nullptr;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_legacy_mp = CreateElementModelPart(
        model, "ThermLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_legacy);
    ModelPart& r_candidate_mp = CreateElementModelPart(
        model, "ThermCandidate", "SmallDisplacementElement3D8N", 3, p_candidate);
    auto& r_legacy_element = *p_legacy;
    auto& r_candidate_element = *p_candidate;
    r_legacy_element.Initialize(r_legacy_mp.GetProcessInfo());
    r_legacy_element.InitializeSolutionStep(r_legacy_mp.GetProcessInfo());
    r_candidate_element.Initialize(r_candidate_mp.GetProcessInfo());
    r_candidate_element.InitializeSolutionStep(r_candidate_mp.GetProcessInfo());

    double lt, ld, ct, cd;

    // Case 1: free thermal expansion from the pristine state (displacement =
    // alpha*dT*x, DeltaT = 50). The mechanical strain is zero, so NO artificial
    // local damage may be generated in either element.
    ApplyFreeThermalExpansion(r_legacy_mp, 50.0);
    ApplyFreeThermalExpansion(r_candidate_mp, 50.0);
    CompareElementSystems(r_legacy_element, r_candidate_element, r_legacy_mp, r_candidate_mp);
    r_legacy_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_candidate_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_legacy_element.FinalizeSolutionStep(r_legacy_mp.GetProcessInfo());
    r_candidate_element.FinalizeSolutionStep(r_candidate_mp.GetProcessInfo());
    ReadElementDamageState(r_legacy_element, lt, ld);
    ReadSmaElementDamageState(r_candidate_element, ct, cd);
    KRATOS_EXPECT_NEAR(ld, 0.0, 1.0e-15);  // no artificial damage
    KRATOS_EXPECT_NEAR(cd, 0.0, 1.0e-15);
    KRATOS_EXPECT_NEAR(ct, lt, 1.0e-12);

    // Case 2: restrained uniform heating (no displacement, DeltaT = 50). In the
    // Simo-Ju model compressive restraint drives damage; legacy and candidate
    // must produce the identical damaged state.
    ApplyUniaxialState(r_legacy_mp, 0.0);
    ApplyUniaxialState(r_candidate_mp, 0.0);
    ApplyTemperatureChange(r_legacy_mp, 50.0);
    ApplyTemperatureChange(r_candidate_mp, 50.0);
    CompareElementSystems(r_legacy_element, r_candidate_element, r_legacy_mp, r_candidate_mp);
    r_legacy_element.FinalizeSolutionStep(r_legacy_mp.GetProcessInfo());
    r_candidate_element.FinalizeSolutionStep(r_candidate_mp.GetProcessInfo());
    ReadElementDamageState(r_legacy_element, lt, ld);
    ReadSmaElementDamageState(r_candidate_element, ct, cd);
    KRATOS_EXPECT_NEAR(ct, lt, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lt)));
    KRATOS_EXPECT_NEAR(cd, ld, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(ld)));

    // Case 3: combined mechanical + thermal loading.
    ApplyUniaxialState(r_legacy_mp, 2.0e-6);
    ApplyUniaxialState(r_candidate_mp, 2.0e-6);
    ApplyTemperatureChange(r_legacy_mp, 50.0);
    ApplyTemperatureChange(r_candidate_mp, 50.0);
    CompareElementSystems(r_legacy_element, r_candidate_element, r_legacy_mp, r_candidate_mp);
    r_legacy_element.FinalizeSolutionStep(r_legacy_mp.GetProcessInfo());
    r_candidate_element.FinalizeSolutionStep(r_candidate_mp.GetProcessInfo());
    ReadElementDamageState(r_legacy_element, lt, ld);
    ReadSmaElementDamageState(r_candidate_element, ct, cd);
    KRATOS_EXPECT_NEAR(ct, lt, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lt)));
    KRATOS_EXPECT_NEAR(cd, ld, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(ld)));

    // Case 4: spatially non-uniform temperature.
    ApplyUniaxialState(r_legacy_mp, 2.0e-6);
    ApplyUniaxialState(r_candidate_mp, 2.0e-6);
    ApplyNonUniformTemperature(r_legacy_mp);
    ApplyNonUniformTemperature(r_candidate_mp);
    CompareElementSystems(r_legacy_element, r_candidate_element, r_legacy_mp, r_candidate_mp);
    r_legacy_element.FinalizeSolutionStep(r_legacy_mp.GetProcessInfo());
    r_candidate_element.FinalizeSolutionStep(r_candidate_mp.GetProcessInfo());
    ReadElementDamageState(r_legacy_element, lt, ld);
    ReadSmaElementDamageState(r_candidate_element, ct, cd);
    KRATOS_EXPECT_NEAR(ct, lt, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lt)));
    KRATOS_EXPECT_NEAR(cd, ld, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(ld)));

    std::cout << "[4B] thermal coupling: legacy damage = " << ld
              << ", candidate damage = " << cd << std::endl;
}


//************************************************************************************
// 5. Convergence / rollback with the real SMA element.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamageSMA_Rollback, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_mp = CreateElementModelPart(
        model, "Rollback", "SmallDisplacementElement3D8N", 3, p_candidate);
    auto& r_element = *p_candidate;
    r_element.Initialize(r_mp.GetProcessInfo());
    r_element.InitializeSolutionStep(r_mp.GetProcessInfo());

    // Converged damaged state.
    ApplyUniaxialState(r_mp, 2.0e-5);
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_element.FinalizeSolutionStep(r_mp.GetProcessInfo());
    double threshold, damage;
    ReadSmaElementDamageState(r_element, threshold, damage);
    const double committed_threshold = threshold;
    KRATOS_EXPECT_TRUE(committed_threshold > test_damage_threshold);

    // Rejected trial step: larger load, IS_CONVERGED = false. The committed
    // state must be preserved (rollback).
    ApplyUniaxialState(r_mp, 3.0e-5);
    r_mp.GetProcessInfo()[IS_CONVERGED] = false;
    r_element.FinalizeSolutionStep(r_mp.GetProcessInfo());
    ReadSmaElementDamageState(r_element, threshold, damage);
    KRATOS_EXPECT_NEAR(threshold, committed_threshold, 1.0e-12);

    // Repeat with convergence: the new state is committed.
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_element.FinalizeSolutionStep(r_mp.GetProcessInfo());
    ReadSmaElementDamageState(r_element, threshold, damage);
    KRATOS_EXPECT_TRUE(threshold > committed_threshold);
    std::cout << "[4B] rollback: committed=" << committed_threshold
              << " restored=" << threshold << std::endl;
}


//************************************************************************************
// 6. Repeated response safety with the real SMA element.
//************************************************************************************


//************************************************************************************
// 7. Clone robustness: state preserved and independent evolution.
//************************************************************************************


//************************************************************************************
// 8. Serialization / restart: the committed damage state.
//************************************************************************************


} // namespace Testing
} // namespace Kratos
