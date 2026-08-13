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

/// Test-only element subclass exposing the constitutive-law vector (protected
/// in the Element base class) for non-invasive diagnostic access.
class TestElement : public SmallDisplacement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TestElement);
    using BaseType = SmallDisplacement;

    TestElement(IndexType NewId, GeometryType::Pointer pGeometry,
                PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties) {}

    ConstitutiveLaw& GetConstitutiveLaw(std::size_t i) { return *mConstitutiveLawVector[i]; }
};

/// Reads the committed threshold (kappa) and damage (d) of an element.
void ReadDamageState(Element& rElement, double& rThreshold, double& rDamage)
{
    const auto& r_diagnostic = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(
        static_cast<TestElement&>(rElement).GetConstitutiveLaw(0));
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

/// Compares the LHS/RHS of the historical thermo alias and the direct SMA name.
void CompareElementSystems(Element& rHistAlias, Element& rDirectSma,
                           ModelPart& rHistAliasMp, ModelPart& rDirectSmaMp)
{
    Matrix lhs_hist_alias, lhs_direct_sma;
    Vector rhs_hist_alias, rhs_direct_sma;
    rHistAlias.CalculateLocalSystem(lhs_hist_alias, rhs_hist_alias, rHistAliasMp.GetProcessInfo());
    rDirectSma.CalculateLocalSystem(lhs_direct_sma, rhs_direct_sma, rDirectSmaMp.GetProcessInfo());
    KRATOS_EXPECT_EQ(lhs_hist_alias.size1(), lhs_direct_sma.size1());
    KRATOS_EXPECT_EQ(rhs_hist_alias.size(), rhs_direct_sma.size());
    for (std::size_t i = 0; i < rhs_hist_alias.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs_direct_sma[i], rhs_hist_alias[i],
                           std::max(1.0e-9,
                                    comparison_relative_tolerance * std::abs(rhs_hist_alias[i])));
    }
    for (std::size_t i = 0; i < lhs_hist_alias.size1(); ++i) {
        for (std::size_t j = 0; j < lhs_hist_alias.size2(); ++j) {
            KRATOS_EXPECT_NEAR(lhs_direct_sma(i, j), lhs_hist_alias(i, j),
                               std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lhs_hist_alias(i, j))));
        }
    }
}

} // namespace

//************************************************************************************
// 3D lifecycle A-F: historical alias and direct SMA agree.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage3DLifecycle, KratosDamFastSuite)
{
    Model model;
    TestElement* p_hist_alias = nullptr;
    TestElement* p_direct_sma = nullptr;
    ModelPart& r_hist_alias_mp = CreateElementModelPart(
        model, "HistAlias", "SmallDisplacementThermoMechanicElement3D8N", 3, p_hist_alias);
    ModelPart& r_direct_sma_mp = CreateElementModelPart(
        model, "DirectSma", "SmallDisplacementElement3D8N", 3, p_direct_sma);

    auto& r_hist_alias_element = *p_hist_alias;
    auto& r_direct_sma_element = *p_direct_sma;
    r_hist_alias_element.Initialize(r_hist_alias_mp.GetProcessInfo());
    r_hist_alias_element.InitializeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.Initialize(r_direct_sma_mp.GetProcessInfo());
    r_direct_sma_element.InitializeSolutionStep(r_direct_sma_mp.GetProcessInfo());

    double hist_alias_threshold = 0.0, hist_alias_damage = 0.0;
    double direct_sma_threshold = 0.0, direct_sma_damage = 0.0;
    double last_hist_alias_threshold = test_damage_threshold;
    double last_hist_alias_damage = 0.0;

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
        ApplyUniaxialState(r_hist_alias_mp, eps);
        ApplyUniaxialState(r_direct_sma_mp, eps);

        // Trial LHS/RHS must agree between the historical alias and the direct SMA name.
        CompareElementSystems(r_hist_alias_element, r_direct_sma_element, r_hist_alias_mp, r_direct_sma_mp);

        // Converged finalize commits the trial state.
        r_hist_alias_mp.GetProcessInfo()[IS_CONVERGED] = true;
        r_direct_sma_mp.GetProcessInfo()[IS_CONVERGED] = true;
        r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
        r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());

        ReadDamageState(r_hist_alias_element, hist_alias_threshold, hist_alias_damage);
        ReadDamageState(r_direct_sma_element, direct_sma_threshold, direct_sma_damage);

        KRATOS_EXPECT_NEAR(direct_sma_threshold, hist_alias_threshold,
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(hist_alias_threshold)));
        KRATOS_EXPECT_NEAR(direct_sma_damage, hist_alias_damage,
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(hist_alias_damage)));

        if (step == 0) {
            // A: elastic, no damage.
            KRATOS_EXPECT_NEAR(hist_alias_threshold, test_damage_threshold, 1.0e-12);
            KRATOS_EXPECT_NEAR(hist_alias_damage, 0.0, 1.0e-12);
        } else if (step == 1) {
            // B: damage initiation.
            KRATOS_EXPECT_TRUE(hist_alias_threshold > test_damage_threshold);
            KRATOS_EXPECT_TRUE(hist_alias_damage > 0.0);
        } else if (step >= 2 && step <= 4) {
            // C: progressive damage growth.
            KRATOS_EXPECT_TRUE(hist_alias_damage >= last_hist_alias_damage);
            KRATOS_EXPECT_TRUE(hist_alias_threshold >= last_hist_alias_threshold);
        } else if (step == 5) {
            // D: unload - irreversible damage.
            KRATOS_EXPECT_NEAR(hist_alias_threshold, last_hist_alias_threshold, 1.0e-12);
            KRATOS_EXPECT_NEAR(hist_alias_damage, last_hist_alias_damage, 1.0e-12);
        } else if (step == 6) {
            // E: reload below maximum - no growth.
            KRATOS_EXPECT_NEAR(hist_alias_threshold, last_hist_alias_threshold, 1.0e-12);
            KRATOS_EXPECT_NEAR(hist_alias_damage, last_hist_alias_damage, 1.0e-12);
        } else {
            // F: reload beyond maximum - growth.
            KRATOS_EXPECT_TRUE(hist_alias_threshold > last_hist_alias_threshold);
            KRATOS_EXPECT_TRUE(hist_alias_damage >= last_hist_alias_damage);
        }
        last_hist_alias_threshold = hist_alias_threshold;
        last_hist_alias_damage = hist_alias_damage;
    }

    // Rollback: a rejected larger trial must preserve the committed state, and
    // a converged larger trial must commit a new (higher) state.
    const double committed_alias_threshold = hist_alias_threshold;
    ApplyUniaxialState(r_hist_alias_mp, 3.0e-6);
    ApplyUniaxialState(r_direct_sma_mp, 3.0e-6);
    r_hist_alias_mp.GetProcessInfo()[IS_CONVERGED] = false;
    r_direct_sma_mp.GetProcessInfo()[IS_CONVERGED] = false;
    r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());
    ReadDamageState(r_hist_alias_element, hist_alias_threshold, hist_alias_damage);
    ReadDamageState(r_direct_sma_element, direct_sma_threshold, direct_sma_damage);
    KRATOS_EXPECT_NEAR(hist_alias_threshold, committed_alias_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(direct_sma_threshold, committed_alias_threshold, 1.0e-12);

    r_hist_alias_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_direct_sma_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());
    ReadDamageState(r_hist_alias_element, hist_alias_threshold, hist_alias_damage);
    ReadDamageState(r_direct_sma_element, direct_sma_threshold, direct_sma_damage);
    KRATOS_EXPECT_TRUE(hist_alias_threshold > committed_alias_threshold);
    KRATOS_EXPECT_NEAR(direct_sma_threshold, hist_alias_threshold,
                       std::max(comparison_absolute_tolerance,
                                comparison_relative_tolerance * std::abs(hist_alias_threshold)));

    std::cout << "[damage] 3D lifecycle: alias damage=" << hist_alias_damage
              << " direct damage=" << direct_sma_damage
              << " threshold=" << hist_alias_threshold << std::endl;
}

//************************************************************************************
// Thermal coupling.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamageThermalCoupling, KratosDamFastSuite)
{
    Model model;
    TestElement* p_hist_alias = nullptr;
    TestElement* p_direct_sma = nullptr;
    ModelPart& r_hist_alias_mp = CreateElementModelPart(
        model, "ThermHistAlias", "SmallDisplacementThermoMechanicElement3D8N", 3, p_hist_alias);
    ModelPart& r_direct_sma_mp = CreateElementModelPart(
        model, "ThermDirectSma", "SmallDisplacementElement3D8N", 3, p_direct_sma);
    auto& r_hist_alias_element = *p_hist_alias;
    auto& r_direct_sma_element = *p_direct_sma;
    r_hist_alias_element.Initialize(r_hist_alias_mp.GetProcessInfo());
    r_hist_alias_element.InitializeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.Initialize(r_direct_sma_mp.GetProcessInfo());
    r_direct_sma_element.InitializeSolutionStep(r_direct_sma_mp.GetProcessInfo());

    double lt, ld, ct, cd;

    // Case 1: free thermal expansion from the pristine state (displacement =
    // alpha*dT*x, DeltaT = 50). The mechanical strain is zero, so NO artificial
    // local damage may be generated in either element.
    ApplyFreeThermalExpansion(r_hist_alias_mp, 50.0);
    ApplyFreeThermalExpansion(r_direct_sma_mp, 50.0);
    CompareElementSystems(r_hist_alias_element, r_direct_sma_element, r_hist_alias_mp, r_direct_sma_mp);
    r_hist_alias_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_direct_sma_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());
    ReadDamageState(r_hist_alias_element, lt, ld);
    ReadDamageState(r_direct_sma_element, ct, cd);
    KRATOS_EXPECT_NEAR(ld, 0.0, 1.0e-15);  // no artificial damage
    KRATOS_EXPECT_NEAR(cd, 0.0, 1.0e-15);
    KRATOS_EXPECT_NEAR(ct, lt, 1.0e-12);

    // Case 2: restrained uniform heating (no displacement, DeltaT = 50). In the
    // Simo-Ju model compressive restraint drives damage; the historical alias
    // and the direct SMA name must produce the identical damaged state.
    ApplyUniaxialState(r_hist_alias_mp, 0.0);
    ApplyUniaxialState(r_direct_sma_mp, 0.0);
    ApplyTemperatureChange(r_hist_alias_mp, 50.0);
    ApplyTemperatureChange(r_direct_sma_mp, 50.0);
    CompareElementSystems(r_hist_alias_element, r_direct_sma_element, r_hist_alias_mp, r_direct_sma_mp);
    r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());
    ReadDamageState(r_hist_alias_element, lt, ld);
    ReadDamageState(r_direct_sma_element, ct, cd);
    KRATOS_EXPECT_NEAR(ct, lt, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lt)));
    KRATOS_EXPECT_NEAR(cd, ld, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(ld)));

    // Case 3: combined mechanical + thermal loading.
    ApplyUniaxialState(r_hist_alias_mp, 2.0e-6);
    ApplyUniaxialState(r_direct_sma_mp, 2.0e-6);
    ApplyTemperatureChange(r_hist_alias_mp, 50.0);
    ApplyTemperatureChange(r_direct_sma_mp, 50.0);
    CompareElementSystems(r_hist_alias_element, r_direct_sma_element, r_hist_alias_mp, r_direct_sma_mp);
    r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());
    ReadDamageState(r_hist_alias_element, lt, ld);
    ReadDamageState(r_direct_sma_element, ct, cd);
    KRATOS_EXPECT_NEAR(ct, lt, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lt)));
    KRATOS_EXPECT_NEAR(cd, ld, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(ld)));

    // Case 4: spatially non-uniform temperature.
    ApplyUniaxialState(r_hist_alias_mp, 2.0e-6);
    ApplyUniaxialState(r_direct_sma_mp, 2.0e-6);
    ApplyNonUniformTemperature(r_hist_alias_mp);
    ApplyNonUniformTemperature(r_direct_sma_mp);
    CompareElementSystems(r_hist_alias_element, r_direct_sma_element, r_hist_alias_mp, r_direct_sma_mp);
    r_hist_alias_element.FinalizeSolutionStep(r_hist_alias_mp.GetProcessInfo());
    r_direct_sma_element.FinalizeSolutionStep(r_direct_sma_mp.GetProcessInfo());
    ReadDamageState(r_hist_alias_element, lt, ld);
    ReadDamageState(r_direct_sma_element, ct, cd);
    KRATOS_EXPECT_NEAR(ct, lt, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lt)));
    KRATOS_EXPECT_NEAR(cd, ld, std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(ld)));

    std::cout << "[damage] thermal coupling: alias damage = " << ld
              << ", direct damage = " << cd << std::endl;
}

} // namespace Testing
} // namespace Kratos
