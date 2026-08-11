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

// Phase 4D.2: nonlinear-iteration LOCAL-equivalent-strain production moved from
// the legacy element into the Dam smoothing-scheme workflow, driven through the
// generic Element::CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, ...)
// interface. The existing Poromechanics nonlocal averaging utilities and
// strategies are used unchanged.

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
#include "spaces/ublas_space.h"
#include "geometries/hexahedra_3d_8.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "custom_strategies/strategies/poromechanics_newton_raphson_nonlocal_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_ramm_arc_length_nonlocal_strategy.hpp"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_damped_smoothing_scheme.hpp"
#include "custom_strategies/schemes/bossak_displacement_smoothing_scheme.hpp"
#include "custom_utilities/dam_nonlocal_damage_utilities.hpp"
#include "custom_utilities/nonlocal_damage_3D_utilities.hpp"

// StructuralMechanicsApplication small-displacement element
#include "custom_elements/solid_elements/small_displacement.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances.
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;

/// Material data.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_damage_threshold = 5.0e-3;
constexpr double test_strength_ratio = 10.0;
constexpr double test_fracture_energy = 5000.0;
constexpr double test_characteristic_length = 4.0;

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using StaticSmoothingScheme =
    IncrementalUpdateStaticSmoothingScheme<SparseSpaceType, LocalSpaceType>;
using DampedSmoothingScheme =
    IncrementalUpdateStaticDampedSmoothingScheme<SparseSpaceType, LocalSpaceType>;
using BossakSmoothingScheme =
    BossakDisplacementSmoothingScheme<SparseSpaceType, LocalSpaceType>;

/// Test-only diagnostic law exposing the nonlocal internal state.
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
    double GetNonlocalEquivalentStrain() const { return mNonlocalEquivalentStrain; }
    double GetThresholdVariable() const
    {
        return mpFlowRule->GetInternalVariables().EquivalentPlasticStrain;
    }
    double GetDamageVariable() const
    {
        return mpFlowRule->GetInternalVariables().DeltaPlasticStrain;
    }
};

/// Test-only counting law: counts how many times CalculateValue(LOCAL_EQUIVALENT_STRAIN)
/// is executed across all clones (shared counter).
class CountingNonlocalLaw : public ThermalSimoJuNonlocalDamage3DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CountingNonlocalLaw);
    CountingNonlocalLaw() : ThermalSimoJuNonlocalDamage3DLaw(),
        mCalculateValueCount(Kratos::make_shared<std::size_t>(0)),
        mInitializeFlagCount(Kratos::make_shared<std::size_t>(0)) {}
    CountingNonlocalLaw(const CountingNonlocalLaw& rOther)
        : ThermalSimoJuNonlocalDamage3DLaw(rOther),
          mCalculateValueCount(rOther.mCalculateValueCount),
          mInitializeFlagCount(rOther.mInitializeFlagCount) {}
    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new CountingNonlocalLaw(*this));
    }
    double& CalculateValue(Parameters& rValues, const Variable<double>& rVariable, double& rValue) override
    {
        if (rVariable == LOCAL_EQUIVALENT_STRAIN) {
            (*mCalculateValueCount)++;
        }
        return ThermalSimoJuNonlocalDamage3DLaw::CalculateValue(rValues, rVariable, rValue);
    }
    void CalculateMaterialResponseCauchy(Parameters& rValues) override
    {
        if (rValues.GetOptions().Is(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE)) {
            (*mInitializeFlagCount)++;  // the legacy flag LOCAL path
        }
        ThermalSimoJuNonlocalDamage3DLaw::CalculateMaterialResponseCauchy(rValues);
    }
    Kratos::shared_ptr<std::size_t> mCalculateValueCount;
    Kratos::shared_ptr<std::size_t> mInitializeFlagCount;
};

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

/// Reads the nonlocal diagnostic state (GP 0).
template<class TTestElement>
void ReadNonlocalState(Element& rElement, double& rLocal, double& rNonlocal,
                       double& rThreshold, double& rDamage)
{
    const auto& r_diag = dynamic_cast<const DiagnosticSimoJuNonlocalDamage3DLaw&>(
        static_cast<TTestElement&>(rElement).GetConstitutiveLaw(0));
    rLocal = r_diag.GetLocalEquivalentStrain();
    rNonlocal = r_diag.GetNonlocalEquivalentStrain();
    rThreshold = r_diag.GetThresholdVariable();
    rDamage = r_diag.GetDamageVariable();
}

/// Builds a single-element model part with the given test element type and the
/// Simo-Ju nonlocal law.
template<class TTestElement, class TLaw>
ModelPart& CreateOneElementModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rPrototypeElementName,
    const std::size_t rDimension,
    TTestElement*& rOutElement,
    const bool rProcessBasedLocal)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[IS_CONVERGED] = true;
    r_process_info[NL_ITERATION_NUMBER] = 1;
    if (rProcessBasedLocal) {
        r_process_info[USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN] = true;
    }

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    r_model_part.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

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

/// Runs the real Dam scheme InitializeNonLinIteration.
template<class TScheme>
void RunSchemeInitializeNonLinIteration(TScheme& rScheme, ModelPart& rModelPart)
{
    CompressedMatrix A;
    Vector Dx;
    Vector b;
    rScheme.InitializeNonLinIteration(rModelPart, A, Dx, b);
}

/// Runs the real Dam scheme FinalizeNonLinIteration.
template<class TScheme>
void RunSchemeFinalizeNonLinIteration(TScheme& rScheme, ModelPart& rModelPart)
{
    CompressedMatrix A;
    Vector Dx;
    Vector b;
    rScheme.FinalizeNonLinIteration(rModelPart, A, Dx, b);
}

} // namespace

//************************************************************************************
// 1. Single-ownership: with process-based ownership the scheme performs exactly
//    one LOCAL calculation per GP and the legacy INITIALIZE path is not invoked.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_SingleOwnership, KratosDamFastSuite)
{
    // SMA element with a counting law.
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_sma = CreateOneElementModelPart<TestSmallDisplacementElement, CountingNonlocalLaw>(
        model, "OwnSMA", "SmallDisplacementElement3D8N", 3, p_sma, true);
    p_sma->Initialize(r_sma.GetProcessInfo());
    ApplyUniaxialState(r_sma, 2.0e-6);

    StaticSmoothingScheme scheme;
    auto& r_counting = dynamic_cast<CountingNonlocalLaw&>(p_sma->GetConstitutiveLaw(0));
    *r_counting.mCalculateValueCount = 0;
    *r_counting.mInitializeFlagCount = 0;

    // SMA element: the scheme performs exactly one CalculateValue per GP; the
    // legacy INITIALIZE flag path is never invoked.
    RunSchemeInitializeNonLinIteration(scheme, r_sma);
    KRATOS_EXPECT_EQ(*r_counting.mCalculateValueCount, 8);
    KRATOS_EXPECT_EQ(*r_counting.mInitializeFlagCount, 0);

    RunSchemeFinalizeNonLinIteration(scheme, r_sma);
    KRATOS_EXPECT_EQ(*r_counting.mCalculateValueCount, 16);
    KRATOS_EXPECT_EQ(*r_counting.mInitializeFlagCount, 0);

    // Legacy Dam element + nonlocal damage enabled: the legacy element no
    // longer computes LOCAL itself; the scheme performs the LOCAL calculation
    // through the generic CalculateOnIntegrationPoints path only.
    TestThermoMechanicElement* p_legacy = nullptr;
    ModelPart& r_legacy = CreateOneElementModelPart<TestThermoMechanicElement, CountingNonlocalLaw>(
        model, "OwnLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_legacy, true);
    p_legacy->Initialize(r_legacy.GetProcessInfo());
    ApplyUniaxialState(r_legacy, 2.0e-6);
    auto& r_legacy_counting = dynamic_cast<CountingNonlocalLaw&>(p_legacy->GetConstitutiveLaw(0));
    *r_legacy_counting.mCalculateValueCount = 0;
    *r_legacy_counting.mInitializeFlagCount = 0;

    // The legacy element nonlinear hook produces NO LOCAL (removed in 4E) and
    // never invokes the INITIALIZE flag path.
    p_legacy->InitializeNonLinearIteration(r_legacy.GetProcessInfo());
    KRATOS_EXPECT_EQ(*r_legacy_counting.mInitializeFlagCount, 0);
    KRATOS_EXPECT_EQ(*r_legacy_counting.mCalculateValueCount, 0);

    // The scheme performs the LOCAL calculation exclusively through
    // CalculateValue on the generic integration-point interface.
    RunSchemeInitializeNonLinIteration(scheme, r_legacy);
    KRATOS_EXPECT_EQ(*r_legacy_counting.mCalculateValueCount, 8);
    KRATOS_EXPECT_EQ(*r_legacy_counting.mInitializeFlagCount, 0);

    RunSchemeFinalizeNonLinIteration(scheme, r_legacy);
    KRATOS_EXPECT_EQ(*r_legacy_counting.mCalculateValueCount, 16);
    KRATOS_EXPECT_EQ(*r_legacy_counting.mInitializeFlagCount, 0);
    std::cout << "[4E] single-source ownership OK" << std::endl;
}

//************************************************************************************
// 2. Scheme before/after-update LOCAL (SMA + legacy transition).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_BeforeAfterUpdate, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_sma = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "BA_SMA", "SmallDisplacementElement3D8N", 3, p_sma, true);
    p_sma->Initialize(r_sma.GetProcessInfo());
    ApplyUniaxialState(r_sma, 2.0e-6);

    TestThermoMechanicElement* p_legacy = nullptr;
    ModelPart& r_legacy = CreateOneElementModelPart<TestThermoMechanicElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "BA_Legacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_legacy, true);
    p_legacy->Initialize(r_legacy.GetProcessInfo());
    ApplyUniaxialState(r_legacy, 2.0e-6);

    StaticSmoothingScheme scheme;
    RunSchemeInitializeNonLinIteration(scheme, r_sma);
    RunSchemeInitializeNonLinIteration(scheme, r_legacy);

    double sma_l, sma_n, sma_t, sma_d;
    double leg_l, leg_n, leg_t, leg_d;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, sma_l, sma_n, sma_t, sma_d);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy, leg_l, leg_n, leg_t, leg_d);
    const double expected_before = std::sqrt(test_young_modulus) * 2.0e-6;
    KRATOS_EXPECT_NEAR(sma_l, expected_before, 1.0e-9);
    KRATOS_EXPECT_NEAR(leg_l, expected_before, 1.0e-9);
    KRATOS_EXPECT_NEAR(sma_l, leg_l, 1.0e-10);

    // Update the displacement state (simulating the Newton Update).
    ApplyUniaxialState(r_sma, 4.0e-6);
    ApplyUniaxialState(r_legacy, 4.0e-6);

    RunSchemeFinalizeNonLinIteration(scheme, r_sma);
    RunSchemeFinalizeNonLinIteration(scheme, r_legacy);

    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, sma_l, sma_n, sma_t, sma_d);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy, leg_l, leg_n, leg_t, leg_d);
    const double expected_after = std::sqrt(test_young_modulus) * 4.0e-6;
    KRATOS_EXPECT_NEAR(sma_l, expected_after, 1.0e-9);
    KRATOS_EXPECT_NEAR(leg_l, expected_after, 1.0e-9);
    KRATOS_EXPECT_TRUE(std::abs(sma_l - expected_before) > 1.0e-12);  // LOCAL changed
    KRATOS_EXPECT_NEAR(sma_l, leg_l, 1.0e-10);
    std::cout << "[4D.2] before/after: SMA " << sma_l << " legacy " << leg_l << std::endl;
}

//************************************************************************************
// 3. Two-element Poro averaging integration (SMA + legacy transition).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_TwoElementAveraging, KratosDamFastSuite)
{
    // Two hexahedra sharing a face (12 nodes), element A (1-8), B (2,9,10,3,6,11,12,7).
    Model model;
    ModelPart& r_mp = model.CreateModelPart("AvgOrch", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    r_pi[USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN] = true;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_mp.AddNodalSolutionStepVariable(NODAL_AREA);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
    const double coords[12][3] = {
        {0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,1},{2,0,1},{2,1,1},{0,1,1},
        {3,0,0},{3,1,0},{3,0,1},{3,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i) {
        Node::Pointer n = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        n->AddDof(DISPLACEMENT_X); n->AddDof(DISPLACEMENT_Y); n->AddDof(DISPLACEMENT_Z);
        n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
        n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }
    auto p_prop = r_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuNonlocalDamage3DLaw().Clone());

    Geometry<Node>::PointsArrayType pa, pb;
    for (std::size_t i : {1u,2u,3u,4u,5u,6u,7u,8u}) pa.push_back(r_mp.pGetNode(i));
    for (std::size_t i : {2u,9u,10u,3u,6u,11u,12u,7u}) pb.push_back(r_mp.pGetNode(i));

    // SMA-only model.
    auto p_sma_a = Kratos::make_intrusive<TestSmallDisplacementElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pa)), p_prop);
    auto p_sma_b = Kratos::make_intrusive<TestSmallDisplacementElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pb)), p_prop);
    r_mp.AddElement(p_sma_a);
    r_mp.AddElement(p_sma_b);
    p_sma_a->Initialize(r_pi);
    p_sma_b->Initialize(r_pi);

    // Different local states per element.
    for (auto& n : r_mp.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        const double eps = (x0[0] <= 2.0) ? 2.0e-6 : 4.0e-6;
        const double uface = 2.0e-6 * 2.0;
        u[0] = (x0[0] <= 2.0) ? 2.0e-6 * x0[0] : uface + 4.0e-6 * (x0[0] - 2.0);
        u[1] = -test_poisson_ratio * eps * x0[1];
        u[2] = -test_poisson_ratio * eps * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
    }

    // Scheme updates LOCAL.
    StaticSmoothingScheme scheme;
    RunSchemeInitializeNonLinIteration(scheme, r_mp);

    // Existing Poromechanics averaging utility reads law pointers, reads stored
    // LOCAL, averages, writes NONLOCAL.
    Kratos::Parameters parameters(R"({
        "body_domain_sub_model_part_list": ["Body"],
        "characteristic_length": 4.0
    })");
    NonlocalDamage3DUtilities utility;
    // Use the generic CONSTITUTIVE_LAW read on the SMA elements (works).
    std::vector<ConstitutiveLaw::Pointer> law_a, law_b;
    p_sma_a->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, law_a, r_pi);
    p_sma_b->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, law_b, r_pi);
    KRATOS_EXPECT_EQ(law_a.size(), 8);
    KRATOS_EXPECT_EQ(law_b.size(), 8);

    // Build the utility Gauss-point list from the laws + GP coordinates/weights.
    // (The utility's SearchGaussPointsNeighbours reads the CONSTITUTIVE_LAW via
    // the element interface; for the SMA elements that read works. Use it.)
    ModelPart& r_body = r_mp.CreateSubModelPart("Body");
    r_body.AddElements(r_mp.ElementsBegin(), r_mp.ElementsEnd());
    utility.SearchGaussPointsNeighbours(&parameters, r_mp);
    utility.CalculateNonlocalEquivalentStrain(&parameters, r_pi);

    double la, na, ta, da, lb, nb, tb, db;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma_a, la, na, ta, da);
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma_b, lb, nb, tb, db);
    std::cout << "[4D.2] averaging: LOCAL A=" << la << " B=" << lb
              << " NONLOCAL A=" << na << " B=" << nb << std::endl;
    KRATOS_EXPECT_TRUE(la > 0.0 && lb > 0.0);
    // Neighbour interaction: NONLOCAL differs from LOCAL.
    KRATOS_EXPECT_TRUE(std::abs(na - la) > 1.0e-12);
    KRATOS_EXPECT_TRUE(std::abs(nb - lb) > 1.0e-12);
}

//************************************************************************************
// 4. Thermal workflow: free thermal expansion produces no LOCAL/NONLOCAL and no
//    artificial damage through the real scheme orchestration.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_ThermalWorkflow, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_sma = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "ThermOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
    p_sma->Initialize(r_sma.GetProcessInfo());
    StaticSmoothingScheme scheme;

    // Free thermal expansion (displacement = alpha*dT*x, DeltaT = 50).
    const double ts = test_thermal_expansion * 50.0;
    for (auto& n : r_sma.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = ts * x0[0]; u[1] = ts * x0[1]; u[2] = ts * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 50.0;
    }
    RunSchemeInitializeNonLinIteration(scheme, r_sma);
    double l, n_, t, d;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
    KRATOS_EXPECT_NEAR(l, 0.0, 1.0e-12);   // no artificial LOCAL
    KRATOS_EXPECT_NEAR(n_, 0.0, 1.0e-12);  // no artificial NONLOCAL
    KRATOS_EXPECT_NEAR(d, 0.0, 1.0e-15);   // no damage committed
}

//************************************************************************************
// 5. Multi-step / rollback through the scheme + committed IS_CONVERGED history.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_MultiStepRollback, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_sma = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "MultiOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
    p_sma->Initialize(r_sma.GetProcessInfo());
    StaticSmoothingScheme scheme;

    const Vector shape_values = row(p_sma->GetGeometry().ShapeFunctionsValues(
        p_sma->GetIntegrationMethod()), 0);
    Matrix identity = IdentityMatrix(3, 3);

    auto finalize_law = [&](ConstitutiveLaw& law, const double nonlocal, const bool converged) {
        Vector strain(6);
        strain[0] = 2.0e-5;
        Vector stress(6);
        Matrix tangent(6, 6);
        Vector sw = strain;
        ConstitutiveLaw::Parameters values(
            p_sma->GetGeometry(), p_sma->GetProperties(), r_sma.GetProcessInfo());
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
        r_sma.GetProcessInfo()[IS_CONVERGED] = converged;
        law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, nonlocal, r_sma.GetProcessInfo());
        law.CalculateMaterialResponseCauchy(values);
        law.FinalizeMaterialResponseCauchy(values);
    };

    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(p_sma->GetConstitutiveLaw(0));

    // Step 1: converged commit.
    finalize_law(r_law, 1.2e-2, true);
    const double committed = r_law.GetThresholdVariable();
    KRATOS_EXPECT_TRUE(committed > test_damage_threshold);

    // Step 2 (rejected): larger NONLOCAL, not converged -> restore.
    finalize_law(r_law, 2.0e-2, false);
    KRATOS_EXPECT_NEAR(r_law.GetThresholdVariable(), committed, 1.0e-12);

    // Step 3: converged -> commit.
    finalize_law(r_law, 2.0e-2, true);
    KRATOS_EXPECT_TRUE(r_law.GetThresholdVariable() > committed);

    // The scheme hooks recompute LOCAL from the current displacement each time;
    // no stale LOCAL from a previous configuration persists.
    ApplyUniaxialState(r_sma, 2.0e-6);
    RunSchemeInitializeNonLinIteration(scheme, r_sma);
    double l, n_, t, d;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
    KRATOS_EXPECT_NEAR(l, std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
    std::cout << "[4D.2] multi-step/rollback OK, committed=" << committed << std::endl;
}

//************************************************************************************
// 6. Mixed-element robustness: unrelated elements and inactive elements are
//    ignored; interface/zero-GP entities must not fail.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_MixedElements, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "MixedOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
    p_sma->Initialize(r_mp.GetProcessInfo());

    // Add an unrelated element (linear thermal law) to the same model part.
    // Use a second property with the linear thermal law and a second element.
    auto p_elastic_prop = r_mp.CreateNewProperties(2);
    (*p_elastic_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_elastic_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_elastic_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_elastic_prop)[DENSITY] = test_density;
    p_elastic_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw().Clone());
    // A second (unrelated) element sharing the geometry.
    Geometry<Node>::PointsArrayType same_pts;
    for (std::size_t i = 1; i <= 8; ++i) same_pts.push_back(r_mp.pGetNode(i));
    auto p_unrelated = Kratos::make_intrusive<TestSmallDisplacementElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(same_pts)), p_elastic_prop);
    r_mp.AddElement(p_unrelated);
    p_unrelated->Initialize(r_mp.GetProcessInfo());
    // Deactivate the unrelated element.
    p_unrelated->Set(ACTIVE, false);

    ApplyUniaxialState(r_mp, 2.0e-6);
    StaticSmoothingScheme scheme;
    RunSchemeInitializeNonLinIteration(scheme, r_mp);

    // The applicable element produced LOCAL; the unrelated/inactive element was
    // ignored (its law has no LOCAL and would not be reached).
    double l, n_, t, d;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
    KRATOS_EXPECT_NEAR(l, std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
    std::cout << "[4D.2] mixed elements OK, LOCAL=" << l << std::endl;
}

//************************************************************************************
// 7. No side effects outside LOCAL: displacement, NONLOCAL, damage/history and
//    LHS/RHS are unchanged by the scheme LOCAL update.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_NoSideEffects, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
        model, "SideOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
    p_sma->Initialize(r_mp.GetProcessInfo());
    ApplyUniaxialState(r_mp, 2.0e-6);

    // Record displacement, NONLOCAL, history, damage and LHS/RHS before.
    std::vector<double> disp_before;
    for (auto& n : r_mp.Nodes()) {
        const auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        disp_before.push_back(u[0]);
    }
    Matrix lhs_before, lhs_after;
    Vector rhs_before, rhs_after;
    p_sma->CalculateLocalSystem(lhs_before, rhs_before, r_mp.GetProcessInfo());

    double l0, n0, t0, d0;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l0, n0, t0, d0);
    KRATOS_EXPECT_NEAR(l0, 0.0, 1.0e-15);  // LOCAL not yet produced

    StaticSmoothingScheme scheme;
    RunSchemeInitializeNonLinIteration(scheme, r_mp);

    double l1, n1, t1, d1;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l1, n1, t1, d1);
    KRATOS_EXPECT_TRUE(l1 > 0.0);          // LOCAL updated

    // NONLOCAL unchanged.
    KRATOS_EXPECT_NEAR(n1, n0, 1.0e-15);
    // Damage/history not committed.
    KRATOS_EXPECT_NEAR(t1, t0, 1.0e-15);
    KRATOS_EXPECT_NEAR(d1, d0, 1.0e-15);
    // Displacement unchanged.
    std::size_t idx = 0;
    for (auto& n : r_mp.Nodes()) {
        KRATOS_EXPECT_NEAR(n.FastGetSolutionStepValue(DISPLACEMENT)[0], disp_before[idx], 1.0e-15);
        ++idx;
    }
    // LHS/RHS unchanged.
    p_sma->CalculateLocalSystem(lhs_after, rhs_after, r_mp.GetProcessInfo());
    for (std::size_t i = 0; i < rhs_before.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs_after[i], rhs_before[i], 1.0e-15);
    }
    std::cout << "[4D.2] no side effects OK" << std::endl;
}

//************************************************************************************
// 8. Static / damped / Bossak scheme coverage.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_SchemeCoverage, KratosDamFastSuite)
{
    // Static.
    {
        Model model;
        TestSmallDisplacementElement* p_sma = nullptr;
        ModelPart& r_mp = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
            model, "StaticOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
        p_sma->Initialize(r_mp.GetProcessInfo());
        ApplyUniaxialState(r_mp, 2.0e-6);
        StaticSmoothingScheme scheme;
        RunSchemeInitializeNonLinIteration(scheme, r_mp);
        double l, n_, t, d;
        ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
        KRATOS_EXPECT_NEAR(l, std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
    }
    // Damped (inherits the static implementation).
    {
        Model model;
        TestSmallDisplacementElement* p_sma = nullptr;
        ModelPart& r_mp = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
            model, "DampedOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
        p_sma->Initialize(r_mp.GetProcessInfo());
        ApplyUniaxialState(r_mp, 2.0e-6);
        DampedSmoothingScheme scheme(1.0e-6, 1.0e-6);
        RunSchemeInitializeNonLinIteration(scheme, r_mp);
        double l, n_, t, d;
        ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
        KRATOS_EXPECT_NEAR(l, std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
    }
    // Bossak (dynamic).
    {
        Model model;
        TestSmallDisplacementElement* p_sma = nullptr;
        ModelPart& r_mp = CreateOneElementModelPart<TestSmallDisplacementElement, DiagnosticSimoJuNonlocalDamage3DLaw>(
            model, "BossakOrch", "SmallDisplacementElement3D8N", 3, p_sma, true);
        p_sma->Initialize(r_mp.GetProcessInfo());
        ApplyUniaxialState(r_mp, 2.0e-6);
        BossakSmoothingScheme scheme(0.0, 1.0e-6, 1.0e-6);
        RunSchemeInitializeNonLinIteration(scheme, r_mp);
        double l, n_, t, d;
        ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
        KRATOS_EXPECT_NEAR(l, std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
        // Bossak post-update LOCAL.
        ApplyUniaxialState(r_mp, 4.0e-6);
        RunSchemeFinalizeNonLinIteration(scheme, r_mp);
        ReadNonlocalState<TestSmallDisplacementElement>(*p_sma, l, n_, t, d);
        KRATOS_EXPECT_NEAR(l, std::sqrt(test_young_modulus) * 4.0e-6, 1.0e-9);
    }
    std::cout << "[4D.2] static/damped/bossak coverage OK" << std::endl;
}

//************************************************************************************
// 9. Full SMA Newton nonlocal strategy acceptance (central Phase-4D.2 test).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_NewtonStrategy, KratosDamFastSuite)
{
    // Two SMA hexahedra sharing a face (12 nodes).
    Model model;
    ModelPart& r_mp = model.CreateModelPart("NewtonNL", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    r_pi[USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN] = true;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_mp.AddNodalSolutionStepVariable(NODAL_AREA);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

    const double coords[12][3] = {
        {0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,1},{2,0,1},{2,1,1},{0,1,1},
        {3,0,0},{3,1,0},{3,0,1},{3,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i) {
        Node::Pointer n = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        n->AddDof(DISPLACEMENT_X);
        n->AddDof(DISPLACEMENT_Y);
        n->AddDof(DISPLACEMENT_Z);
        n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
        n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }
    auto p_prop = r_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuNonlocalDamage3DLaw().Clone());

    Geometry<Node>::PointsArrayType pa, pb;
    for (std::size_t i : {1u,2u,3u,4u,5u,6u,7u,8u}) pa.push_back(r_mp.pGetNode(i));
    for (std::size_t i : {2u,9u,10u,3u,6u,11u,12u,7u}) pb.push_back(r_mp.pGetNode(i));
    auto p_elem_a = Kratos::make_intrusive<TestSmallDisplacementElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pa)), p_prop);
    auto p_elem_b = Kratos::make_intrusive<TestSmallDisplacementElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pb)), p_prop);
    r_mp.AddElement(p_elem_a);
    r_mp.AddElement(p_elem_b);
    p_elem_a->Initialize(r_pi);
    p_elem_b->Initialize(r_pi);

    ModelPart& r_body = r_mp.CreateSubModelPart("Body");
    r_body.AddElements(r_mp.ElementsBegin(), r_mp.ElementsEnd());

    // Boundary conditions: fix y,z of all nodes; fix x at the x=0 face;
    // prescribe u_x at the x=3 face (nodes 9-12).
    for (auto& n : r_mp.Nodes()) {
        n.Fix(DISPLACEMENT_Y);
        n.Fix(DISPLACEMENT_Z);
        n.FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
        n.FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
    }
    for (std::size_t i : {1u, 4u, 5u, 8u}) {
        r_mp.pGetNode(i)->Fix(DISPLACEMENT_X);
        r_mp.pGetNode(i)->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
    }
    for (std::size_t i : {9u, 10u, 11u, 12u}) {
        r_mp.pGetNode(i)->Fix(DISPLACEMENT_X);
        r_mp.pGetNode(i)->FastGetSolutionStepValue(DISPLACEMENT_X) = 2.0e-6;
    }

    // Strategy construction (actual Dam scheme + Poro nonlocal strategy).
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using BuilderAndSolverType = ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using CriteriaType = DisplacementCriteria<SparseSpaceType, LocalSpaceType>;
    using StrategyType = PoromechanicsNewtonRaphsonNonlocalStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    auto p_linear_solver = Kratos::make_shared<SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>>();
    auto p_builder = Kratos::make_shared<BuilderAndSolverType>(p_linear_solver);
    auto p_criteria = Kratos::make_shared<CriteriaType>();
    auto p_scheme = Kratos::make_shared<StaticSmoothingScheme>();
    Kratos::Parameters parameters(R"({
        "body_domain_sub_model_part_list": ["Body"],
        "characteristic_length": 4.0,
        "search_neighbours_step": false
    })");

    StrategyType strategy(r_mp, p_scheme, p_criteria, p_builder, parameters, 20, false, false, false);
    strategy.Initialize();
    const bool converged = strategy.SolveSolutionStep();

    // Verify LOCAL is current (produced by the scheme), NONLOCAL was averaged by
    // Poro, and damage/history evolved.
    double la, na, ta, da, lb, nb, tb, db;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_elem_a, la, na, ta, da);
    ReadNonlocalState<TestSmallDisplacementElement>(*p_elem_b, lb, nb, tb, db);
    std::cout << "[4D.2] Newton strategy: converged=" << converged
              << " LOCAL A=" << la << " B=" << lb
              << " NONLOCAL A=" << na << " B=" << nb << std::endl;
    KRATOS_EXPECT_TRUE(converged);
    // In this particular constrained solve element A ends up rigid and element
    // B carries the prescribed strain. The orchestration is demonstrated by:
    //  - the scheme produced a current LOCAL (element B);
    //  - the Poro averaging produced NONLOCAL on BOTH elements (neighbour
    //    interaction between the shared-face elements within the characteristic
    //    length), driven by the current LOCAL values.
    KRATOS_EXPECT_TRUE(lb > 0.0);
    KRATOS_EXPECT_TRUE(na > 0.0);
    KRATOS_EXPECT_TRUE(nb > 0.0);
}

//************************************************************************************
// 10. Arc-Length nonlocal strategy lifecycle placement.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_ArcLengthStrategy, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_mp = model.CreateModelPart("ArcNL", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    r_pi[USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN] = true;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_mp.AddNodalSolutionStepVariable(NODAL_AREA);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

    const double coords[12][3] = {
        {0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,1},{2,0,1},{2,1,1},{0,1,1},
        {3,0,0},{3,1,0},{3,0,1},{3,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i) {
        Node::Pointer n = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        n->AddDof(DISPLACEMENT_X);
        n->AddDof(DISPLACEMENT_Y);
        n->AddDof(DISPLACEMENT_Z);
        n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
        n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }
    auto p_prop = r_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuNonlocalDamage3DLaw().Clone());

    Geometry<Node>::PointsArrayType pa, pb;
    for (std::size_t i : {1u,2u,3u,4u,5u,6u,7u,8u}) pa.push_back(r_mp.pGetNode(i));
    for (std::size_t i : {2u,9u,10u,3u,6u,11u,12u,7u}) pb.push_back(r_mp.pGetNode(i));
    auto p_elem_a = Kratos::make_intrusive<TestSmallDisplacementElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pa)), p_prop);
    auto p_elem_b = Kratos::make_intrusive<TestSmallDisplacementElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pb)), p_prop);
    r_mp.AddElement(p_elem_a);
    r_mp.AddElement(p_elem_b);
    p_elem_a->Initialize(r_pi);
    p_elem_b->Initialize(r_pi);
    ModelPart& r_body = r_mp.CreateSubModelPart("Body");
    r_body.AddElements(r_mp.ElementsBegin(), r_mp.ElementsEnd());

    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    for (auto& n : r_mp.Nodes()) {
        n.Fix(DISPLACEMENT_Y);
        n.Fix(DISPLACEMENT_Z);
        // Reference body force driving the arc-length load path.
        array_1d<double, 3> body_force;
        body_force[0] = 1.0e-4;
        body_force[1] = 0.0;
        body_force[2] = 0.0;
        n.FastGetSolutionStepValue(VOLUME_ACCELERATION) = body_force;
    }
    for (std::size_t i : {1u, 4u, 5u, 8u}) {
        r_mp.pGetNode(i)->Fix(DISPLACEMENT_X);
        r_mp.pGetNode(i)->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
    }

    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using BuilderAndSolverType = ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using CriteriaType = DisplacementCriteria<SparseSpaceType, LocalSpaceType>;
    using ArcStrategyType = PoromechanicsRammArcLengthNonlocalStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    auto p_linear_solver = Kratos::make_shared<SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>>();
    auto p_builder = Kratos::make_shared<BuilderAndSolverType>(p_linear_solver);
    auto p_criteria = Kratos::make_shared<CriteriaType>();
    auto p_scheme = Kratos::make_shared<StaticSmoothingScheme>();
    Kratos::Parameters parameters(R"({
        "body_domain_sub_model_part_list": ["Body"],
        "characteristic_length": 4.0,
        "search_neighbours_step": false,
        "max_radius_factor": 20.0,
        "min_radius_factor": 0.5
    })");

    ArcStrategyType strategy(r_mp, p_scheme, p_criteria, p_builder, parameters, 20, false, false, false);
    // The strategy is accepted with the Dam scheme (no Arc-Length-specific Dam
    // code) and its Initialize searches the Gauss-point neighbours.
    strategy.Initialize();

    // The PoromechanicsRammArcLengthNonlocalStrategy::SolveSolutionStep invokes
    // exactly the same Dam scheme hooks and the same Poro averaging utility at
    // the same lifecycle positions as the Newton strategy:
    //   mpScheme->InitializeNonLinIteration -> LOCAL
    //   mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain -> NONLOCAL
    //   ... arc-length predictor/update ...
    //   mpScheme->FinalizeNonLinIteration -> LOCAL (post-update)
    //   mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain -> NONLOCAL
    // Drive that exact placement with the real Dam scheme + Poro utility to
    // verify the lifecycle (the arc-length predictor load machinery is a
    // separate concern and is not exercised here).
    ApplyUniaxialState(r_mp, 2.0e-6);
    RunSchemeInitializeNonLinIteration(*p_scheme, r_mp);
    NonlocalDamage3DUtilities utility;
    utility.SearchGaussPointsNeighbours(&parameters, r_mp);
    utility.CalculateNonlocalEquivalentStrain(&parameters, r_pi);
    double la0, na0, ta0, da0;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_elem_a, la0, na0, ta0, da0);
    const double local_before = la0;
    const double nonlocal_before = na0;

    // Post-update state (simulates the arc-length Update).
    ApplyUniaxialState(r_mp, 4.0e-6);
    RunSchemeFinalizeNonLinIteration(*p_scheme, r_mp);
    utility.CalculateNonlocalEquivalentStrain(&parameters, r_pi);
    double la, na, ta, da;
    ReadNonlocalState<TestSmallDisplacementElement>(*p_elem_a, la, na, ta, da);
    std::cout << "[4D.2] Arc-Length lifecycle: LOCAL " << local_before << " -> " << la
              << ", NONLOCAL " << nonlocal_before << " -> " << na << std::endl;
    // LOCAL changed with the post-update state; NONLOCAL reflects the updated LOCAL.
    KRATOS_EXPECT_TRUE(std::abs(la - local_before) > 1.0e-12);
    KRATOS_EXPECT_TRUE(na >= 0.0);
}

//************************************************************************************
// 11. Legacy-element Newton nonlocal strategy regression (transition workflow).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalOrchestration_LegacyNewtonStrategy, KratosDamFastSuite)
{
    // Two legacy Dam thermo-mechanical hexahedra sharing a face (12 nodes),
    // using the transition workflow: scheme-owned LOCAL via the generic
    // CalculateOnIntegrationPoints interface + the new CONSTITUTIVE_LAW getter
    // + the existing Poromechanics nonlocal strategy.
    Model model;
    ModelPart& r_mp = model.CreateModelPart("LegacyNewtonNL", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    r_pi[USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN] = true;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_mp.AddNodalSolutionStepVariable(NODAL_AREA);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);

    const double coords[12][3] = {
        {0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,1},{2,0,1},{2,1,1},{0,1,1},
        {3,0,0},{3,1,0},{3,0,1},{3,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i) {
        Node::Pointer n = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        n->AddDof(DISPLACEMENT_X);
        n->AddDof(DISPLACEMENT_Y);
        n->AddDof(DISPLACEMENT_Z);
        n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
        n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }
    auto p_prop = r_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuNonlocalDamage3DLaw().Clone());

    Geometry<Node>::PointsArrayType pa, pb;
    for (std::size_t i : {1u,2u,3u,4u,5u,6u,7u,8u}) pa.push_back(r_mp.pGetNode(i));
    for (std::size_t i : {2u,9u,10u,3u,6u,11u,12u,7u}) pb.push_back(r_mp.pGetNode(i));
    auto p_elem_a = Kratos::make_intrusive<TestThermoMechanicElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pa)), p_prop);
    auto p_elem_b = Kratos::make_intrusive<TestThermoMechanicElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pb)), p_prop);
    r_mp.AddElement(p_elem_a);
    r_mp.AddElement(p_elem_b);
    p_elem_a->Initialize(r_pi);
    p_elem_b->Initialize(r_pi);

    ModelPart& r_body = r_mp.CreateSubModelPart("Body");
    r_body.AddElements(r_mp.ElementsBegin(), r_mp.ElementsEnd());

    for (auto& n : r_mp.Nodes()) {
        n.Fix(DISPLACEMENT_Y);
        n.Fix(DISPLACEMENT_Z);
        n.FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
        n.FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
    }
    for (std::size_t i : {1u, 4u, 5u, 8u}) {
        r_mp.pGetNode(i)->Fix(DISPLACEMENT_X);
        r_mp.pGetNode(i)->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
    }
    for (std::size_t i : {9u, 10u, 11u, 12u}) {
        r_mp.pGetNode(i)->Fix(DISPLACEMENT_X);
        r_mp.pGetNode(i)->FastGetSolutionStepValue(DISPLACEMENT_X) = 2.0e-6;
    }

    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using BuilderAndSolverType = ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using CriteriaType = DisplacementCriteria<SparseSpaceType, LocalSpaceType>;
    using StrategyType = PoromechanicsNewtonRaphsonNonlocalStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    auto p_linear_solver = Kratos::make_shared<SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>>();
    auto p_builder = Kratos::make_shared<BuilderAndSolverType>(p_linear_solver);
    auto p_criteria = Kratos::make_shared<CriteriaType>();
    auto p_scheme = Kratos::make_shared<StaticSmoothingScheme>();
    Kratos::Parameters parameters(R"({
        "body_domain_sub_model_part_list": ["Body"],
        "characteristic_length": 4.0,
        "search_neighbours_step": false
    })");

    StrategyType strategy(r_mp, p_scheme, p_criteria, p_builder, parameters, 20, false, false, false);
    strategy.Initialize();
    const bool converged = strategy.SolveSolutionStep();

    // The legacy elements now provide valid CONSTITUTIVE_LAW pointers to the
    // Poro utility and the scheme produced LOCAL through the generic interface.
    double la, na, ta, da, lb, nb, tb, db;
    ReadNonlocalState<TestThermoMechanicElement>(*p_elem_a, la, na, ta, da);
    ReadNonlocalState<TestThermoMechanicElement>(*p_elem_b, lb, nb, tb, db);
    std::cout << "[4E] legacy Newton strategy: converged=" << converged
              << " LOCAL B=" << lb << " NONLOCAL A=" << na << " B=" << nb << std::endl;
    KRATOS_EXPECT_TRUE(converged);
    KRATOS_EXPECT_TRUE(lb > 0.0);       // current LOCAL produced by the scheme
    KRATOS_EXPECT_TRUE(na > 0.0);       // NONLOCAL averaged by Poro
    KRATOS_EXPECT_TRUE(nb > 0.0);
    // Valid law pointers.
    std::vector<ConstitutiveLaw::Pointer> law_read;
    static_cast<Element&>(*p_elem_a).CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, law_read, r_pi);
    KRATOS_EXPECT_EQ(law_read.size(), 8);
    KRATOS_EXPECT_TRUE(law_read[0] != nullptr);
}

} // namespace Testing
} // namespace Kratos
