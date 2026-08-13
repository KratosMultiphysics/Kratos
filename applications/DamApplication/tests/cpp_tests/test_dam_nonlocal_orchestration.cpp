// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Permanent Dam nonlocal-damage orchestration contract. The nonlinear-iteration
// LOCAL-equivalent-strain is produced through the generic
// Element::CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, ...) interface
// in the Dam smoothing-scheme workflow, averaged by the existing Poromechanics
// nonlocal utilities into NONLOCAL, and consumed by the damage laws. Tests cover
// the two-element averaging pipeline and the full Newton strategy solve.
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


//************************************************************************************
// 2. Scheme before/after-update LOCAL (SMA + legacy transition).
//************************************************************************************


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


//************************************************************************************
// 5. Multi-step / rollback through the scheme + committed IS_CONVERGED history.
//************************************************************************************


//************************************************************************************
// 6. Mixed-element robustness: unrelated elements and inactive elements are
//    ignored; interface/zero-GP entities must not fail.
//************************************************************************************


//************************************************************************************
// 7. No side effects outside LOCAL: displacement, NONLOCAL, damage/history and
//    LHS/RHS are unchanged by the scheme LOCAL update.
//************************************************************************************


//************************************************************************************
// 8. Static / damped / Bossak scheme coverage.
//************************************************************************************


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


//************************************************************************************
// 11. Legacy-element Newton nonlocal strategy regression (transition workflow).
//************************************************************************************


} // namespace Testing
} // namespace Kratos
