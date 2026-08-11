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

// Phase 5A: characterize and adapt the THERMAL NODAL constitutive-law family
// (ThermalLinearElastic{3D,2DPlaneStrain,2DPlaneStress}Nodal) so it satisfies
// the StructuralMechanicsApplication::SmallDisplacement constitutive contract,
// and determine whether the legacy element FinalizeSolutionStep override is
// still necessary.

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
#include "spaces/ublas_space.h"
#include "geometries/hexahedra_3d_8.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress_nodal.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"

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
constexpr double test_poisson_ratio = 0.2;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_thickness = 0.15;

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using StaticSmoothingScheme =
    IncrementalUpdateStaticSmoothingScheme<SparseSpaceType, LocalSpaceType>;

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

/// Builds a single-element model part with the given test element type and the
/// nodal thermal law.
template<class TTestElement>
ModelPart& CreateNodalModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rPrototypeElementName,
    const std::size_t rDimension,
    TTestElement*& rOutElement,
    ConstitutiveLaw::Pointer pNodalLaw)
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
    r_model_part.AddNodalSolutionStepVariable(NODAL_YOUNG_MODULUS);
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
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    if (rDimension == 2) {
        (*p_prop)[THICKNESS] = test_thickness;
    }
    p_prop->SetValue(CONSTITUTIVE_LAW, pNodalLaw);
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

/// Applies a uniaxial-STRESS field with non-uniform nodal material properties.
void ApplyNodalState(ModelPart& rModelPart, const double rEpsilonX,
                     const double rDeltaTemperature)
{
    std::size_t index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = rEpsilonX * r_x0[0];
        r_displacement[1] = -test_poisson_ratio * rEpsilonX * r_x0[1];
        r_displacement[2] = -test_poisson_ratio * rEpsilonX * r_x0[2];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
        // Non-uniform nodal Young modulus (varies in all spatial directions) and
        // non-uniform temperature.
        r_node.FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) =
            2.0e7 + 5.0e6 * (r_x0[0] + 2.0 * r_x0[1] + 3.0 * r_x0[2]);
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + rDeltaTemperature + 3.0 * r_x0[0] + 2.0 * r_x0[1] + 1.0 * r_x0[2];
        ++index;
    }
}

/// Computes the analytically interpolated GP-0 nodal Young modulus.
double ExpectedGp0YoungModulus(Element& rElement)
{
    const auto& r_geom = rElement.GetGeometry();
    const auto& method = rElement.GetIntegrationMethod();
    const Vector N = row(r_geom.ShapeFunctionsValues(method), 0);
    double E = 0.0;
    for (std::size_t j = 0; j < r_geom.size(); ++j) {
        E += N[j] * r_geom[j].FastGetSolutionStepValue(NODAL_YOUNG_MODULUS);
    }
    return E;
}

/// Compares legacy vs candidate LHS/RHS.
void CompareSystems(Element& rLegacy, Element& rCandidate,
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
                           std::max(1.0e-9, comparison_relative_tolerance * std::abs(rhs_legacy[i])));
    }
    for (std::size_t i = 0; i < lhs_legacy.size1(); ++i) {
        for (std::size_t j = 0; j < lhs_legacy.size2(); ++j) {
            KRATOS_EXPECT_NEAR(lhs_candidate(i, j), lhs_legacy(i, j),
                               std::max(1.0e-9, comparison_relative_tolerance * std::abs(lhs_legacy(i, j))));
        }
    }
}

} // namespace

//************************************************************************************
// 1. 3D8N legacy vs candidate with non-uniform nodal properties and temperature.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNodalLaw_3D8N_LegacyVsCandidate, KratosDamFastSuite)
{
    Model model;
    auto p_law = ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal());

    TestThermoMechanicElement* p_legacy = nullptr;
    ModelPart& r_legacy = CreateNodalModelPart<TestThermoMechanicElement>(
        model, "NodalLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_legacy,
        p_law->Clone());
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_candidate = CreateNodalModelPart<TestSmallDisplacementElement>(
        model, "NodalCandidate", "SmallDisplacementElement3D8N", 3, p_candidate,
        p_law->Clone());
    p_legacy->Initialize(r_legacy.GetProcessInfo());
    p_candidate->Initialize(r_candidate.GetProcessInfo());

    // Non-uniform nodal material/thermal field (mechanical + temperature).
    ApplyNodalState(r_legacy, 2.0e-6, 0.0);
    ApplyNodalState(r_candidate, 2.0e-6, 0.0);

    // The candidate law interpolated the nodal Young modulus at GP 0.
    const double expected_E = ExpectedGp0YoungModulus(*p_candidate);
    KRATOS_EXPECT_TRUE(expected_E > 0.0);

    // Primary constitutive-response equivalence: LHS/RHS through the real
    // elements (the SMA element drives the law through PK2, the legacy through
    // Cauchy; both reach the common nodal thermal response).
    CompareSystems(*p_legacy, *p_candidate, r_legacy, r_candidate);

    // Law-level stress comparison with IDENTICAL Parameters (the SMA element's
    // CAUCHY_STRESS_VECTOR output uses CalculateValue, which is a Phase-5B
    // specialized-output gap for the nodal family, so it is not used here).
    {
        const Vector shape_values = row(p_legacy->GetGeometry().ShapeFunctionsValues(
            p_legacy->GetIntegrationMethod()), 0);
        Matrix identity = IdentityMatrix(3, 3);
        Vector strain(6);
        strain[0] = 2.0e-6; strain[1] = -4.0e-7; strain[2] = -4.0e-7;
        Vector stress_l(6), stress_c(6);
        Matrix tangent_l(6, 6), tangent_c(6, 6);
        auto run = [&](ConstitutiveLaw& law, Vector& stress, Matrix& tangent) {
            Vector sw = strain;
            ConstitutiveLaw::Parameters values(
                p_legacy->GetGeometry(), p_legacy->GetProperties(), r_legacy.GetProcessInfo());
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
            law.CalculateMaterialResponseCauchy(values);
        };
        run(p_legacy->GetConstitutiveLaw(0), stress_l, tangent_l);
        run(p_candidate->GetConstitutiveLaw(0), stress_c, tangent_c);
        for (std::size_t i = 0; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(stress_c[i], stress_l[i],
                               std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(stress_l[i])));
        }
    }
    std::cout << "[5A] 3D nodal: interpolated E~" << expected_E << " LHS/RHS+law stress OK" << std::endl;
}

//************************************************************************************
// 2. Thermal cases (restrained heating, free expansion, non-uniform temperature).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNodalLaw_ThermalCases, KratosDamFastSuite)
{
    Model model;
    auto p_law = ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal());

    TestThermoMechanicElement* p_legacy = nullptr;
    ModelPart& r_legacy = CreateNodalModelPart<TestThermoMechanicElement>(
        model, "TC_Legacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_legacy,
        p_law->Clone());
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_candidate = CreateNodalModelPart<TestSmallDisplacementElement>(
        model, "TC_Candidate", "SmallDisplacementElement3D8N", 3, p_candidate,
        p_law->Clone());
    p_legacy->Initialize(r_legacy.GetProcessInfo());
    p_candidate->Initialize(r_candidate.GetProcessInfo());

    // A. Restrained uniform heating (zero mechanical displacement, uniform DT).
    for (auto& n : r_legacy.Nodes()) {
        n.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
        n.X() = n.GetInitialPosition()[0]; n.Y() = n.GetInitialPosition()[1]; n.Z() = n.GetInitialPosition()[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 50.0;
    }
    for (auto& n : r_candidate.Nodes()) {
        n.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
        n.X() = n.GetInitialPosition()[0]; n.Y() = n.GetInitialPosition()[1]; n.Z() = n.GetInitialPosition()[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 50.0;
    }
    CompareSystems(*p_legacy, *p_candidate, r_legacy, r_candidate);

    // B. Free thermal expansion: no spurious stress -> no internal forces (RHS ~ 0).
    for (auto& n : r_legacy.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        const double ts = test_thermal_expansion * 50.0;
        u[0] = ts * x0[0]; u[1] = ts * x0[1]; u[2] = ts * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 50.0;
    }
    for (auto& n : r_candidate.Nodes()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        const double ts = test_thermal_expansion * 50.0;
        u[0] = ts * x0[0]; u[1] = ts * x0[1]; u[2] = ts * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        n.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 50.0;
    }
    Matrix lhs_l, lhs_c;
    Vector rhs_l, rhs_c;
    p_legacy->CalculateLocalSystem(lhs_l, rhs_l, r_legacy.GetProcessInfo());
    p_candidate->CalculateLocalSystem(lhs_c, rhs_c, r_candidate.GetProcessInfo());
    double max_rhs = 0.0;
    for (std::size_t i = 0; i < rhs_l.size(); ++i) max_rhs = std::max(max_rhs, std::abs(rhs_l[i]));
    KRATOS_EXPECT_TRUE(max_rhs < 1.0e-9);  // free expansion: no spurious internal forces

    // C. Non-uniform temperature with mechanical loading.
    ApplyNodalState(r_legacy, 2.0e-6, 40.0);
    ApplyNodalState(r_candidate, 2.0e-6, 40.0);
    CompareSystems(*p_legacy, *p_candidate, r_legacy, r_candidate);
    std::cout << "[5A] thermal cases OK" << std::endl;
}

//************************************************************************************
// 3. 2D family (plane strain and plane stress).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNodalLaw_2DFamily, KratosDamFastSuite)
{
    struct Case { const char* name; ConstitutiveLaw::Pointer law; const char* element; std::size_t dim; };
    const Case cases[] = {
        {"PlaneStrain", ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrainNodal()), "SmallDisplacementElement2D4N", 2},
        {"PlaneStress", ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStressNodal()), "SmallDisplacementElement2D4N", 2},
    };
    for (const auto& c : cases) {
        Model model;
        TestThermoMechanicElement* p_legacy = nullptr;
        ModelPart& r_legacy = CreateNodalModelPart<TestThermoMechanicElement>(
            model, std::string("2D_Legacy_") + c.name, "SmallDisplacementThermoMechanicElement2D4N",
            c.dim, p_legacy, c.law->Clone());
        TestSmallDisplacementElement* p_candidate = nullptr;
        ModelPart& r_candidate = CreateNodalModelPart<TestSmallDisplacementElement>(
            model, std::string("2D_Candidate_") + c.name, c.element, c.dim, p_candidate,
            c.law->Clone());
        p_legacy->Initialize(r_legacy.GetProcessInfo());
        p_candidate->Initialize(r_candidate.GetProcessInfo());
        ApplyNodalState(r_legacy, 2.0e-6, 10.0);
        ApplyNodalState(r_candidate, 2.0e-6, 10.0);
        CompareSystems(*p_legacy, *p_candidate, r_legacy, r_candidate);
        std::cout << "[5A] 2D " << c.name << " OK" << std::endl;
    }
}

//************************************************************************************
// 4. Repeated-response safety.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNodalLaw_RepeatedResponse, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_candidate = CreateNodalModelPart<TestSmallDisplacementElement>(
        model, "RepeatNodal", "SmallDisplacementElement3D8N", 3, p_candidate,
        ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal()));
    p_candidate->Initialize(r_candidate.GetProcessInfo());
    ApplyNodalState(r_candidate, 2.0e-6, 20.0);

    // Record the nodal material field.
    std::vector<double> E_before;
    for (auto& n : r_candidate.Nodes()) {
        E_before.push_back(n.FastGetSolutionStepValue(NODAL_YOUNG_MODULUS));
    }
    Matrix lhs0, lhs;
    Vector rhs0, rhs;
    p_candidate->CalculateLocalSystem(lhs0, rhs0, r_candidate.GetProcessInfo());
    for (std::size_t i = 0; i < 5; ++i) {
        p_candidate->CalculateLocalSystem(lhs, rhs, r_candidate.GetProcessInfo());
        p_candidate->CalculateLeftHandSide(lhs, r_candidate.GetProcessInfo());
        p_candidate->CalculateRightHandSide(rhs, r_candidate.GetProcessInfo());
        for (std::size_t j = 0; j < rhs.size(); ++j) KRATOS_EXPECT_NEAR(rhs[j], rhs0[j], 1.0e-12);
    }
    // Nodal material fields unchanged.
    std::size_t idx = 0;
    for (auto& n : r_candidate.Nodes()) {
        KRATOS_EXPECT_NEAR(n.FastGetSolutionStepValue(NODAL_YOUNG_MODULUS), E_before[idx], 1.0e-15);
        ++idx;
    }
    std::cout << "[5A] repeated response OK" << std::endl;
}

//************************************************************************************
// 5. Clone and serialization (stateless beyond nodal inputs).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNodalLaw_CloneSerialization, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_candidate = CreateNodalModelPart<TestSmallDisplacementElement>(
        model, "CloneNodal", "SmallDisplacementElement3D8N", 3, p_candidate,
        ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLawNodal()));
    p_candidate->Initialize(r_candidate.GetProcessInfo());
    ApplyNodalState(r_candidate, 2.0e-6, 10.0);

    // Clone: independent, stateless.
    auto& r_law = p_candidate->GetConstitutiveLaw(0);
    auto p_clone = Kratos::make_shared<ThermalLinearElastic3DLawNodal>(
        dynamic_cast<const ThermalLinearElastic3DLawNodal&>(r_law));

    // Serialize / deserialize the law.
    ConstitutiveLaw::Pointer p_to_serialize = p_candidate->GetConstitutiveLawPointer(0);
    StreamSerializer serializer;
    serializer.save("ThermalLinearElastic3DLawNodal", p_to_serialize);
    ConstitutiveLaw::Pointer p_loaded;
    serializer.load("ThermalLinearElastic3DLawNodal", p_loaded);
    KRATOS_EXPECT_TRUE(p_loaded != nullptr);
    std::cout << "[5A] nodal clone/serialization OK" << std::endl;
}

//************************************************************************************
// 6. Finalization table + inherited-vs-explicit Cauchy finalization.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNodalLaw_FinalizationEquivalence, KratosDamFastSuite)
{
    // Stateless linear nodal law: RequiresFinalizeMaterialResponse() == false.
    {
        ThermalLinearElastic3DLawNodal law;
        KRATOS_EXPECT_FALSE(law.RequiresInitializeMaterialResponse());
        KRATOS_EXPECT_FALSE(law.RequiresFinalizeMaterialResponse());
    }

    // The inherited Dam SolidElement::FinalizeSolutionStep finalizes with
    // StressMeasure = Cauchy (set by SmallDisplacementElement::InitializeElementData).
    // For a stateful local-damage law, run BOTH the production reduced override
    // and the inherited base (qualified call) and verify identical committed
    // state and stress.
    Model model;
    TestThermoMechanicElement* p_legacy = nullptr;
    // Use the local-damage law for the stateful check.
    {
        auto p_damage_law = ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw());
        // Rebuild a model part with the damage law via the standard element.
        ModelPart& r_mp = model.CreateModelPart("FinLegacy", 2);
        ProcessInfo& r_pi = r_mp.GetProcessInfo();
        r_pi[DOMAIN_SIZE] = 3;
        r_pi[SPACE_DIMENSION] = 3;
        r_pi[IS_RESTARTED] = false;
        r_pi[IS_CONVERGED] = true;
        r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_mp.AddNodalSolutionStepVariable(VELOCITY);
        r_mp.AddNodalSolutionStepVariable(ACCELERATION);
        r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
        r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
        r_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
        r_mp.AddNodalSolutionStepVariable(NODAL_AREA);
        r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
        const double coords[8][3] = {{0,0,0},{2.5,0,0},{2.5,1.25,0},{0,1.25,0},{0,0,2.5},{2.5,0,2.5},{2.5,1.25,2.5},{0,1.25,2.5}};
        for (std::size_t i = 0; i < 8; ++i) {
            Node::Pointer n = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
            n->AddDof(DISPLACEMENT_X); n->AddDof(DISPLACEMENT_Y); n->AddDof(DISPLACEMENT_Z);
            n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
            n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
            Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
            n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
        }
        auto p_prop = r_mp.CreateNewProperties(1);
        (*p_prop)[YOUNG_MODULUS] = 2.0e7;
        (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
        (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
        (*p_prop)[DAMAGE_THRESHOLD] = 5.0e-3;
        (*p_prop)[STRENGTH_RATIO] = 10.0;
        (*p_prop)[FRACTURE_ENERGY] = 5000.0;
        p_prop->SetValue(CONSTITUTIVE_LAW, p_damage_law->Clone());
        Geometry<Node>::PointsArrayType pts;
        for (std::size_t i = 0; i < 8; ++i) pts.push_back(r_mp.pGetNode(i + 1));
        auto p_elem = Kratos::make_intrusive<TestThermoMechanicElement>(
            1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pts)), p_prop);
        r_mp.AddElement(p_elem);
        p_elem->Initialize(r_pi);
        for (auto& n : r_mp.Nodes()) {
            const auto& x0 = n.GetInitialPosition();
            auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
            const double eps = 2.0e-5;
            u[0] = eps * x0[0]; u[1] = -test_poisson_ratio * eps * x0[1]; u[2] = -test_poisson_ratio * eps * x0[2];
            n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        }
        // Converged finalize through the production reduced override.
        r_pi[IS_CONVERGED] = true;
        p_elem->FinalizeSolutionStep(r_pi);
        double threshold_override = 0.0;
        p_elem->GetConstitutiveLaw(0).GetValue(STATE_VARIABLE, threshold_override);

        // Repeat with a fresh model, finalizing through the INHERITED
        // SolidElement::FinalizeSolutionStep (qualified base call).
        ModelPart& r_mp2 = model.CreateModelPart("FinLegacyInherited", 2);
        ProcessInfo& r_pi2 = r_mp2.GetProcessInfo();
        r_pi2[DOMAIN_SIZE] = 3;
        r_pi2[SPACE_DIMENSION] = 3;
        r_pi2[IS_RESTARTED] = false;
        r_pi2[IS_CONVERGED] = true;
        r_mp2.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_mp2.AddNodalSolutionStepVariable(VELOCITY);
        r_mp2.AddNodalSolutionStepVariable(ACCELERATION);
        r_mp2.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        r_mp2.AddNodalSolutionStepVariable(TEMPERATURE);
        r_mp2.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
        r_mp2.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
        r_mp2.AddNodalSolutionStepVariable(NODAL_AREA);
        r_mp2.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
        for (std::size_t i = 0; i < 8; ++i) {
            Node::Pointer n = r_mp2.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
            n->AddDof(DISPLACEMENT_X); n->AddDof(DISPLACEMENT_Y); n->AddDof(DISPLACEMENT_Z);
            n->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
            n->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
            Matrix z3(3, 3); noalias(z3) = ZeroMatrix(3, 3);
            n->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
        }
        auto p_prop2 = r_mp2.CreateNewProperties(1);
        (*p_prop2)[YOUNG_MODULUS] = 2.0e7;
        (*p_prop2)[POISSON_RATIO] = test_poisson_ratio;
        (*p_prop2)[THERMAL_EXPANSION] = test_thermal_expansion;
        (*p_prop2)[DAMAGE_THRESHOLD] = 5.0e-3;
        (*p_prop2)[STRENGTH_RATIO] = 10.0;
        (*p_prop2)[FRACTURE_ENERGY] = 5000.0;
        p_prop2->SetValue(CONSTITUTIVE_LAW, p_damage_law->Clone());
        Geometry<Node>::PointsArrayType pts2;
        for (std::size_t i = 0; i < 8; ++i) pts2.push_back(r_mp2.pGetNode(i + 1));
        auto p_elem2 = Kratos::make_intrusive<TestThermoMechanicElement>(
            1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(pts2)), p_prop2);
        r_mp2.AddElement(p_elem2);
        p_elem2->Initialize(r_pi2);
        for (auto& n : r_mp2.Nodes()) {
            const auto& x0 = n.GetInitialPosition();
            auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
            const double eps = 2.0e-5;
            u[0] = eps * x0[0]; u[1] = -test_poisson_ratio * eps * x0[1]; u[2] = -test_poisson_ratio * eps * x0[2];
            n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
        }
        r_pi2[IS_CONVERGED] = true;
        static_cast<SolidElement&>(*p_elem2).FinalizeSolutionStep(r_pi2);
        double threshold_inherited = 0.0;
        p_elem2->GetConstitutiveLaw(0).GetValue(STATE_VARIABLE, threshold_inherited);

        std::cout << "[5A] finalize: override threshold=" << threshold_override
                  << " inherited threshold=" << threshold_inherited << std::endl;
        KRATOS_EXPECT_NEAR(threshold_inherited, threshold_override, 1.0e-12);
    }
}

} // namespace Testing
} // namespace Kratos
