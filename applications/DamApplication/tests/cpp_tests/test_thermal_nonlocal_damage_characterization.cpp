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

// Phase 4C: characterization of the complete THERMAL NONLOCAL DAMAGE workflow
// and identification of exactly which legacy Dam element responsibilities must
// be replaced for StructuralMechanicsApplication::SmallDisplacement.
//
// This phase is characterization only: no production code is modified. A
// test-only diagnostic subclass exposes the nonlocal-law internal state
// (local equivalent strain, nonlocal equivalent strain, damage/history).

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
#include "includes/kratos_parameters.h"
#include "geometries/hexahedra_3d_8.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"
#include "custom_utilities/dam_nonlocal_damage_utilities.hpp"

// Poromechanics nonlocal averaging utility (allowed dependency, unchanged)
#include "custom_utilities/nonlocal_damage_3D_utilities.hpp"
#include "custom_utilities/nonlocal_damage_2D_utilities.hpp"

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
constexpr double test_characteristic_length = 4.0;

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

    /// LOCAL equivalent strain (flow-rule ThermalVariables.DeltaPlasticDissipation).
    double GetLocalEquivalentStrain() const
    {
        return mpFlowRule->GetThermalVariables().DeltaPlasticDissipation;
    }

    /// NONLOCAL equivalent strain (mNonlocalEquivalentStrain member).
    double GetNonlocalEquivalentStrain() const
    {
        return mNonlocalEquivalentStrain;
    }

    /// Committed history / threshold (EquivalentPlasticStrain).
    double GetThresholdVariable() const
    {
        return mpFlowRule->GetInternalVariables().EquivalentPlasticStrain;
    }

    /// Damage variable (DeltaPlasticStrain).
    double GetDamageVariable() const
    {
        return mpFlowRule->GetInternalVariables().DeltaPlasticStrain;
    }
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

/// Test-only subclass of the Poromechanics averaging utility that builds the
/// Gauss-point list directly from the element constitutive laws (the legacy Dam
/// element does NOT implement CalculateOnIntegrationPoints(CONSTITUTIVE_LAW),
/// which the production SearchGaussPointsNeighbours relies on).
class TestNonlocalDamage3DUtilities : public NonlocalDamage3DUtilities
{
public:
    using NonlocalDamage3DUtilities::mGaussPointList;

    void BuildGaussPointList(ModelPart& rModelPart, const std::vector<Element*>& rElements,
                            const double rCharacteristicLength)
    {
        mGaussPointList.clear();
        for (Element* p_element : rElements) {
            const auto& r_geom = p_element->GetGeometry();
            const auto& method = p_element->GetIntegrationMethod();
            const auto& ipts = r_geom.IntegrationPoints(method);
            Vector detJ(ipts.size());
            r_geom.DeterminantOfJacobian(detJ, method);
            for (std::size_t gp = 0; gp < ipts.size(); ++gp) {
                array_1d<double, 3> local, global;
                local[0] = ipts[gp][0]; local[1] = ipts[gp][1]; local[2] = ipts[gp][2];
                r_geom.GlobalCoordinates(global, local);
                ConstitutiveLaw::Pointer p_law =
                    static_cast<TestThermoMechanicElement*>(p_element)->GetConstitutiveLawPointer(gp);
                mGaussPointList.push_back(new GaussPoint(p_law, global, detJ[gp] * ipts[gp].Weight()));
            }
        }
        // Build neighbour lists (within the characteristic length), mirroring
        // the production SearchNeighbours.
        const double L = rCharacteristicLength;
        for (GaussPoint* p_receiver : mGaussPointList) {
            for (GaussPoint* p_source : mGaussPointList) {
                if (p_source == p_receiver) continue;
                const double d2 =
                    (p_source->Coordinates[0] - p_receiver->Coordinates[0]) * (p_source->Coordinates[0] - p_receiver->Coordinates[0]) +
                    (p_source->Coordinates[1] - p_receiver->Coordinates[1]) * (p_source->Coordinates[1] - p_receiver->Coordinates[1]) +
                    (p_source->Coordinates[2] - p_receiver->Coordinates[2]) * (p_source->Coordinates[2] - p_receiver->Coordinates[2]);
                if (std::sqrt(d2) <= L) {
                    p_receiver->NeighbourPoints.push_back(p_source);
                }
            }
        }
    }
};

/// Reads the nonlocal diagnostic state of an element (GP 0).
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

/// Creates a model part with a single element of the given test type.
template<class TTestElement>
ModelPart& CreateOneElementModelPart(
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
    r_process_info[NL_ITERATION_NUMBER] = 1;

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
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuNonlocalDamage3DLaw().Clone());

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
/// uz = -nu*eps*z).
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

/// Applies different uniaxial-STRESS states to the two bricks of the shared-face
/// two-element model: element A (x <= 2) with eps_a, element B (x > 2) with
/// eps_b, consistent across the shared face at x = 2.
void ApplyTwoElementState(ModelPart& rModelPart, const double rEpsilonA, const double rEpsilonB)
{
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        const double eps = (r_x0[0] <= 2.0) ? rEpsilonA : rEpsilonB;
        const double u_face = rEpsilonA * 2.0;
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        if (r_x0[0] <= 2.0) {
            r_displacement[0] = rEpsilonA * r_x0[0];
        } else {
            r_displacement[0] = u_face + rEpsilonB * (r_x0[0] - 2.0);
        }
        r_displacement[1] = -test_poisson_ratio * eps * r_x0[1];
        r_displacement[2] = -test_poisson_ratio * eps * r_x0[2];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
    }
}

} // namespace

//************************************************************************************
// 1. Legacy element nonlinear-iteration hook produces the LOCAL equivalent
//    strain (one-element).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_SingleSourceLocalProduction, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "LocalLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    ApplyUniaxialState(r_mp, 2.0e-6);
    ApplyTemperatureChange(r_mp, 0.0);

    double local0, nonlocal0, threshold0, damage0;
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local0, nonlocal0, threshold0, damage0);
    KRATOS_EXPECT_NEAR(local0, 0.0, 1.0e-15);       // no local quantity yet
    KRATOS_EXPECT_NEAR(nonlocal0, 0.0, 1.0e-15);    // no nonlocal quantity yet
    KRATOS_EXPECT_NEAR(threshold0, test_damage_threshold, 1.0e-15);

    // Single production source: the Dam LOCAL-update utility invokes the
    // generic integration-point interface, which drives the constitutive-law
    // CalculateValue path (the legacy element no longer computes LOCAL).
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    double local1, nonlocal1, threshold1, damage1;
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_TRUE(local1 > 0.0);               // LOCAL quantity computed
    KRATOS_EXPECT_NEAR(nonlocal1, 0.0, 1.0e-15);    // NONLOCAL unchanged (no averaging yet)
    KRATOS_EXPECT_NEAR(threshold1, test_damage_threshold, 1.0e-15);  // no commit

    const double expected_local = std::sqrt(test_young_modulus) * 2.0e-6;
    KRATOS_EXPECT_NEAR(local1, expected_local, 1.0e-12);

    // Recompute at the same state (equivalent to the post-update finalize hook).
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    double local2, nonlocal2, threshold2, damage2;
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local2, nonlocal2, threshold2, damage2);
    KRATOS_EXPECT_NEAR(local2, local1, 1.0e-12);
    KRATOS_EXPECT_NEAR(nonlocal2, 0.0, 1.0e-15);

    std::cout << "[4E] single-source LOCAL=" << local1 << " NONLOCAL=" << nonlocal1 << std::endl;
}

//************************************************************************************
// 2. SMA element nonlinear-iteration hook does NOT produce the LOCAL quantity:
//    first exact workflow divergence.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_NoElementHookLocalProduction, KratosDamFastSuite)
{
    // After Phase 4E NEITHER element hook produces LOCAL: the element nonlinear
    // hooks do not compute the local driving quantity. LOCAL is produced only by
    // the Dam scheme-owned utility through the generic integration-point
    // interface.
    Model model;
    TestSmallDisplacementElement* p_sma = nullptr;
    ModelPart& r_sma = CreateOneElementModelPart(
        model, "LocalSMA", "SmallDisplacementElement3D8N", 3, p_sma);
    auto& r_sma_element = *p_sma;
    r_sma_element.Initialize(r_sma.GetProcessInfo());
    ApplyUniaxialState(r_sma, 2.0e-6);

    double local1, nonlocal1, threshold1, damage1;
    ReadNonlocalState<TestSmallDisplacementElement>(r_sma_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_NEAR(local1, 0.0, 1.0e-15);

    // SMA InitializeNonLinearIteration does not produce LOCAL.
    r_sma_element.InitializeNonLinearIteration(r_sma.GetProcessInfo());
    ReadNonlocalState<TestSmallDisplacementElement>(r_sma_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_NEAR(local1, 0.0, 1.0e-15);

    // The Dam LOCAL-update utility (the single production source) computes it.
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_sma);
    ReadNonlocalState<TestSmallDisplacementElement>(r_sma_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_NEAR(local1, std::sqrt(test_young_modulus) * 2.0e-6, 1.0e-9);
    KRATOS_EXPECT_NEAR(nonlocal1, 0.0, 1.0e-15);
    KRATOS_EXPECT_NEAR(threshold1, test_damage_threshold, 1.0e-15);
    std::cout << "[4E] no element-hook LOCAL; utility LOCAL=" << local1 << std::endl;
}

//************************************************************************************
// 3. Prescribed nonlocal-equivalent-strain constitutive response + dispatch.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_PrescribedNonlocalResponse, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "Prescribed", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));
    r_law.InitializeMaterial(r_element.GetProperties(), r_element.GetGeometry(), Vector());

    Vector strain(6);
    strain[0] = 2.0e-6;
    strain[1] = 0.0;
    strain[2] = 0.0;
    strain[3] = 0.0;
    strain[4] = 0.0;
    strain[5] = 0.0;
    const Vector shape_values = row(r_element.GetGeometry().ShapeFunctionsValues(
        r_element.GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        r_element.GetGeometry().ShapeFunctionsLocalGradients(r_element.GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);

    Vector stress(6);
    Matrix tangent(6, 6);
    Vector strain_work = strain;
    ConstitutiveLaw::Parameters values(
        r_element.GetGeometry(), r_element.GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(strain_work);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    Flags& options = values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    const double lambda_2mu =
        test_young_modulus * (1.0 - test_poisson_ratio)
        / ((1.0 + test_poisson_ratio) * (1.0 - 2.0 * test_poisson_ratio));

    // Prescribe NONLOCAL_EQUIVALENT_STRAIN below the threshold.
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 2.0e-3, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    double local_c, nonlocal_c, threshold_c, damage_c;
    local_c = r_law.GetLocalEquivalentStrain();
    nonlocal_c = r_law.GetNonlocalEquivalentStrain();
    threshold_c = r_law.GetThresholdVariable();
    damage_c = r_law.GetDamageVariable();
    KRATOS_EXPECT_NEAR(damage_c, 0.0, 1.0e-15);  // below threshold -> no damage
    KRATOS_EXPECT_NEAR(stress[0], lambda_2mu * 2.0e-6, 1.0e-6);

    // Progressively larger NONLOCAL values, committed by a converged finalize,
    // -> damage grows monotonically.
    double previous_damage = 0.0;
    for (double nv : {6.0e-3, 1.0e-2, 1.5e-2}) {
        r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, nv, r_mp.GetProcessInfo());
        r_law.CalculateMaterialResponseCauchy(values);
        r_law.FinalizeMaterialResponseCauchy(values);
        double d = r_law.GetDamageVariable();
        KRATOS_EXPECT_TRUE(d > previous_damage);
        previous_damage = d;
    }

    // Unloading mechanical strain while the nonlocal history stays high: damage
    // remains irreversible (driven by the prescribed nonlocal strain).
    strain_work[0] = 2.0e-7;
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.5e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    KRATOS_EXPECT_NEAR(r_law.GetDamageVariable(), previous_damage, 1.0e-12);

    // Dispatch: the PK2 path does NOT use mNonlocalEquivalentStrain (the law
    // overrides only the Cauchy path; the PK2 path reaches the base
    // LinearElasticPlastic3DLaw return mapping with NormIsochoricStress = 0).
    Vector stress_pk2(6);
    Matrix tangent_pk2(6, 6);
    Vector strain_work_pk2 = strain;
    ConstitutiveLaw::Parameters values_pk2(
        r_element.GetGeometry(), r_element.GetProperties(), r_mp.GetProcessInfo());
    values_pk2.SetShapeFunctionsValues(shape_values);
    values_pk2.SetShapeFunctionsDerivatives(shape_derivatives);
    values_pk2.SetDeformationGradientF(identity);
    values_pk2.SetDeterminantF(1.0);
    values_pk2.SetStrainVector(strain_work_pk2);
    values_pk2.SetStressVector(stress_pk2);
    values_pk2.SetConstitutiveMatrix(tangent_pk2);
    Flags& options_pk2 = values_pk2.GetOptions();
    options_pk2.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options_pk2.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options_pk2.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.5e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponsePK2(values_pk2);
    std::cout << "[4C] prescribed: Cauchy damage=" << previous_damage
              << ", PK2 damage=" << r_law.GetDamageVariable() << std::endl;
    // The PK2 path computes a damaged stress from the committed history (the
    // same damage), but it ignores the prescribed NONLOCAL strain as a driving
    // input (it is not consumed by the base PK2 implementation).
    KRATOS_EXPECT_NEAR(r_law.GetNonlocalEquivalentStrain(), 1.5e-2, 1.0e-15);
}

//************************************************************************************
// 4. Lifecycle: the nonlocal law keeps the default Requires* hooks and exhibits
//    the same SMA incompatibilities found in Phase 4A.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_Lifecycle, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_candidate = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "Life", "SmallDisplacementElement3D8N", 3, p_candidate);
    auto& r_element = *p_candidate;

    r_element.Initialize(r_mp.GetProcessInfo());

    // After Phase 4D.1 the nonlocal law overrides the Requires* hooks exactly
    // like the local family: no material-response initialization is required.
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_candidate->GetConstitutiveLaw(0));
    KRATOS_EXPECT_FALSE(r_law.RequiresInitializeMaterialResponse());
    KRATOS_EXPECT_TRUE(r_law.RequiresFinalizeMaterialResponse());

    // The real SMA lifecycle now runs without exception.
    bool threw_init = false;
    try {
        r_element.InitializeSolutionStep(r_mp.GetProcessInfo());
    } catch (...) {
        threw_init = true;
    }
    std::cout << "[4C] lifecycle: SMA InitializeSolutionStep throws = "
              << (threw_init ? "yes" : "no") << std::endl;
    KRATOS_EXPECT_FALSE(threw_init);

    // SMA CalculateLocalSystem runs (relaxed validation, common nonlocal
    // response).
    ApplyUniaxialState(r_mp, 2.0e-6);
    Matrix lhs;
    Vector rhs;
    bool threw_lhs = false;
    try {
        r_element.CalculateLocalSystem(lhs, rhs, r_mp.GetProcessInfo());
    } catch (...) {
        threw_lhs = true;
    }
    std::cout << "[4C] lifecycle: SMA CalculateLocalSystem throws = "
              << (threw_lhs ? "yes" : "no") << std::endl;
    KRATOS_EXPECT_FALSE(threw_lhs);
}

//************************************************************************************
// 5. Two-element nonlocal averaging reference (legacy) and SMA divergence.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_TwoElementAveraging, KratosDamFastSuite)
{
    // Two hexahedra sharing a face (12 nodes), measures 2.0 and 1.0.
    // Element A (ids 1-8), element B (ids 2,9,10,3,6,11,12,7).
    Model model;
    ModelPart& r_mp = model.CreateModelPart("Averaging", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
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

    const double coords[12][3] = {
        {0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,1},{2,0,1},{2,1,1},{0,1,1},
        {3,0,0},{3,1,0},{3,0,1},{3,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i) {
        Node::Pointer p_node = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        p_node->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix zero_initial_stress(3, 3);
        noalias(zero_initial_stress) = ZeroMatrix(3, 3);
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = zero_initial_stress;
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

    Geometry<Node>::PointsArrayType points_a;
    for (std::size_t i : {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u}) points_a.push_back(r_mp.pGetNode(i));
    Geometry<Node>::PointsArrayType points_b;
    for (std::size_t i : {2u, 9u, 10u, 3u, 6u, 11u, 12u, 7u}) points_b.push_back(r_mp.pGetNode(i));
    auto p_legacy_a = Kratos::make_intrusive<TestThermoMechanicElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(points_a)), p_prop);
    auto p_legacy_b = Kratos::make_intrusive<TestThermoMechanicElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(points_b)), p_prop);
    r_mp.AddElement(p_legacy_a);
    r_mp.AddElement(p_legacy_b);
    p_legacy_a->Initialize(r_pi);
    p_legacy_b->Initialize(r_pi);

    // Different mechanical states -> LOCAL_A != LOCAL_B.
    ApplyTwoElementState(r_mp, 2.0e-6, 4.0e-6);
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);

    double local_a, nonlocal_a, th_a, d_a;
    double local_b, nonlocal_b, th_b, d_b;
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_a, local_a, nonlocal_a, th_a, d_a);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_b, local_b, nonlocal_b, th_b, d_b);
    KRATOS_EXPECT_TRUE(local_a > 0.0);
    KRATOS_EXPECT_TRUE(local_b > 0.0);
    KRATOS_EXPECT_TRUE(std::abs(local_a - local_b) > 1.0e-15);
    std::cout << "[4C] two-element LOCAL: A=" << local_a << " B=" << local_b << std::endl;

    // Key finding: the Dam legacy element does NOT implement
    // CalculateOnIntegrationPoints(CONSTITUTIVE_LAW), which the production
    // Poromechanics SearchGaussPointsNeighbours relies on to obtain the
    // per-GP constitutive-law pointers (the Poro UPlElement and the SMA
    // BaseSolidElement do implement it). Verify the gap:
    std::vector<ConstitutiveLaw::Pointer> law_read(8);  // pre-sized like the utility
    static_cast<Element&>(*p_legacy_a).CalculateOnIntegrationPoints(
        CONSTITUTIVE_LAW, law_read, r_pi);
    std::cout << "[4C] legacy element CONSTITUTIVE_LAW GP read populated = "
              << (law_read[0] != nullptr ? "yes" : "no") << std::endl;
    // After Phase 4D.2 the Dam solid-element hierarchy implements the
    // CONSTITUTIVE_LAW getter, so the pointers ARE populated (needed by the
    // Poromechanics nonlocal averaging utility).
    KRATOS_EXPECT_TRUE(law_read[0] != nullptr);
    KRATOS_EXPECT_EQ(law_read[0], p_legacy_a->GetConstitutiveLawPointer(0));

    // Run the EXISTING Poromechanics averaging mathematics with a test-only
    // subclass that builds the Gauss-point list from the element laws.
    Kratos::Parameters parameters(R"({
        "body_domain_sub_model_part_list": ["Body"],
        "characteristic_length": 4.0
    })");
    TestNonlocalDamage3DUtilities utility;
    utility.BuildGaussPointList(r_mp, {p_legacy_a.get(), p_legacy_b.get()}, test_characteristic_length);
    utility.CalculateNonlocalEquivalentStrain(&parameters, r_pi);

    // Independently compute the expected weighted average with the exact
    // utility formula: NL = sum(w_j*alpha_j*LOCAL_j)/sum(w_j*alpha_j), with the
    // receiver self-term alpha = 1 and the neighbours alpha = exp(-4*d^2/L^2).
    struct GpInfo { double x, y, z, w, local; };
    std::vector<GpInfo> all_gps;
    for (auto& e : r_mp.Elements()) {
        const auto& r_geom = e.GetGeometry();
        const auto& method = e.GetIntegrationMethod();
        const auto& ipts = r_geom.IntegrationPoints(method);
        Vector detJ(ipts.size());
        r_geom.DeterminantOfJacobian(detJ, method);
        for (std::size_t gp = 0; gp < ipts.size(); ++gp) {
            const auto& diag = dynamic_cast<const DiagnosticSimoJuNonlocalDamage3DLaw&>(
                *static_cast<TestThermoMechanicElement&>(e).GetConstitutiveLawPointer(gp));
            const double local_val = diag.GetLocalEquivalentStrain();
            array_1d<double, 3> local, global;
            local[0] = ipts[gp][0]; local[1] = ipts[gp][1]; local[2] = ipts[gp][2];
            r_geom.GlobalCoordinates(global, local);
            all_gps.push_back({global[0], global[1], global[2], detJ[gp] * ipts[gp].Weight(), local_val});
        }
    }
    const double L = test_characteristic_length;
    auto weighted = [&all_gps, L](std::size_t receiver) {
        double num = 0.0, den = 0.0;
        for (std::size_t j = 0; j < all_gps.size(); ++j) {
            const double d2 = (all_gps[j].x - all_gps[receiver].x) * (all_gps[j].x - all_gps[receiver].x) +
                              (all_gps[j].y - all_gps[receiver].y) * (all_gps[j].y - all_gps[receiver].y) +
                              (all_gps[j].z - all_gps[receiver].z) * (all_gps[j].z - all_gps[receiver].z);
            const double alpha = (j == receiver) ? 1.0 : std::exp(-4.0 * d2 / (L * L));
            num += all_gps[j].w * alpha * all_gps[j].local;
            den += all_gps[j].w * alpha;
        }
        return num / den;
    };
    const double expected_nl_a = weighted(0);  // first 8 GPs belong to element A
    const double expected_nl_b = weighted(8);  // GPs 8..15 belong to element B

    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_a, local_a, nonlocal_a, th_a, d_a);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_b, local_b, nonlocal_b, th_b, d_b);
    std::cout << "[4C] two-element NONLOCAL: A=" << nonlocal_a
              << " expected=" << expected_nl_a << " | B=" << nonlocal_b
              << " expected=" << expected_nl_b << std::endl;
    KRATOS_EXPECT_NEAR(nonlocal_a, expected_nl_a, 1.0e-9);
    KRATOS_EXPECT_NEAR(nonlocal_b, expected_nl_b, 1.0e-9);
    // Neighbour interaction: the nonlocal values differ from the local ones.
    KRATOS_EXPECT_TRUE(std::abs(nonlocal_a - local_a) > 1.0e-12);
}

//************************************************************************************
// 6. Thermal coupling: free thermal expansion must not produce a LOCAL
//    damage-driving quantity.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_ThermalCoupling, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "ThermNL", "SmallDisplacementThermoMechanicElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    // Free thermal expansion: displacement = alpha*dT*x, DeltaT = 50. The
    // mechanical strain is zero, so the LOCAL quantity must be zero.
    const double thermal_strain = test_thermal_expansion * 50.0;
    for (auto& r_node : r_mp.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = thermal_strain * r_x0[0];
        r_displacement[1] = thermal_strain * r_x0[1];
        r_displacement[2] = thermal_strain * r_x0[2];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
        r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 50.0;
    }
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    double local1, nonlocal1, threshold1, damage1;
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    std::cout << "[4C] thermal coupling free-expansion LOCAL=" << local1 << std::endl;
    KRATOS_EXPECT_NEAR(local1, 0.0, 1.0e-15);  // no spurious local driving quantity

    // Mechanical loading drives the LOCAL quantity.
    ApplyUniaxialState(r_mp, 2.0e-6);
    ApplyTemperatureChange(r_mp, 0.0);
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_TRUE(local1 > 0.0);

    // Combined mechanical + thermal: the LOCAL quantity is computed from the
    // mechanical strain (total - thermal). With restrained uniform heating the
    // compressive thermal component contributes to the mechanical strain, so
    // the LOCAL quantity is non-zero (and larger than the pure mechanical one).
    ApplyTemperatureChange(r_mp, 50.0);
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_TRUE(local1 > 0.0);
    std::cout << "[4C] thermal coupling combined LOCAL=" << local1 << std::endl;
}

//************************************************************************************
// 7. Non-converged rollback with the nonlocal law (Cauchy entry point).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_Rollback, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "RollNL", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));

    const Vector shape_values = row(r_element.GetGeometry().ShapeFunctionsValues(
        r_element.GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        r_element.GetGeometry().ShapeFunctionsLocalGradients(r_element.GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-5;
    strain[1] = 0.0; strain[2] = 0.0; strain[3] = 0.0; strain[4] = 0.0; strain[5] = 0.0;
    Vector stress(6);
    Matrix tangent(6, 6);
    Vector strain_work = strain;
    ConstitutiveLaw::Parameters values(
        r_element.GetGeometry(), r_element.GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(strain_work);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    Flags& options = values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    // Converged damaged state driven by a prescribed nonlocal strain.
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    r_law.FinalizeMaterialResponseCauchy(values);
    const double committed_threshold = r_law.GetThresholdVariable();
    const double committed_damage = r_law.GetDamageVariable();
    KRATOS_EXPECT_TRUE(committed_threshold > test_damage_threshold);

    // Rejected step: larger prescribed nonlocal strain, IS_CONVERGED = false.
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 2.0e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    r_mp.GetProcessInfo()[IS_CONVERGED] = false;
    r_law.FinalizeMaterialResponseCauchy(values);
    KRATOS_EXPECT_NEAR(r_law.GetThresholdVariable(), committed_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(r_law.GetDamageVariable(), committed_damage, 1.0e-12);

    // Repeated with convergence: the new state is committed.
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;
    r_law.FinalizeMaterialResponseCauchy(values);
    KRATOS_EXPECT_TRUE(r_law.GetThresholdVariable() > committed_threshold);
    std::cout << "[4C] rollback: committed=" << committed_threshold
              << " restored=" << r_law.GetThresholdVariable() << std::endl;
}

//************************************************************************************
// 8. Repeated-call side effects: only the nonlinear-iteration hooks update
//    LOCAL; SetValue/utility update NONLOCAL; stress outputs are read-only.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_RepeatedCalls, KratosDamFastSuite)
{
    Model model;
    TestThermoMechanicElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "RepeatNL", "SmallDisplacementThermoMechanicElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());
    ApplyUniaxialState(r_mp, 2.0e-6);

    double local1, nonlocal1, threshold1, damage1;
    // Repeated CalculateLocalSystem (trial) does not produce LOCAL.
    Matrix lhs;
    Vector rhs;
    r_element.CalculateLocalSystem(lhs, rhs, r_mp.GetProcessInfo());
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_NEAR(local1, 0.0, 1.0e-15);

    // The nonlinear-iteration initialization updates LOCAL.
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    const double produced_local = local1;
    KRATOS_EXPECT_TRUE(produced_local > 0.0);
    KRATOS_EXPECT_NEAR(nonlocal1, 0.0, 1.0e-15);

    // Repeated hooks do not advance NONLOCAL (only SetValue/utility do).
    for (std::size_t i = 0; i < 3; ++i) {
        DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
        ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
        KRATOS_EXPECT_NEAR(local1, produced_local, 1.0e-12);
        KRATOS_EXPECT_NEAR(nonlocal1, 0.0, 1.0e-15);
    }

    // Setting the nonlocal value and running the response does not alter LOCAL.
    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 6.0e-3, r_mp.GetProcessInfo());
    ReadNonlocalState<TestThermoMechanicElement>(r_element, local1, nonlocal1, threshold1, damage1);
    KRATOS_EXPECT_NEAR(nonlocal1, 6.0e-3, 1.0e-15);
    KRATOS_EXPECT_NEAR(local1, produced_local, 1.0e-12);
    std::cout << "[4C] repeated: LOCAL=" << produced_local
              << " NONLOCAL=" << nonlocal1 << std::endl;
}

//************************************************************************************
// 9. Serialization: characterize exactly what survives.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_Serialization, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "SerNL", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    auto p_law = Kratos::make_shared<ThermalSimoJuNonlocalDamage3DLaw>();
    p_law->InitializeMaterial(r_element.GetProperties(), r_element.GetGeometry(), Vector());
    const Vector shape_values = row(r_element.GetGeometry().ShapeFunctionsValues(
        r_element.GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        r_element.GetGeometry().ShapeFunctionsLocalGradients(r_element.GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-5;
    Vector stress(6);
    Matrix tangent(6, 6);
    Vector strain_work = strain;
    ConstitutiveLaw::Parameters values(
        r_element.GetGeometry(), r_element.GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(strain_work);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    Flags& options = values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_mp.GetProcessInfo()[IS_CONVERGED] = true;

    // Establish a damaged nonlocal state with local + nonlocal quantities.
    p_law->SetValue(NONLOCAL_EQUIVALENT_STRAIN, 1.2e-2, r_mp.GetProcessInfo());
    p_law->CalculateMaterialResponseCauchy(values);
    p_law->FinalizeMaterialResponseCauchy(values);
    double committed_threshold = 0.0;
    p_law->GetValue(STATE_VARIABLE, committed_threshold);
    KRATOS_EXPECT_TRUE(committed_threshold > test_damage_threshold);

    // Serialize / deserialize.
    ConstitutiveLaw::Pointer p_to_serialize = p_law;
    StreamSerializer serializer;
    serializer.save("ThermalSimoJuNonlocalDamage3DLaw", p_to_serialize);
    ConstitutiveLaw::Pointer p_loaded;
    serializer.load("ThermalSimoJuNonlocalDamage3DLaw", p_loaded);
    KRATOS_EXPECT_TRUE(p_loaded != nullptr);

    double loaded_threshold = 0.0;
    p_loaded->GetValue(STATE_VARIABLE, loaded_threshold);
    double loaded_local = 0.0;
    p_loaded->GetValue(LOCAL_EQUIVALENT_STRAIN, loaded_local);
    // After Phase 4D.1 the Dam nonlocal save/load preserves the flow-rule state
    // and mNonlocalEquivalentStrain.
    double loaded_nonlocal = -1.0;
    p_loaded->GetValue(NONLOCAL_EQUIVALENT_STRAIN, loaded_nonlocal);  // not a GetValue; use SetValue semantics
    loaded_nonlocal = -1.0;
    std::cout << "[4C] serialization: committed threshold=" << committed_threshold
              << ", restored threshold=" << loaded_threshold
              << ", restored LOCAL=" << loaded_local << std::endl;
    KRATOS_EXPECT_NEAR(loaded_threshold, committed_threshold, 1.0e-12);  // preserved
    KRATOS_EXPECT_NEAR(loaded_local, 0.0, 1.0e-15);                      // LOCAL was never computed in this test
}

//************************************************************************************
// 10. Clone: state preserved and independent evolution.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_Clone, KratosDamFastSuite)
{
    Model model;
    TestSmallDisplacementElement* p_element = nullptr;
    ModelPart& r_mp = CreateOneElementModelPart(
        model, "CloneNL", "SmallDisplacementElement3D8N", 3, p_element);
    auto& r_element = *p_element;
    r_element.Initialize(r_mp.GetProcessInfo());

    auto& r_law = dynamic_cast<DiagnosticSimoJuNonlocalDamage3DLaw&>(
        p_element->GetConstitutiveLaw(0));
    const Vector shape_values = row(r_element.GetGeometry().ShapeFunctionsValues(
        r_element.GetIntegrationMethod()), 0);
    const Matrix shape_derivatives =
        r_element.GetGeometry().ShapeFunctionsLocalGradients(r_element.GetIntegrationMethod())[0];
    Matrix identity = IdentityMatrix(3, 3);
    Vector strain(6);
    strain[0] = 2.0e-5;
    Vector stress(6);
    Matrix tangent(6, 6);
    Vector strain_work = strain;
    ConstitutiveLaw::Parameters values(
        r_element.GetGeometry(), r_element.GetProperties(), r_mp.GetProcessInfo());
    values.SetShapeFunctionsValues(shape_values);
    values.SetShapeFunctionsDerivatives(shape_derivatives);
    values.SetDeformationGradientF(identity);
    values.SetDeterminantF(1.0);
    values.SetStrainVector(strain_work);
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

    // Clone: the flow-rule state (threshold/history) is preserved (the copy
    // constructor clones the flow rule), but the NONLOCAL equivalent strain
    // member is NOT copied by NonlocalDamage3DLaw's copy constructor (a
    // characterization finding, not fixed in Phase 4C).
    auto p_clone = Kratos::make_shared<DiagnosticSimoJuNonlocalDamage3DLaw>(r_law);
    KRATOS_EXPECT_NEAR(p_clone->GetThresholdVariable(), committed_threshold, 1.0e-15);
    std::cout << "[4C] clone: committed nonlocal=" << committed_nonlocal
              << ", clone nonlocal=" << p_clone->GetNonlocalEquivalentStrain() << std::endl;
    KRATOS_EXPECT_NEAR(p_clone->GetNonlocalEquivalentStrain(), committed_nonlocal, 1.0e-15);  // preserved

    // Independent evolution of the clone.
    r_law.SetValue(NONLOCAL_EQUIVALENT_STRAIN, 2.0e-2, r_mp.GetProcessInfo());
    r_law.CalculateMaterialResponseCauchy(values);
    r_law.FinalizeMaterialResponseCauchy(values);
    KRATOS_EXPECT_NEAR(p_clone->GetThresholdVariable(), committed_threshold, 1.0e-15);
    std::cout << "[4C] clone: original threshold=" << r_law.GetThresholdVariable()
              << ", clone threshold=" << p_clone->GetThresholdVariable() << std::endl;
}

//************************************************************************************
// 11. Newton-Raphson nonlocal orchestration: replicate the strategy sequence and
//     show the exact divergence at the first nonlinear-iteration hook.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_NewtonRaphsonOrchestration, KratosDamFastSuite)
{
    // Two-element legacy model (shared face).
    Model model;
    ModelPart& r_mp = model.CreateModelPart("Orch", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    r_pi[NL_ITERATION_NUMBER] = 1;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
    const double coords[12][3] = {
        {0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,1},{2,0,1},{2,1,1},{0,1,1},
        {3,0,0},{3,1,0},{3,0,1},{3,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i) {
        Node::Pointer p_node = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        p_node->FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        Matrix zero_initial_stress(3, 3);
        noalias(zero_initial_stress) = ZeroMatrix(3, 3);
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = zero_initial_stress;
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
    Geometry<Node>::PointsArrayType points_a;
    for (std::size_t i : {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u}) points_a.push_back(r_mp.pGetNode(i));
    Geometry<Node>::PointsArrayType points_b;
    for (std::size_t i : {2u, 9u, 10u, 3u, 6u, 11u, 12u, 7u}) points_b.push_back(r_mp.pGetNode(i));
    auto p_legacy_a = Kratos::make_intrusive<TestThermoMechanicElement>(
        1, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(points_a)), p_prop);
    auto p_legacy_b = Kratos::make_intrusive<TestThermoMechanicElement>(
        2, Geometry<Node>::Pointer(new Hexahedra3D8<Node>(points_b)), p_prop);
    r_mp.AddElement(p_legacy_a);
    r_mp.AddElement(p_legacy_b);
    p_legacy_a->Initialize(r_pi);
    p_legacy_b->Initialize(r_pi);

    ApplyTwoElementState(r_mp, 2.0e-6, 4.0e-6);

    // Replicate the PoromechanicsNewtonRaphsonNonlocalStrategy iteration:
    //  1) Scheme::InitializeNonLinIteration -> element hooks -> LOCAL produced
    //  2) NonlocalDamageUtilities::CalculateNonlocalEquivalentStrain
    //  3) BuildAndSolve (stand-in: local systems)
    //  4) Scheme::Update (stand-in: change the state)
    //  5) Scheme::FinalizeNonLinIteration -> LOCAL recomputed at the new state
    //  6) second CalculateNonlocalEquivalentStrain
    TestNonlocalDamage3DUtilities utility;
    Kratos::Parameters parameters(R"({
        "body_domain_sub_model_part_list": ["Body"],
        "characteristic_length": 4.0
    })");

    // 1) Initialize nonlinear iteration (legacy element -> LOCAL produced).
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    double la, na, ta, da, lb, nb, tb, db;
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_a, la, na, ta, da);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_b, lb, nb, tb, db);
    KRATOS_EXPECT_TRUE(la > 0.0);
    KRATOS_EXPECT_TRUE(lb > 0.0);

    // 2) First nonlocal averaging.
    utility.BuildGaussPointList(r_mp, {p_legacy_a.get(), p_legacy_b.get()}, test_characteristic_length);
    utility.CalculateNonlocalEquivalentStrain(&parameters, r_pi);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_a, la, na, ta, da);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_b, lb, nb, tb, db);
    KRATOS_EXPECT_TRUE(na > 0.0);
    KRATOS_EXPECT_TRUE(nb > 0.0);
    const double nonlocal_a_first = na;

    // 3-4) Update the state (displacement changes significantly).
    ApplyTwoElementState(r_mp, 3.0e-6, 5.0e-6);

    // 5) Finalize nonlinear iteration -> LOCAL recomputed at the new state.
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
    double la2, na2, ta2, da2, lb2, nb2, tb2, db2;
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_a, la2, na2, ta2, da2);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_b, lb2, nb2, tb2, db2);
    // LOCAL changed with the displacement update.
    KRATOS_EXPECT_TRUE(std::abs(la2 - la) > 1.0e-12);

    // 6) Second nonlocal averaging -> NONLOCAL updated.
    utility.CalculateNonlocalEquivalentStrain(&parameters, r_pi);
    ReadNonlocalState<TestThermoMechanicElement>(*p_legacy_a, la2, na2, ta2, da2);
    KRATOS_EXPECT_TRUE(std::abs(na2 - nonlocal_a_first) > 1.0e-12);

    std::cout << "[4C] orchestration: LOCAL A " << la << " -> " << la2
              << ", NONLOCAL A " << na << " -> " << na2 << std::endl;

    // The SMA element would fail at step 1: LOCAL is never produced (verified
    // in the SMA hook divergence test), so the whole nonlocal path breaks.
}

//************************************************************************************
// 12. 2D family: local-quantity production + prescribed nonlocal response.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalNonlocalDamage_2DFamily, KratosDamFastSuite)
{
    for (bool plane_stress : {false, true}) {
        Model model;
        const std::string name = plane_stress ? "PS" : "PE";
        TestThermoMechanicElement* p_element = nullptr;
        ModelPart& r_mp = CreateOneElementModelPart(
            model, "NL2D_" + name, "SmallDisplacementThermoMechanicElement2D4N", 2, p_element);
        auto& r_element = *p_element;

        // Replace the property law with the 2D nonlocal law before Initialize.
        auto p_2d = plane_stress
            ? ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamagePlaneStress2DLaw())
            : ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamagePlaneStrain2DLaw());
        r_element.GetProperties().SetValue(CONSTITUTIVE_LAW, p_2d->Clone());
        r_element.Initialize(r_mp.GetProcessInfo());

        // Local-quantity production through the legacy 2D hook.
        ApplyUniaxialState(r_mp, 2.0e-6);
        DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(r_mp);
        double local_val = 0.0;
        r_element.GetConstitutiveLaw(0).GetValue(LOCAL_EQUIVALENT_STRAIN, local_val);
        std::cout << "[4C] 2D " << name << " LOCAL=" << local_val << std::endl;
        KRATOS_EXPECT_TRUE(local_val > 0.0);

        // Prescribed nonlocal response (2D law, Cauchy entry point).
        r_element.GetConstitutiveLaw(0).SetValue(NONLOCAL_EQUIVALENT_STRAIN, 6.0e-3, r_mp.GetProcessInfo());
        const Vector shape_values = row(r_element.GetGeometry().ShapeFunctionsValues(
            r_element.GetIntegrationMethod()), 0);
        const Matrix shape_derivatives =
            r_element.GetGeometry().ShapeFunctionsLocalGradients(r_element.GetIntegrationMethod())[0];
        Matrix identity = IdentityMatrix(3, 3);
        Vector strain(3);
        strain[0] = 2.0e-6; strain[1] = 0.0; strain[2] = 0.0;
        Vector stress(3);
        Matrix tangent(3, 3);
        Vector strain_work = strain;
        ConstitutiveLaw::Parameters values(
            r_element.GetGeometry(), r_element.GetProperties(), r_mp.GetProcessInfo());
        values.SetShapeFunctionsValues(shape_values);
        values.SetShapeFunctionsDerivatives(shape_derivatives);
        values.SetDeformationGradientF(identity);
        values.SetDeterminantF(1.0);
        values.SetStrainVector(strain_work);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(tangent);
        Flags& options = values.GetOptions();
        options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_mp.GetProcessInfo()[IS_CONVERGED] = true;
        r_element.GetConstitutiveLaw(0).CalculateMaterialResponseCauchy(values);
        r_element.GetConstitutiveLaw(0).FinalizeMaterialResponseCauchy(values);
        double damage_val = 0.0;
        r_element.GetConstitutiveLaw(0).GetValue(DAMAGE_VARIABLE, damage_val);
        KRATOS_EXPECT_TRUE(damage_val > 0.0);
    }
}

} // namespace Testing
} // namespace Kratos
