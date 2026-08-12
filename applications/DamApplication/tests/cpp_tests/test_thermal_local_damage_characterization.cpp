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

// Phase 4A: characterization of the thermal LOCAL damage constitutive-law
// family with the standard StructuralMechanicsApplication::SmallDisplacement
// element.
//
// This phase is characterization only: no production constitutive law or
// element is modified. The diagnostics below use a test-only subclass of the
// Simo-Ju law to read the internal damage state (non-invasive diagnostic
// access).
//
// The legacy element (SmallDisplacementThermoMechanicElement) drives the law
// through CalculateMaterialResponseCauchy / FinalizeMaterialResponseCauchy,
// while StructuralMechanicsApplication::SmallDisplacement drives the law
// through CalculateMaterialResponsePK2 / FinalizeMaterialResponsePK2
// (GetStressMeasure() returns StressMeasure_PK2). The thermal damage laws only
// override the CAUCHY entry points, so these tests establish what the PK2 path
// actually reaches.

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

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_stress_2D_law.hpp"

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

/// Test-only diagnostic subclass exposing the internal damage state
/// (non-invasive diagnostic access; no production code is modified).
class DiagnosticSimoJuLocalDamage3DLaw : public ThermalSimoJuLocalDamage3DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DiagnosticSimoJuLocalDamage3DLaw);

    DiagnosticSimoJuLocalDamage3DLaw() : ThermalSimoJuLocalDamage3DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new DiagnosticSimoJuLocalDamage3DLaw(*this));
    }

    /// Damage threshold / history variable (maximum historical equivalent
    /// strain).
    double GetThresholdVariable() const
    {
        return mpFlowRule->GetInternalVariables().EquivalentPlasticStrain;
    }

    /// Damage variable (d, the printed DeltaPlasticStrain).
    double GetDamageVariable() const
    {
        return mpFlowRule->GetInternalVariables().DeltaPlasticStrain;
    }
};

class DiagnosticSimoJuLocalDamagePlaneStrain2DLaw : public ThermalSimoJuLocalDamagePlaneStrain2DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DiagnosticSimoJuLocalDamagePlaneStrain2DLaw);

    DiagnosticSimoJuLocalDamagePlaneStrain2DLaw() : ThermalSimoJuLocalDamagePlaneStrain2DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new DiagnosticSimoJuLocalDamagePlaneStrain2DLaw(*this));
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

class DiagnosticSimoJuLocalDamagePlaneStress2DLaw : public ThermalSimoJuLocalDamagePlaneStress2DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DiagnosticSimoJuLocalDamagePlaneStress2DLaw);

    DiagnosticSimoJuLocalDamagePlaneStress2DLaw() : ThermalSimoJuLocalDamagePlaneStress2DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return ConstitutiveLaw::Pointer(new DiagnosticSimoJuLocalDamagePlaneStress2DLaw(*this));
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

/// Creates a model part with a single element of the given type and the
/// Simo-Ju local-damage law.
ModelPart& CreateDamageModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::size_t rDimension)
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

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    (*p_prop)[DAMAGE_THRESHOLD] = test_damage_threshold;
    (*p_prop)[STRENGTH_RATIO] = test_strength_ratio;
    (*p_prop)[FRACTURE_ENERGY] = test_fracture_energy;
    p_prop->SetValue(CONSTITUTIVE_LAW, DiagnosticSimoJuLocalDamage3DLaw().Clone());

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        element_nodes[i] = i + 1;
    }
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_prop);

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

/// Holds a ConstitutiveLaw::Parameters together with the storage it references
/// (shape-function values and deformation gradient), so that the references
/// remain valid for the duration of the law call.
struct LawContext
{
    Vector shape_values;
    Matrix deformation_gradient;
    ConstitutiveLaw::Parameters values;

    LawContext(
        ModelPart& rModelPart,
        Vector& rStrainVector,
        Vector& rStressVector,
        Matrix& rConstitutiveMatrix)
        : values(rModelPart.pGetElement(1)->GetGeometry(),
                 rModelPart.pGetElement(1)->GetProperties(),
                 rModelPart.GetProcessInfo())
    {
        Element& r_element = *rModelPart.pGetElement(1);
        const auto& r_geometry = r_element.GetGeometry();
        const auto& r_integration_method = r_element.GetIntegrationMethod();

        const auto& r_shape_functions = r_geometry.ShapeFunctionsValues(r_integration_method);
        const auto& r_shape_gradients = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);
        shape_values = row(r_shape_functions, 0);

        values.SetShapeFunctionsValues(shape_values);
        values.SetShapeFunctionsDerivatives(r_shape_gradients[0]);

        deformation_gradient = IdentityMatrix(3, 3);
        values.SetDeformationGradientF(deformation_gradient);
        values.SetDeterminantF(1.0);

        values.SetStrainVector(rStrainVector);
        values.SetStressVector(rStressVector);
        values.SetConstitutiveMatrix(rConstitutiveMatrix);
    }
};

/// Runs a Cauchy trial response (what the legacy element does during the
/// iterations) and returns the damaged stress, tangent and internal state.
template<class TDiagnosticLaw>
void RunCauchyTrial(
    ConstitutiveLaw& rLaw,
    ModelPart& rModelPart,
    const Vector& rTotalStrain,
    TDiagnosticLaw& rDiagnostic,
    Vector& rStress,
    Matrix& rTangent,
    double& rThreshold,
    double& rDamage)
{
    Vector strain_vector = rTotalStrain;
    rStress.resize(strain_vector.size(), false);
    rTangent.resize(strain_vector.size(), strain_vector.size(), false);
    LawContext context(rModelPart, strain_vector, rStress, rTangent);

    Flags& options = context.values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    rLaw.CalculateMaterialResponseCauchy(context.values);

    rThreshold = rDiagnostic.GetThresholdVariable();
    rDamage = rDiagnostic.GetDamageVariable();
}

/// Runs a Cauchy finalize (what the legacy element does at the end of the
/// step). IS_CONVERGED in the ProcessInfo selects commit / restore.
template<class TDiagnosticLaw>
void RunCauchyFinalize(
    ConstitutiveLaw& rLaw,
    ModelPart& rModelPart,
    const Vector& rTotalStrain,
    const bool rConverged,
    TDiagnosticLaw& rDiagnostic,
    double& rThreshold,
    double& rDamage)
{
    rModelPart.GetProcessInfo()[IS_CONVERGED] = rConverged;

    Vector strain_vector = rTotalStrain;
    Vector stress_vector;
    Matrix tangent;
    stress_vector.resize(strain_vector.size(), false);
    tangent.resize(strain_vector.size(), strain_vector.size(), false);
    LawContext context(rModelPart, strain_vector, stress_vector, tangent);

    Flags& options = context.values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    rLaw.FinalizeMaterialResponseCauchy(context.values);

    rThreshold = rDiagnostic.GetThresholdVariable();
    rDamage = rDiagnostic.GetDamageVariable();
}

/// Runs a PK2 trial response (what StructuralMechanics::SmallDisplacement does
/// during the iterations) and returns the damaged stress, tangent and state.
template<class TDiagnosticLaw>
void RunPK2Trial(
    ConstitutiveLaw& rLaw,
    ModelPart& rModelPart,
    const Vector& rTotalStrain,
    TDiagnosticLaw& rDiagnostic,
    Vector& rStress,
    Matrix& rTangent,
    double& rThreshold,
    double& rDamage)
{
    Vector strain_vector = rTotalStrain;
    rStress.resize(strain_vector.size(), false);
    rTangent.resize(strain_vector.size(), strain_vector.size(), false);
    LawContext context(rModelPart, strain_vector, rStress, rTangent);

    Flags& options = context.values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    rLaw.CalculateMaterialResponsePK2(context.values);

    rThreshold = rDiagnostic.GetThresholdVariable();
    rDamage = rDiagnostic.GetDamageVariable();
}

/// Runs a PK2 finalize (what StructuralMechanics::SmallDisplacement does at the
/// end of the step).
template<class TDiagnosticLaw>
void RunPK2Finalize(
    ConstitutiveLaw& rLaw,
    ModelPart& rModelPart,
    const Vector& rTotalStrain,
    const bool rConverged,
    TDiagnosticLaw& rDiagnostic,
    double& rThreshold,
    double& rDamage)
{
    rModelPart.GetProcessInfo()[IS_CONVERGED] = rConverged;

    Vector strain_vector = rTotalStrain;
    Vector stress_vector;
    Matrix tangent;
    stress_vector.resize(strain_vector.size(), false);
    tangent.resize(strain_vector.size(), strain_vector.size(), false);
    LawContext context(rModelPart, strain_vector, stress_vector, tangent);

    Flags& options = context.values.GetOptions();
    options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    rLaw.FinalizeMaterialResponsePK2(context.values);

    rThreshold = rDiagnostic.GetThresholdVariable();
    rDamage = rDiagnostic.GetDamageVariable();
}

/// Uniaxial strain vector used for the characterization states (3D).
Vector UniaxialStrain3D(const double rEpsilonX)
{
    Vector strain(6);
    strain[0] = rEpsilonX;
    strain[1] = 0.0;
    strain[2] = 0.0;
    strain[3] = 0.0;
    strain[4] = 0.0;
    strain[5] = 0.0;
    return strain;
}

/// Uniaxial strain vector used for the characterization states (2D).
Vector UniaxialStrain2D(const double rEpsilonX)
{
    Vector strain(3);
    strain[0] = rEpsilonX;
    strain[1] = 0.0;
    strain[2] = 0.0;
    return strain;
}

/// Applies a uniaxial displacement field (ux = eps*x, uy = -nu*eps*y,
/// uz = -nu*eps*z) to every node of the model part.
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

/// Sets a uniform temperature change on all nodes.
void ApplyTemperatureChange(ModelPart& rModelPart, const double rDeltaTemperature)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + rDeltaTemperature;
    }
}

/// Shared checks for the 2D inherited family: the common response must be
/// reached through both the Cauchy and PK2 entry points and behave identically.
template<class TDiagnosticLaw>
void Run2DInheritedFamilyChecks(
    ConstitutiveLaw& rLaw,
    TDiagnosticLaw& rDiag,
    ModelPart& rModelPart,
    const std::string& rLabel)
{
    // Elastic state: Cauchy and PK2 entry points must agree.
    Vector strain(3);
    strain[0] = 2.0e-6;
    strain[1] = 0.0;
    strain[2] = 0.0;
    Vector stress_c(3), stress_p(3);
    Matrix tangent_c(3, 3), tangent_p(3, 3);
    double tc, dc, tp, dp;
    RunCauchyTrial(rLaw, rModelPart, strain, rDiag, stress_c, tangent_c, tc, dc);
    RunPK2Trial(rLaw, rModelPart, strain, rDiag, stress_p, tangent_p, tp, dp);
    KRATOS_EXPECT_NEAR(tc, test_damage_threshold, 1.0e-15);
    KRATOS_EXPECT_NEAR(tp, test_damage_threshold, 1.0e-15);
    KRATOS_EXPECT_NEAR(dp, dc, 1.0e-15);
    KRATOS_EXPECT_NEAR(stress_p[0], stress_c[0], 1.0e-9);

    // Damage initiation (converged finalize through the Cauchy entry).
    strain[0] = 2.0e-5;
    RunCauchyTrial(rLaw, rModelPart, strain, rDiag, stress_c, tangent_c, tc, dc);
    RunCauchyFinalize(rLaw, rModelPart, strain, true, rDiag, tc, dc);
    const double committed_threshold = tc;
    KRATOS_EXPECT_TRUE(dc > 0.0);

    // PK2 finalize commits the same state.
    RunPK2Finalize(rLaw, rModelPart, strain, true, rDiag, tp, dp);
    KRATOS_EXPECT_NEAR(tp, committed_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(dp, dc, 1.0e-12);

    // Non-converged PK2 finalize with a larger load: restore.
    strain[0] = 3.0e-5;
    RunPK2Finalize(rLaw, rModelPart, strain, false, rDiag, tp, dp);
    KRATOS_EXPECT_NEAR(tp, committed_threshold, 1.0e-12);

    // Unload / reload: irreversible.
    strain[0] = 4.0e-6;
    RunCauchyTrial(rLaw, rModelPart, strain, rDiag, stress_c, tangent_c, tc, dc);
    KRATOS_EXPECT_NEAR(tc, committed_threshold, 1.0e-12);
    strain[0] = 2.8e-5;
    RunCauchyTrial(rLaw, rModelPart, strain, rDiag, stress_c, tangent_c, tc, dc);
    RunCauchyFinalize(rLaw, rModelPart, strain, true, rDiag, tc, dc);
    KRATOS_EXPECT_TRUE(tc > committed_threshold);

    // Thermal coupling: free thermal total strain must not damage in the
    // PK2 path (inherited common response).
    ApplyTemperatureChange(rModelPart, 50.0);
    Vector thermal_total(3);
    thermal_total[0] = test_thermal_expansion * 50.0;
    thermal_total[1] = test_thermal_expansion * 50.0;
    thermal_total[2] = 0.0;
    RunPK2Trial(rLaw, rModelPart, thermal_total, rDiag, stress_p, tangent_p, tp, dp);
    RunPK2Finalize(rLaw, rModelPart, thermal_total, true, rDiag, tp, dp);
    KRATOS_EXPECT_NEAR(dp, dc, 1.0e-12);

    std::cout << "[4B] 2D " << rLabel << ": d_damaged=" << dc << std::endl;
}

} // namespace

//************************************************************************************
// 1. Constitutive-response characterization: elastic -> initiation -> damaged
//    -> unload -> reload (legacy Cauchy path, 3D Simo-Ju).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_ResponseCharacterization3D, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "Char3D", "SmallDisplacementSolidElement3D8N", 3);

    ConstitutiveLaw& r_law = *r_model_part.pGetElement(1)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(r_law);
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    // Uniaxial-strain elastic modulus (lambda + 2 mu) for epsilon_y = epsilon_z = 0.
    const double lambda_2mu =
        test_young_modulus * (1.0 - test_poisson_ratio)
        / ((1.0 + test_poisson_ratio) * (1.0 - 2.0 * test_poisson_ratio));

    // State A: comfortably below damage initiation (elastic trial, d = 0,
    // threshold = material threshold, stress = (lambda+2mu)*epsilon).
    double stress_a, threshold_a, damage_a;
    Matrix tangent_a;
    Vector stress_vector_a;
    RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(2.0e-6), r_diag, stress_vector_a, tangent_a, threshold_a, damage_a);
    KRATOS_EXPECT_NEAR(threshold_a, test_damage_threshold, 1.0e-15);
    KRATOS_EXPECT_NEAR(damage_a, 0.0, 1.0e-15);
    KRATOS_EXPECT_NEAR(stress_vector_a[0], lambda_2mu * 2.0e-6, 1.0e-6);

    // State B: just beyond damage initiation. The TRIAL still uses the
    // committed history (d = 0); the damage is committed at the finalize.
    double stress_b, threshold_b, damage_b;
    Matrix tangent_b;
    Vector stress_vector_b;
    RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(1.2e-5), r_diag, stress_vector_b, tangent_b, threshold_b, damage_b);
    KRATOS_EXPECT_NEAR(threshold_b, test_damage_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(damage_b, 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(stress_vector_b[0], lambda_2mu * 1.2e-5, 1.0e-4);

    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(1.2e-5), true, r_diag, threshold_b, damage_b);
    KRATOS_EXPECT_TRUE(threshold_b > test_damage_threshold);
    KRATOS_EXPECT_TRUE(damage_b > 0.0);

    // State C: clearly damaged state (monotonic loading, committed at each
    // converged finalize).
    for (double eps : {1.6e-5, 2.0e-5, 2.4e-5}) {
        double stress_c, threshold_c, damage_c;
        Matrix tangent_c;
        Vector stress_vector_c;
        RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(eps), r_diag, stress_vector_c, tangent_c, threshold_c, damage_c);
        RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(eps), true, r_diag, threshold_c, damage_c);
        KRATOS_EXPECT_TRUE(damage_c > 0.0);
    }
    double stress_c, threshold_c, damage_c;
    Matrix tangent_c;
    Vector stress_vector_c;
    RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(2.4e-5), r_diag, stress_vector_c, tangent_c, threshold_c, damage_c);
    const double last_threshold = threshold_c;
    const double last_damage = damage_c;
    KRATOS_EXPECT_TRUE(damage_c > 0.4);

    // State D: unload substantially (below the history). Damage is
    // irreversible: threshold and damage stay.
    double stress_d, threshold_d, damage_d;
    Matrix tangent_d;
    Vector stress_vector_d;
    RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(4.0e-6), r_diag, stress_vector_d, tangent_d, threshold_d, damage_d);
    KRATOS_EXPECT_NEAR(threshold_d, last_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(damage_d, last_damage, 1.0e-12);

    // State E: reload below the previous maximum -> no damage growth.
    double stress_e, threshold_e, damage_e;
    Matrix tangent_e;
    Vector stress_vector_e;
    RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(2.0e-5), r_diag, stress_vector_e, tangent_e, threshold_e, damage_e);
    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(2.0e-5), true, r_diag, threshold_e, damage_e);
    KRATOS_EXPECT_NEAR(threshold_e, last_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(damage_e, last_damage, 1.0e-12);

    // State F: reload beyond the previous maximum -> damage evolves.
    double stress_f, threshold_f, damage_f;
    Matrix tangent_f;
    Vector stress_vector_f;
    RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(2.8e-5), r_diag, stress_vector_f, tangent_f, threshold_f, damage_f);
    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(2.8e-5), true, r_diag, threshold_f, damage_f);
    KRATOS_EXPECT_TRUE(threshold_f > last_threshold);
    KRATOS_EXPECT_TRUE(damage_f > last_damage);

    std::cout << "[4A] 3D response: threshold A=" << threshold_a
              << " dA=" << damage_a << " dB=" << damage_b
              << " dC=" << damage_c << " dD=" << damage_d
              << " dE=" << damage_e << " dF=" << damage_f << std::endl;
}

//************************************************************************************
// 2. Dispatch: Cauchy (thermal, legacy) vs PK2 (candidate) under thermal load.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_DispatchCauchyVsPK2Thermal, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "Dispatch3D", "SmallDisplacementSolidElement3D8N", 3);
    ConstitutiveLaw& r_law = *r_model_part.pGetElement(1)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(r_law);
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    ApplyUniaxialState(r_model_part, 0.0);
    ApplyTemperatureChange(r_model_part, 50.0);

    const double thermal_strain = test_thermal_expansion * 50.0;  // alpha * dT = 5e-4

    // Case 1: restrained uniform thermal loading. The element supplies a ZERO
    // total strain (no displacement).
    Vector zero_strain = ZeroVector(6);

    // Cauchy path (legacy): the law subtracts the thermal strain, so the
    // mechanical strain is -thermal (compression). Restrained heating is
    // correctly captured as compression.
    double threshold_cauchy, damage_cauchy;
    Vector stress_cauchy;
    Matrix tangent_cauchy;
    RunCauchyTrial(r_law, r_model_part, zero_strain, r_diag, stress_cauchy, tangent_cauchy, threshold_cauchy, damage_cauchy);
    // Pure isotropic thermal compression: sigma = -3K*alpha*dT on the diagonal
    // (K = bulk modulus; 3K = E/(1-2nu)).
    const double three_bulk_modulus = test_young_modulus / (1.0 - 2.0 * test_poisson_ratio);
    KRATOS_EXPECT_NEAR(stress_cauchy[0], -three_bulk_modulus * thermal_strain, 1.0e-6);
    KRATOS_EXPECT_NEAR(damage_cauchy, 0.0, 1.0e-15);

    // After Phase 4B the PK2 path executes the SAME thermal damage response:
    // the temperature is used, the thermal strain is subtracted and the
    // restrained-heating compression matches the Cauchy path exactly.
    double threshold_pk2, damage_pk2;
    Vector stress_pk2;
    Matrix tangent_pk2;
    RunPK2Trial(r_law, r_model_part, zero_strain, r_diag, stress_pk2, tangent_pk2, threshold_pk2, damage_pk2);
    KRATOS_EXPECT_NEAR(stress_pk2[0], stress_cauchy[0],
                       std::max(comparison_absolute_tolerance,
                                comparison_relative_tolerance * std::abs(stress_cauchy[0])));
    KRATOS_EXPECT_NEAR(damage_pk2, 0.0, 1.0e-15);

    // Case 2: free thermal expansion supplied as the TOTAL strain
    // (alpha*dT on the diagonal, as an SMA element would compute from the
    // displacement field of a free thermal expansion).
    Vector thermal_total_strain(6);
    thermal_total_strain[0] = thermal_strain;
    thermal_total_strain[1] = thermal_strain;
    thermal_total_strain[2] = thermal_strain;
    thermal_total_strain[3] = 0.0;
    thermal_total_strain[4] = 0.0;
    thermal_total_strain[5] = 0.0;

    // Cauchy path (legacy): the thermal strain is subtracted -> mechanical
    // strain zero -> no stress and no damage (correct: free thermal expansion
    // produces no stress).
    RunCauchyTrial(r_law, r_model_part, thermal_total_strain, r_diag, stress_cauchy, tangent_cauchy, threshold_cauchy, damage_cauchy);
    for (double val : stress_cauchy) {
        KRATOS_EXPECT_NEAR(val, 0.0, 1.0e-9);
    }
    KRATOS_EXPECT_NEAR(damage_cauchy, 0.0, 1.0e-15);

    // PK2 path (candidate): identical result - no artificial thermal stress or
    // damage (the thermal strain is subtracted in both paths).
    RunPK2Trial(r_law, r_model_part, thermal_total_strain, r_diag, stress_pk2, tangent_pk2, threshold_pk2, damage_pk2);
    for (double val : stress_pk2) {
        KRATOS_EXPECT_NEAR(val, 0.0, 1.0e-9);
    }
    KRATOS_EXPECT_NEAR(damage_pk2, 0.0, 1.0e-15);

    // After Phase 4B, the PK2 path executes the SAME thermal damage response as
    // the Cauchy path: a converged FINALIZE must NOT commit artificial damage
    // for a free thermal expansion (mechanical strain zero in both paths).
    double threshold_f_c, damage_f_c;
    RunCauchyFinalize(r_law, r_model_part, thermal_total_strain, true, r_diag, threshold_f_c, damage_f_c);
    KRATOS_EXPECT_NEAR(damage_f_c, 0.0, 1.0e-15);
    KRATOS_EXPECT_NEAR(threshold_f_c, test_damage_threshold, 1.0e-15);

    double threshold_f_p, damage_f_p;
    RunPK2Finalize(r_law, r_model_part, thermal_total_strain, true, r_diag, threshold_f_p, damage_f_p);
    std::cout << "[4B] dispatch: Cauchy damage after finalize = " << damage_f_c
              << ", PK2 damage after finalize = " << damage_f_p << std::endl;
    KRATOS_EXPECT_NEAR(damage_f_p, 0.0, 1.0e-15);
    KRATOS_EXPECT_NEAR(threshold_f_p, test_damage_threshold, 1.0e-15);
    // PK2 and Cauchy produce the same stress for the same input.
    KRATOS_EXPECT_NEAR(stress_pk2[0], stress_cauchy[0],
                       std::max(comparison_absolute_tolerance,
                                comparison_relative_tolerance * std::abs(stress_cauchy[0])));
}

//************************************************************************************
// 3. Lifecycle: Requires* hooks and the SMA element InitializeSolutionStep.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_LifecycleRequiresHooks, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_legacy = CreateDamageModelPart(
        model, "LifeLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3);

    // The damage law is stateful through the FINALIZATION, but no local-damage
    // state is initialized through the InitializeMaterialResponse hooks: the
    // history/threshold is initialized in InitializeMaterial and managed in the
    // finalization. So the initialization hooks are not required.
    ConstitutiveLaw& r_law = *r_legacy.pGetElement(1)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    KRATOS_EXPECT_FALSE(r_law.RequiresInitializeMaterialResponse());
    KRATOS_EXPECT_TRUE(r_law.RequiresFinalizeMaterialResponse());

    // The SMA element calls InitializeMaterialResponse(Values, PK2) from its
    // InitializeSolutionStep only when RequiresInitializeMaterialResponse()
    // returns true. With the override to false, the SMA lifecycle now runs
    // without exception and does not change the committed damage state.
    ModelPart& r_candidate = CreateDamageModelPart(
        model, "LifeCandidate", "SmallDisplacementElement3D8N", 3);
    auto& r_candidate_element = *r_candidate.pGetElement(1);
    r_candidate_element.Initialize(r_candidate.GetProcessInfo());

    bool threw = false;
    try {
        r_candidate_element.InitializeSolutionStep(r_candidate.GetProcessInfo());
    } catch (...) {
        threw = true;
    }
    std::cout << "[4B] lifecycle: SMA element InitializeSolutionStep throws with "
                 "the damage law: " << (threw ? "yes" : "no") << std::endl;
    KRATOS_EXPECT_FALSE(threw);
}
//************************************************************************************
// 4. State commit/restore: Cauchy finalize honours IS_CONVERGED.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_ISConvergedCommitRestoreCauchy, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "Conv3D", "SmallDisplacementSolidElement3D8N", 3);
    ConstitutiveLaw& r_law = *r_model_part.pGetElement(1)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(r_law);
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    // Converged damaged step: commit the trial history.
    double threshold, damage;
    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(2.0e-5), true, r_diag, threshold, damage);
    const double committed_threshold = threshold;
    KRATOS_EXPECT_TRUE(committed_threshold > test_damage_threshold);

    // Larger trial load, NOT converged: the history must be restored (kept at
    // the previous equilibrium value).
    double threshold_restored, damage_restored;
    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(3.0e-5), false, r_diag, threshold_restored, damage_restored);
    KRATOS_EXPECT_NEAR(threshold_restored, committed_threshold, 1.0e-12);
    KRATOS_EXPECT_NEAR(damage_restored, damage, 1.0e-12);

    // Repeat the larger trial load with convergence: the new state is
    // committed.
    double threshold_recommit, damage_recommit;
    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(3.0e-5), true, r_diag, threshold_recommit, damage_recommit);
    KRATOS_EXPECT_TRUE(threshold_recommit > committed_threshold);
    KRATOS_EXPECT_TRUE(damage_recommit > damage);
}

//************************************************************************************
// 5. PK2 finalize ignores IS_CONVERGED (candidate lifecycle divergence).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_PK2FinalizeHonoursISConverged, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "PK2Conv3D", "SmallDisplacementSolidElement3D8N", 3);
    ConstitutiveLaw& r_law = *r_model_part.pGetElement(1)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(r_law);
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    // Converged damaged step via PK2 finalize.
    double threshold, damage;
    RunPK2Finalize(r_law, r_model_part, UniaxialStrain3D(2.0e-5), true, r_diag, threshold, damage);
    const double committed_threshold = threshold;

    // Non-converged larger trial via PK2 finalize: the PK2 path now honours
    // IS_CONVERGED and RESTORES the previous equilibrium state (it no longer
    // commits unconditionally as the old inherited HyperElastic3DLaw finalize
    // did).
    double threshold_after, damage_after;
    RunPK2Finalize(r_law, r_model_part, UniaxialStrain3D(3.0e-5), false, r_diag, threshold_after, damage_after);
    std::cout << "[4B] lifecycle PK2: committed=" << committed_threshold
              << " after non-converged finalize=" << threshold_after << std::endl;
    KRATOS_EXPECT_NEAR(threshold_after, committed_threshold, 1.0e-12);

    // Repeated with convergence: the new state is committed.
    double threshold_recommit, damage_recommit;
    RunPK2Finalize(r_law, r_model_part, UniaxialStrain3D(3.0e-5), true, r_diag, threshold_recommit, damage_recommit);
    KRATOS_EXPECT_TRUE(threshold_recommit > committed_threshold);
    KRATOS_EXPECT_TRUE(damage_recommit > damage);
}
//************************************************************************************
// 6. Legacy element full lifecycle with the damage law.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_ElementLegacyLifecycle, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "ElemLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3);
    auto& r_element = *r_model_part.pGetElement(1);
    r_element.Initialize(r_model_part.GetProcessInfo());
    r_element.InitializeSolutionStep(r_model_part.GetProcessInfo());
    r_element.InitializeNonLinearIteration(r_model_part.GetProcessInfo());

    // Apply a loading state beyond damage initiation.
    ApplyUniaxialState(r_model_part, 2.0e-5);

    Matrix lhs;
    Vector rhs;
    r_element.CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());
    KRATOS_EXPECT_EQ(lhs.size1(), rhs.size());
    KRATOS_EXPECT_TRUE(rhs.size() > 0);

    // No damage committed yet (only a trial).
    std::vector<double> thresholds;
    r_element.CalculateOnIntegrationPoints(STATE_VARIABLE, thresholds, r_model_part.GetProcessInfo());
    KRATOS_EXPECT_EQ(thresholds.size(), 8);  // H8 GI_GAUSS_2
    KRATOS_EXPECT_NEAR(thresholds[0], test_damage_threshold, 1.0e-12);

    // Converged finalize: commit the damage state.
    r_model_part.GetProcessInfo()[IS_CONVERGED] = true;
    r_element.FinalizeSolutionStep(r_model_part.GetProcessInfo());
    r_element.CalculateOnIntegrationPoints(STATE_VARIABLE, thresholds, r_model_part.GetProcessInfo());
    const double committed_threshold = thresholds[0];
    KRATOS_EXPECT_TRUE(committed_threshold > test_damage_threshold);

    // Non-converged finalize with a larger load: restore the equilibrium state
    // (the trial history of the rejected step must NOT be committed).
    ApplyUniaxialState(r_model_part, 3.0e-5);
    r_model_part.GetProcessInfo()[IS_CONVERGED] = false;
    r_element.FinalizeSolutionStep(r_model_part.GetProcessInfo());
    r_element.CalculateOnIntegrationPoints(STATE_VARIABLE, thresholds, r_model_part.GetProcessInfo());
    KRATOS_EXPECT_NEAR(thresholds[0], committed_threshold, 1.0e-12);  // restored to the previous equilibrium state
}

//************************************************************************************
// 7. Repeated-response side effects (law level): read/output must not advance
//    the damage history.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_RepeatedResponseNoAdvance, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "Repeat3D", "SmallDisplacementSolidElement3D8N", 3);
    ConstitutiveLaw& r_law = *r_model_part.pGetElement(1)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamage3DLaw&>(r_law);
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    // Damage the state first.
    double threshold, damage;
    RunCauchyFinalize(r_law, r_model_part, UniaxialStrain3D(2.0e-5), true, r_diag, threshold, damage);
    const double committed_threshold = threshold;
    const double committed_damage = damage;

    // Repeated Cauchy trials at a fixed state: the history must NOT advance.
    for (std::size_t i = 0; i < 5; ++i) {
        Vector stress;
        Matrix tangent;
        double t, d;
        RunCauchyTrial(r_law, r_model_part, UniaxialStrain3D(2.0e-5), r_diag, stress, tangent, t, d);
        KRATOS_EXPECT_NEAR(t, committed_threshold, 1.0e-12);
        KRATOS_EXPECT_NEAR(d, committed_damage, 1.0e-12);
    }

    // Repeated PK2 trials: same behaviour (the PK2 trial does not commit).
    for (std::size_t i = 0; i < 5; ++i) {
        Vector stress;
        Matrix tangent;
        double t, d;
        RunPK2Trial(r_law, r_model_part, UniaxialStrain3D(2.0e-5), r_diag, stress, tangent, t, d);
        KRATOS_EXPECT_NEAR(t, committed_threshold, 1.0e-12);
        KRATOS_EXPECT_NEAR(d, committed_damage, 1.0e-12);
    }
}

//************************************************************************************
// 8. Element-level legacy vs candidate core response (pre-damage and damaged).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_ElementLegacyVsCandidateCore, KratosDamFastSuite)
{
    Model model;

    // Legacy element: full lifecycle works with the damage law.
    ModelPart& r_legacy = CreateDamageModelPart(
        model, "CoreLegacy", "SmallDisplacementThermoMechanicElement3D8N", 3);
    auto& r_legacy_element = *r_legacy.pGetElement(1);
    r_legacy_element.Initialize(r_legacy.GetProcessInfo());
    r_legacy_element.InitializeSolutionStep(r_legacy.GetProcessInfo());
    r_legacy_element.InitializeNonLinearIteration(r_legacy.GetProcessInfo());

    // Candidate element: the real SMA lifecycle (InitializeSolutionStep now
    // runs thanks to RequiresInitializeMaterialResponse() == false).
    ModelPart& r_candidate = CreateDamageModelPart(
        model, "CoreCandidate", "SmallDisplacementElement3D8N", 3);
    auto& r_candidate_element = *r_candidate.pGetElement(1);
    r_candidate_element.Initialize(r_candidate.GetProcessInfo());
    r_candidate_element.InitializeSolutionStep(r_candidate.GetProcessInfo());

    // Pre-damage elastic state.
    ApplyUniaxialState(r_legacy, 2.0e-6);
    ApplyUniaxialState(r_candidate, 2.0e-6);
    Matrix lhs_legacy, lhs_candidate;
    Vector rhs_legacy, rhs_candidate;
    r_legacy_element.CalculateLocalSystem(lhs_legacy, rhs_legacy, r_legacy.GetProcessInfo());
    r_candidate_element.CalculateLocalSystem(lhs_candidate, rhs_candidate, r_candidate.GetProcessInfo());
    KRATOS_EXPECT_EQ(lhs_legacy.size1(), lhs_candidate.size1());
    KRATOS_EXPECT_EQ(rhs_legacy.size(), rhs_candidate.size());
    for (std::size_t i = 0; i < rhs_legacy.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs_candidate[i], rhs_legacy[i],
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(rhs_legacy[i])));
    }
    for (std::size_t i = 0; i < lhs_legacy.size1(); ++i) {
        for (std::size_t j = 0; j < lhs_legacy.size2(); ++j) {
            KRATOS_EXPECT_NEAR(lhs_candidate(i, j), lhs_legacy(i, j),
                               std::max(comparison_absolute_tolerance,
                                        comparison_relative_tolerance * std::abs(lhs_legacy(i, j))));
        }
    }

    // Damaged state: trial does not commit, both elements produce the same
    // trial LHS/RHS.
    ApplyUniaxialState(r_legacy, 2.0e-5);
    ApplyUniaxialState(r_candidate, 2.0e-5);
    r_legacy_element.CalculateLocalSystem(lhs_legacy, rhs_legacy, r_legacy.GetProcessInfo());
    r_candidate_element.CalculateLocalSystem(lhs_candidate, rhs_candidate, r_candidate.GetProcessInfo());
    for (std::size_t i = 0; i < rhs_legacy.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs_candidate[i], rhs_legacy[i],
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(rhs_legacy[i])));
    }

    // Commit the damaged state through both elements' FinalizeSolutionStep.
    r_legacy.GetProcessInfo()[IS_CONVERGED] = true;
    r_candidate.GetProcessInfo()[IS_CONVERGED] = true;
    r_legacy_element.FinalizeSolutionStep(r_legacy.GetProcessInfo());
    r_candidate_element.FinalizeSolutionStep(r_candidate.GetProcessInfo());

    std::vector<double> legacy_thresholds, candidate_thresholds;
    r_legacy_element.CalculateOnIntegrationPoints(STATE_VARIABLE, legacy_thresholds, r_legacy.GetProcessInfo());
    r_candidate_element.CalculateOnIntegrationPoints(STATE_VARIABLE, candidate_thresholds, r_candidate.GetProcessInfo());
    KRATOS_EXPECT_EQ(legacy_thresholds.size(), 8);
    KRATOS_EXPECT_EQ(candidate_thresholds.size(), 8);
    for (std::size_t i = 0; i < legacy_thresholds.size(); ++i) {
        KRATOS_EXPECT_NEAR(candidate_thresholds[i], legacy_thresholds[i],
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(legacy_thresholds[i])));
        KRATOS_EXPECT_TRUE(legacy_thresholds[i] > test_damage_threshold);
    }
}

//************************************************************************************
// 9. 2D family: plane strain (Q4).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_2DPlaneStrainFamily, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "PS2D", "SmallDisplacementSolidElement2D4N", 2);

    // The property holds the 3D diagnostic law; replace with the 2D one.
    auto p_2d_law = Kratos::make_shared<DiagnosticSimoJuLocalDamagePlaneStrain2DLaw>();
    r_model_part.pGetElement(1)->GetProperties().SetValue(CONSTITUTIVE_LAW, p_2d_law);
    ConstitutiveLaw& r_law = *p_2d_law;
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    auto& r_2d_diagnostic = dynamic_cast<DiagnosticSimoJuLocalDamagePlaneStrain2DLaw&>(r_law);

    // Elastic.
    Vector strain(3);
    strain[0] = 2.0e-6;
    strain[1] = 0.0;
    strain[2] = 0.0;
    Vector stress(3);
    Matrix tangent(3, 3);
    double threshold, damage;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    KRATOS_EXPECT_NEAR(threshold, test_damage_threshold, 1.0e-15);
    KRATOS_EXPECT_NEAR(damage, 0.0, 1.0e-15);

    // Damaged state.
    strain[0] = 2.0e-5;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    RunCauchyFinalize(r_law, r_model_part, strain, true, r_2d_diagnostic, threshold, damage);
    const double committed_threshold = threshold;
    KRATOS_EXPECT_TRUE(damage > 0.0);

    // Unload.
    strain[0] = 4.0e-6;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    KRATOS_EXPECT_NEAR(threshold, committed_threshold, 1.0e-12);

    // Reload beyond.
    strain[0] = 2.8e-5;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    RunCauchyFinalize(r_law, r_model_part, strain, true, r_2d_diagnostic, threshold, damage);
    KRATOS_EXPECT_TRUE(threshold > committed_threshold);

    std::cout << "[4A] 2D plane strain: d_elastic=0, d_damaged=" << damage << std::endl;
}

//************************************************************************************
// 10. 2D family: plane stress (Q4).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_2DPlaneStressFamily, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = CreateDamageModelPart(
        model, "PV2D", "SmallDisplacementSolidElement2D4N", 2);

    auto p_2d_law = Kratos::make_shared<DiagnosticSimoJuLocalDamagePlaneStress2DLaw>();
    r_model_part.pGetElement(1)->GetProperties().SetValue(CONSTITUTIVE_LAW, p_2d_law);
    ConstitutiveLaw& r_law = *p_2d_law;
    r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                             r_model_part.pGetElement(1)->GetGeometry(),
                             Vector());

    auto& r_2d_diagnostic = dynamic_cast<DiagnosticSimoJuLocalDamagePlaneStress2DLaw&>(r_law);

    Vector strain(3);
    strain[0] = 2.0e-6;
    strain[1] = 0.0;
    strain[2] = 0.0;
    Vector stress(3);
    Matrix tangent(3, 3);
    double threshold, damage;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    KRATOS_EXPECT_NEAR(threshold, test_damage_threshold, 1.0e-15);
    KRATOS_EXPECT_NEAR(damage, 0.0, 1.0e-15);

    strain[0] = 2.0e-5;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    RunCauchyFinalize(r_law, r_model_part, strain, true, r_2d_diagnostic, threshold, damage);
    const double committed_threshold = threshold;
    KRATOS_EXPECT_TRUE(damage > 0.0);

    strain[0] = 4.0e-6;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    KRATOS_EXPECT_NEAR(threshold, committed_threshold, 1.0e-12);

    strain[0] = 2.8e-5;
    RunCauchyTrial(r_law, r_model_part, strain, r_2d_diagnostic, stress, tangent, threshold, damage);
    RunCauchyFinalize(r_law, r_model_part, strain, true, r_2d_diagnostic, threshold, damage);
    KRATOS_EXPECT_TRUE(threshold > committed_threshold);

    std::cout << "[4A] 2D plane stress: d_damaged=" << damage << std::endl;
}

//************************************************************************************
// 11. 2D inherited family: the common response is inherited and PK2 == Cauchy.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalLocalDamage_2DInheritedFamily, KratosDamFastSuite)
{
    for (bool plane_stress : {false, true}) {
        Model model;
        const std::string law_name = plane_stress ? "PlaneStress" : "PlaneStrain";
        ModelPart& r_model_part = CreateDamageModelPart(
            model, "Inherit2D_" + law_name, "SmallDisplacementSolidElement2D4N", 2);

        if (plane_stress) {
            auto p_2d_law = Kratos::make_shared<DiagnosticSimoJuLocalDamagePlaneStress2DLaw>();
            r_model_part.pGetElement(1)->GetProperties().SetValue(CONSTITUTIVE_LAW, p_2d_law);
            ConstitutiveLaw& r_law = *p_2d_law;
            r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                                     r_model_part.pGetElement(1)->GetGeometry(),
                                     Vector());
            auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamagePlaneStress2DLaw&>(r_law);
            Run2DInheritedFamilyChecks(r_law, r_diag, r_model_part, law_name);
        } else {
            auto p_2d_law = Kratos::make_shared<DiagnosticSimoJuLocalDamagePlaneStrain2DLaw>();
            r_model_part.pGetElement(1)->GetProperties().SetValue(CONSTITUTIVE_LAW, p_2d_law);
            ConstitutiveLaw& r_law = *p_2d_law;
            r_law.InitializeMaterial(r_model_part.pGetElement(1)->GetProperties(),
                                     r_model_part.pGetElement(1)->GetGeometry(),
                                     Vector());
            auto& r_diag = dynamic_cast<DiagnosticSimoJuLocalDamagePlaneStrain2DLaw&>(r_law);
            Run2DInheritedFamilyChecks(r_law, r_diag, r_model_part, law_name);
        }
    }
}

} // namespace Testing
} // namespace Kratos
