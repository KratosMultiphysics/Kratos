// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Dam thermal linear-elastic contract. The non-nodal thermal laws (3D, plane
// strain, plane stress) on the StructuralMechanics small-displacement element
// reproduce the analytical restrained-expansion response and serialize
// correctly.
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
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"
#include "custom_utilities/nodal_young_modulus_utilities.h"

// StructuralMechanicsApplication mechanical reference laws (zero thermal strain)
#include "custom_constitutive/flexible_elastic_isotropic_3d.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances.
constexpr double tle_comparison_absolute_tolerance = 1.0e-12;
constexpr double tle_comparison_relative_tolerance = 1.0e-10;
constexpr double tle_machine_precision_allowance = 1.0e-15;

/// Material data shared by all model parts.
constexpr double tle_test_young_modulus = 2.0e7;
constexpr double tle_test_poisson_ratio = 0.2;
constexpr double tle_test_density = 2400.0;
constexpr double tle_test_thermal_expansion = 1.0e-5;
constexpr double tle_test_reference_temperature = 20.0;
constexpr double tle_test_thickness = 0.15;
constexpr double test_delta_temperature = 25.0;

/// Creates the constitutive law for the given configuration.
ConstitutiveLaw::Pointer tle_CreateThermalLaw(const std::string& rLawName)
{
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        return ThermalLinearElastic2DPlaneStrain().Clone();
    } else if (rLawName == "ThermalLinearElastic2DPlaneStress") {
        return ThermalLinearElastic2DPlaneStress().Clone();
    } else {
        return ThermalLinearElastic3DLaw().Clone();
    }
}

/// Element responses compared for every scenario.
struct FamilyResponse
{
    Matrix lhs;
    Vector rhs;
    Matrix independent_lhs;
    Vector independent_rhs;
    std::vector<Vector> strain_vectors;
    std::vector<Vector> cauchy_stress_vectors;
    std::vector<Vector> pk2_stress_vectors;
};

/// Creates one of the two numerically identical model parts for a 2D or 3D
/// configuration.
ModelPart& CreateFamilyModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rLawName,
    const std::size_t rDimension,
    const ModelPart::IndexType rNodeIdOffset)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);

    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    // Optional nodal Young's modulus field for the accessor configuration.
    r_model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS);
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    // Geometry from the registered prototype of the element.
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
        r_model_part.CreateNewNode(rNodeIdOffset + i + 1, x, y, z);
    }

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = tle_test_young_modulus;
    (*p_prop)[POISSON_RATIO] = tle_test_poisson_ratio;
    (*p_prop)[DENSITY] = tle_test_density;
    (*p_prop)[THERMAL_EXPANSION] = tle_test_thermal_expansion;
    if (rDimension == 2) {
        (*p_prop)[THICKNESS] = tle_test_thickness;
    }
    p_prop->SetValue(CONSTITUTIVE_LAW, tle_CreateThermalLaw(rLawName));

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        element_nodes[i] = rNodeIdOffset + i + 1;
    }
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_prop);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = tle_test_reference_temperature;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = tle_test_reference_temperature;
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Prescribes the state of a scenario:
///   scenario 0 (T0): reference temperature, small non-zero displacement field;
///   scenario 1 (T1): zero displacement, uniform temperature increment;
///   scenario 2 (T2): zero displacement, non-uniform nodal temperature field.
void PrescribeFamilyScenario(ModelPart& rModelPart, const std::size_t rScenarioIndex, const std::size_t rDimension)
{
    KRATOS_TRY;

    Matrix displacement_gradient = ZeroMatrix(3, 3);
    array_1d<double, 3> translation;
    if (rScenarioIndex == 0) {
        if (rDimension == 2) {
            displacement_gradient(0, 0) = 2.0e-3;
            displacement_gradient(0, 1) = 3.0e-4;
            displacement_gradient(1, 0) = 5.0e-4;
            displacement_gradient(1, 1) = -1.0e-3;
            translation[0] = 1.0e-2;
            translation[1] = -2.0e-2;
        } else {
            displacement_gradient(0, 0) = 2.0e-3;
            displacement_gradient(0, 1) = 3.0e-4;
            displacement_gradient(0, 2) = 1.0e-4;
            displacement_gradient(1, 0) = 3.0e-4;
            displacement_gradient(1, 1) = -1.0e-3;
            displacement_gradient(1, 2) = 2.0e-4;
            displacement_gradient(2, 0) = 1.0e-4;
            displacement_gradient(2, 1) = 2.0e-4;
            displacement_gradient(2, 2) = 5.0e-4;
            translation[0] = 1.0e-2;
            translation[1] = -2.0e-2;
            translation[2] = 5.0e-3;
        }
    }

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        if (rScenarioIndex == 0) {
            array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            noalias(r_displacement) = prod(displacement_gradient, r_node.GetInitialPosition());
            noalias(r_displacement) += translation;
            // Updated coordinates so that the reference configuration (current
            // position - total displacement) is the initial one.
            r_node.X() = r_node.X0() + r_displacement[0];
            r_node.Y() = r_node.Y0() + r_displacement[1];
            r_node.Z() = r_node.Z0() + r_displacement[2];
        } else {
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
            r_node.X() = r_node.X0();
            r_node.Y() = r_node.Y0();
            r_node.Z() = r_node.Z0();
        }

        double nodal_temperature = tle_test_reference_temperature;
        if (rScenarioIndex == 1) {
            nodal_temperature = tle_test_reference_temperature + test_delta_temperature;
        } else if (rScenarioIndex == 2) {
            nodal_temperature = tle_test_reference_temperature + 5.0 * static_cast<double>(node_index);
        }
        r_node.FastGetSolutionStepValue(TEMPERATURE) = nodal_temperature;

        array_1d<double, 3>& r_volume_acceleration = r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        r_volume_acceleration[0] = 0.0;
        r_volume_acceleration[1] = -9.81;
        r_volume_acceleration[2] = 0.0;

        ++node_index;
    }

    KRATOS_CATCH("");
}



/// Tolerance associated with a reference value (same philosophy as before).
double ComponentTolerance(const double rReferenceValue, const double rReferenceScale)
{
    double tolerance = std::max(tle_comparison_absolute_tolerance,
                                tle_comparison_relative_tolerance * std::abs(rReferenceValue));
    if (rReferenceValue == 0.0) {
        tolerance = std::max(tolerance, tle_machine_precision_allowance * rReferenceScale);
    }
    return tolerance;
}

double MaxAbsEntry(const Matrix& rMatrix)
{
    double max_abs_entry = 0.0;
    for (std::size_t i = 0; i < rMatrix.size1(); ++i) {
        for (std::size_t j = 0; j < rMatrix.size2(); ++j) {
            max_abs_entry = std::max(max_abs_entry, std::abs(rMatrix(i, j)));
        }
    }
    return max_abs_entry;
}

double MaxAbsEntry(const Vector& rVector)
{
    double max_abs_entry = 0.0;
    for (std::size_t i = 0; i < rVector.size(); ++i) {
        max_abs_entry = std::max(max_abs_entry, std::abs(rVector(i)));
    }
    return max_abs_entry;
}

struct ComparisonMetrics
{
    bool pass = true;
    std::size_t component_count = 0;
    std::size_t failed_component_count = 0;
    double max_absolute_difference = 0.0;
    double max_relative_difference = 0.0;
    double max_tolerance_used = 0.0;
    double reference_scale = 0.0;
};

void AccumulateComponentMetrics(
    const double rComputedValue,
    const double rReferenceValue,
    const double rTolerance,
    ComparisonMetrics& rMetrics)
{
    const double absolute_difference = std::abs(rComputedValue - rReferenceValue);
    rMetrics.max_absolute_difference = std::max(rMetrics.max_absolute_difference, absolute_difference);
    rMetrics.max_tolerance_used = std::max(rMetrics.max_tolerance_used, rTolerance);
    rMetrics.reference_scale = std::max(rMetrics.reference_scale, std::abs(rReferenceValue));
    if (rReferenceValue != 0.0) {
        rMetrics.max_relative_difference =
            std::max(rMetrics.max_relative_difference, absolute_difference / std::abs(rReferenceValue));
    }
    if (absolute_difference > rTolerance) {
        ++rMetrics.failed_component_count;
        rMetrics.pass = false;
    }
    ++rMetrics.component_count;
}

void PrintComparisonMetrics(const std::string& rWhat, const ComparisonMetrics& rMetrics)
{
    std::cout << "[CHARACTERIZATION] " << rWhat << ": "
              << (rMetrics.pass ? "PASS" : "FAIL")
              << " | max_abs_diff=" << rMetrics.max_absolute_difference
              << " | max_rel_diff=" << rMetrics.max_relative_difference
              << " | tolerance=" << rMetrics.max_tolerance_used
              << " | scale=" << rMetrics.reference_scale
              << " | failed_components=" << rMetrics.failed_component_count
              << "/" << rMetrics.component_count
              << std::endl;
}

void ExpectVectorComponentsNear(
    const Vector& rComputed,
    const Vector& rReference,
    const std::string& rWhat)
{
    ASSERT_EQ(rComputed.size(), rReference.size()) << "Size mismatch comparing " << rWhat;
    const double reference_scale = MaxAbsEntry(rReference);
    ComparisonMetrics metrics;
    for (std::size_t i = 0; i < rReference.size(); ++i) {
        AccumulateComponentMetrics(
            rComputed(i), rReference(i), ComponentTolerance(rReference(i), reference_scale), metrics);
    }
    PrintComparisonMetrics(rWhat, metrics);
    for (std::size_t i = 0; i < rReference.size(); ++i) {
        KRATOS_EXPECT_NEAR(rComputed(i), rReference(i), ComponentTolerance(rReference(i), reference_scale));
    }
}

void ExpectMatrixComponentsNear(
    const Matrix& rComputed,
    const Matrix& rReference,
    const std::string& rWhat)
{
    ASSERT_EQ(rComputed.size1(), rReference.size1()) << "Row size mismatch comparing " << rWhat;
    ASSERT_EQ(rComputed.size2(), rReference.size2()) << "Column size mismatch comparing " << rWhat;
    const double reference_scale = MaxAbsEntry(rReference);
    ComparisonMetrics metrics;
    for (std::size_t i = 0; i < rReference.size1(); ++i) {
        for (std::size_t j = 0; j < rReference.size2(); ++j) {
            AccumulateComponentMetrics(
                rComputed(i, j), rReference(i, j), ComponentTolerance(rReference(i, j), reference_scale), metrics);
        }
    }
    PrintComparisonMetrics(rWhat, metrics);
    for (std::size_t i = 0; i < rReference.size1(); ++i) {
        for (std::size_t j = 0; j < rReference.size2(); ++j) {
            KRATOS_EXPECT_NEAR(rComputed(i, j), rReference(i, j), ComponentTolerance(rReference(i, j), reference_scale));
        }
    }
}







/// Owns the vectors bound to a ConstitutiveLaw::Parameters (the Parameters
/// stores pointers to them, so they must outlive the response evaluation).
struct FamilyLawParametersBundle
{
    Vector strain;
    Vector stress;
    Matrix constitutive_matrix;
    Vector shape_function_values;
    ConstitutiveLaw::Parameters values;

    FamilyLawParametersBundle(
        const Geometry<Node>& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo,
        const Vector& rStrain)
        : strain(rStrain),
          stress(rStrain.size()),
          constitutive_matrix(rStrain.size(), rStrain.size()),
          shape_function_values(rGeometry.PointsNumber()),
          values(rGeometry, rProperties, rProcessInfo)
    {
        noalias(stress) = ZeroVector(rStrain.size());
        noalias(constitutive_matrix) = ZeroMatrix(rStrain.size(), rStrain.size());
        noalias(shape_function_values) = row(rGeometry.ShapeFunctionsValues(), 0);

        auto& r_options = values.GetOptions();
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        values.SetStrainVector(strain);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(constitutive_matrix);
        values.SetShapeFunctionsValues(shape_function_values);
    }
};

/// The three Dam thermal linear-elastic families covered by these tests.
const std::vector<std::string>& TleFamilyLawNames()
{
    static const std::vector<std::string> law_names = {
        "ThermalLinearElastic3DLaw",
        "ThermalLinearElastic2DPlaneStrain",
        "ThermalLinearElastic2DPlaneStress"};
    return law_names;
}

bool TleIs3D(const std::string& rLawName)
{
    return rLawName == "ThermalLinearElastic3DLaw";
}

std::size_t TleFamilyStrainSize(const std::string& rLawName)
{
    return TleIs3D(rLawName) ? 6 : 3;
}

std::string TleFamilyElementName(const std::string& rLawName)
{
    return TleIs3D(rLawName) ? "SmallDisplacementElement3D4N"
                             : "SmallDisplacementElement2D3N";
}

/// Deterministic, fully populated non-trivial total strain state.
Vector TleFamilyTestStrain(const std::string& rLawName)
{
    const std::size_t strain_size = TleFamilyStrainSize(rLawName);
    Vector strain(strain_size);
    strain[0] = 1.0e-3;
    strain[1] = -5.0e-4;
    strain[2] = 2.0e-4;
    if (strain_size > 3) {
        strain[3] = 3.0e-4;
        strain[4] = -2.0e-4;
        strain[5] = 1.0e-4;
    }
    return strain;
}

/// The mechanical (zero thermal strain) reference law of a Dam thermal family:
/// same Young's modulus, Poisson ratio and strain state, no thermal strain.
ConstitutiveLaw::Pointer tle_CreateMechanicalReferenceLaw(const std::string& rLawName)
{
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        return ConstitutiveLaw::Pointer(new LinearPlaneStrain());
    } else if (rLawName == "ThermalLinearElastic2DPlaneStress") {
        return ConstitutiveLaw::Pointer(new LinearPlaneStress());
    }
    return ConstitutiveLaw::Pointer(new FlexibleElasticIsotropic3D());
}

/// Evaluates the PK2 response (stress vector + constitutive matrix) of a law
/// for the given total strain on the family geometry/properties/process info.
void tle_EvaluatePK2Response(
    ConstitutiveLaw::Pointer pLaw,
    const Geometry<Node>& rGeometry,
    const Properties& rProperties,
    const ProcessInfo& rProcessInfo,
    const Vector& rStrain,
    Vector& rStress,
    Matrix& rConstitutiveMatrix)
{
    FamilyLawParametersBundle values(rGeometry, rProperties, rProcessInfo, rStrain);
    Vector shape_function_values(rGeometry.PointsNumber());
    noalias(shape_function_values) = row(rGeometry.ShapeFunctionsValues(), 0);
    pLaw->InitializeMaterial(rProperties, rGeometry, shape_function_values);
    pLaw->CalculateMaterialResponse(values.values, ConstitutiveLaw::StressMeasure_PK2);
    rStress = values.values.GetStressVector();
    rConstitutiveMatrix = values.values.GetConstitutiveMatrix();
}

} // namespace

//************************************************************************************
// 2D plane strain
//************************************************************************************



//************************************************************************************
// 2D plane stress
//************************************************************************************



//************************************************************************************
// 3D tetrahedron
//************************************************************************************


//************************************************************************************
// Analytical restrained thermal stress
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalFamily_Analytical_PlaneStrainRestrained, KratosDamFastSuite)
{
    // Restrained uniform thermal expansion (T1): zero total strain, uniform
    // delta_T = 25. The plane-strain law applies the effective thermal strain
    // (1+nu)*alpha*delta_T*[1,1,0], which must give the 3D-consistent in-plane
    // restrained stress sigma_xx = sigma_yy = -E*alpha*delta_T/(1-2*nu).
    const double expected_stress =
        -tle_test_young_modulus * tle_test_thermal_expansion * test_delta_temperature / (1.0 - 2.0 * tle_test_poisson_ratio);
    Vector expected(3);
    expected[0] = expected_stress;
    expected[1] = expected_stress;
    expected[2] = 0.0;

    Model model;
    ModelPart& r_model_part = CreateFamilyModelPart(
        model, "AnalyticalPlaneStrain", "SmallDisplacementElement2D3N",
        "ThermalLinearElastic2DPlaneStrain", 2, 0);
    PrescribeFamilyScenario(r_model_part, 1, 2);
    auto p_element = r_model_part.pGetElement(1);
    const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    KRATOS_EXPECT_EQ(p_element->Check(r_pi), 0);
    p_element->Initialize(r_pi);
    p_element->InitializeSolutionStep(r_pi);
    p_element->InitializeNonLinearIteration(r_pi);

    std::vector<Vector> cauchy;
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy, r_pi);
    for (std::size_t gp = 0; gp < cauchy.size(); ++gp) {
        ExpectVectorComponentsNear(cauchy[gp], expected, "plane strain analytical restrained stress");
    }
}


KRATOS_TEST_CASE_IN_SUITE(ThermalFamily_Analytical_PlaneStressRestrained, KratosDamFastSuite)
{
    // Restrained uniform thermal expansion (T1): zero total strain, uniform
    // delta_T = 25. The plane-stress law applies epsilon_th = alpha*delta_T*[1,1,0]
    // giving sigma_xx = sigma_yy = -E*alpha*delta_T/(1-nu).
    const double expected_stress =
        -tle_test_young_modulus * tle_test_thermal_expansion * test_delta_temperature / (1.0 - tle_test_poisson_ratio);
    Vector expected(3);
    expected[0] = expected_stress;
    expected[1] = expected_stress;
    expected[2] = 0.0;

    Model model;
    ModelPart& r_model_part = CreateFamilyModelPart(
        model, "AnalyticalPlaneStress", "SmallDisplacementElement2D3N",
        "ThermalLinearElastic2DPlaneStress", 2, 0);
    PrescribeFamilyScenario(r_model_part, 1, 2);
    auto p_element = r_model_part.pGetElement(1);
    const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    KRATOS_EXPECT_EQ(p_element->Check(r_pi), 0);
    p_element->Initialize(r_pi);
    p_element->InitializeSolutionStep(r_pi);
    p_element->InitializeNonLinearIteration(r_pi);

    std::vector<Vector> cauchy;
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy, r_pi);
    for (std::size_t gp = 0; gp < cauchy.size(); ++gp) {
        ExpectVectorComponentsNear(cauchy[gp], expected, "plane stress analytical restrained stress");
    }
}


KRATOS_TEST_CASE_IN_SUITE(ThermalFamily_Analytical_3DRestrained, KratosDamFastSuite)
{
    // Restrained uniform thermal expansion (T1): zero total strain, uniform
    // delta_T = 25. The 3D law applies epsilon_th = alpha*delta_T*[1,1,1,0,0,0],
    // so epsilon_mech = -epsilon_th and the isotropic confined stress is
    //   sigma_xx = sigma_yy = sigma_zz = -E*alpha*delta_T/(1 - 2*nu),
    // sigma_xy = sigma_yz = sigma_xz = 0  (since 3*lambda + 2*mu = E/(1-2*nu)).
    const double expected_normal_stress =
        -tle_test_young_modulus * tle_test_thermal_expansion * test_delta_temperature /
        (1.0 - 2.0 * tle_test_poisson_ratio);
    Vector expected(6);
    expected[0] = expected_normal_stress;
    expected[1] = expected_normal_stress;
    expected[2] = expected_normal_stress;
    expected[3] = 0.0;
    expected[4] = 0.0;
    expected[5] = 0.0;

    Model model;
    ModelPart& r_model_part = CreateFamilyModelPart(
        model, "Analytical3D", "SmallDisplacementElement3D4N",
        "ThermalLinearElastic3DLaw", 3, 0);
    PrescribeFamilyScenario(r_model_part, 1, 3);
    auto p_element = r_model_part.pGetElement(1);
    const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    KRATOS_EXPECT_EQ(p_element->Check(r_pi), 0);
    p_element->Initialize(r_pi);
    p_element->InitializeSolutionStep(r_pi);
    p_element->InitializeNonLinearIteration(r_pi);

    std::vector<Vector> cauchy;
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy, r_pi);
    KRATOS_EXPECT_FALSE(cauchy.empty());
    for (std::size_t gp = 0; gp < cauchy.size(); ++gp) {
        ExpectVectorComponentsNear(cauchy[gp], expected, "3D analytical restrained stress");
    }
}


//************************************************************************************
// 2D lifecycle (inherited behavior)
//************************************************************************************



//************************************************************************************
// Serialization of the 2D laws (chains through the modified 3D base)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalFamily_Serialization_2DLaws, KratosDamFastSuite)
{
    // Serialize and deserialize the plane-strain and plane-stress laws with the
    // in-memory StreamSerializer. Their save/load chain through
    // ThermalLinearElastic3DLaw (the modified base). Verify that the
    // deserialized laws reproduce the same thermo-elastic response.
    for (const char* p_law_name : {"ThermalLinearElastic2DPlaneStrain",
                                    "ThermalLinearElastic2DPlaneStress"}) {
        const std::string r_law_name(p_law_name);
        Model model;
        ModelPart& r_model_part = CreateFamilyModelPart(
            model, "Serialization2D" + r_law_name, "SmallDisplacementElement2D3N",
            r_law_name, 2, 0);
        auto p_element = r_model_part.pGetElement(1);
        const auto& r_geometry = p_element->GetGeometry();
        const Properties& r_properties = p_element->GetProperties();
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

        ConstitutiveLaw::Pointer p_law = r_properties.GetValue(CONSTITUTIVE_LAW);
        Vector shape_function_values(r_geometry.PointsNumber());
        noalias(shape_function_values) = row(r_geometry.ShapeFunctionsValues(), 0);
        p_law->InitializeMaterial(r_properties, r_geometry, shape_function_values);

        PrescribeFamilyScenario(r_model_part, 1, 2);
        Vector strain(3);
        noalias(strain) = ZeroVector(3);
        strain[0] = 1.0e-3;
        strain[1] = -5.0e-4;
        strain[2] = 2.0e-4;

        FamilyLawParametersBundle values_before(r_geometry, r_properties, r_pi, strain);
        p_law->CalculateMaterialResponse(values_before.values, ConstitutiveLaw::StressMeasure_PK2);

        StreamSerializer serializer;
        serializer.save("Law", p_law);
        ConstitutiveLaw::Pointer p_loaded;
        serializer.load("Law", p_loaded);
        KRATOS_EXPECT_TRUE(p_loaded != nullptr);

        FamilyLawParametersBundle values_after(r_geometry, r_properties, r_pi, strain);
        p_loaded->CalculateMaterialResponse(values_after.values, ConstitutiveLaw::StressMeasure_PK2);

        ExpectVectorComponentsNear(
            values_before.values.GetStressVector(), values_after.values.GetStressVector(),
            r_law_name + " serialization stress");
        ExpectMatrixComponentsNear(
            values_before.values.GetConstitutiveMatrix(), values_after.values.GetConstitutiveMatrix(),
            r_law_name + " serialization constitutive matrix");
    }
}


//************************************************************************************
// Zero temperature increment: thermal adapter == underlying mechanical behavior
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalFamily_ZeroDeltaTemperatureEquivalence, KratosDamFastSuite)
{
    // With Delta_T = T - T_ref = 0 (nodal temperature == nodal reference
    // temperature, as set by CreateFamilyModelPart) the thermal strain vanishes
    // identically. Each Dam thermal adapter must then reproduce the underlying
    // mechanical law exactly for the same Young's modulus, Poisson ratio and
    // strain state, and report null thermal outputs. The comparison is repeated
    // with the standard nodal Young's-modulus accessor configuration (uniform
    // nodal field equal to the property value), where applicable.
    for (const std::string& r_law_name : TleFamilyLawNames()) {
        Model model;
        ModelPart& r_model_part = CreateFamilyModelPart(
            model, "ZeroDeltaT" + r_law_name, TleFamilyElementName(r_law_name),
            r_law_name, TleIs3D(r_law_name) ? 3 : 2, 0);
        auto p_element = r_model_part.pGetElement(1);
        const auto& r_geometry = p_element->GetGeometry();
        Properties& r_properties = p_element->GetProperties();
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        const Vector strain = TleFamilyTestStrain(r_law_name);

        ConstitutiveLaw::Pointer p_thermal_law = r_properties.GetValue(CONSTITUTIVE_LAW);
        ConstitutiveLaw::Pointer p_reference_law = tle_CreateMechanicalReferenceLaw(r_law_name);

        for (const bool use_accessor : {false, true}) {
            if (use_accessor) {
                for (auto& r_node : r_model_part.Nodes()) {
                    r_node.FastGetSolutionStepValue(YOUNG_MODULUS) = tle_test_young_modulus;
                }
                NodalYoungModulusUtilities::InstallDatabaseAccessor(r_properties);
            }
            const std::string label = r_law_name + " zero-dT" +
                                      (use_accessor ? " [nodal-Y accessor]" : "");

            Vector thermal_stress, reference_stress;
            Matrix thermal_matrix, reference_matrix;
            tle_EvaluatePK2Response(p_thermal_law, r_geometry, r_properties, r_pi,
                                    strain, thermal_stress, thermal_matrix);
            tle_EvaluatePK2Response(p_reference_law, r_geometry, r_properties, r_pi,
                                    strain, reference_stress, reference_matrix);
            ExpectVectorComponentsNear(
                thermal_stress, reference_stress, label + " stress vs mechanical reference");
            ExpectMatrixComponentsNear(
                thermal_matrix, reference_matrix,
                label + " constitutive matrix vs mechanical reference");

            // Specialized Dam outputs must show no thermal contribution.
            FamilyLawParametersBundle values(r_geometry, r_properties, r_pi, strain);
            p_thermal_law->CalculateMaterialResponse(values.values, ConstitutiveLaw::StressMeasure_PK2);
            Vector thermal_strain_out, thermal_stress_out, mechanical_stress_out;
            p_thermal_law->CalculateValue(values.values, THERMAL_STRAIN_VECTOR, thermal_strain_out);
            p_thermal_law->CalculateValue(values.values, THERMAL_STRESS_VECTOR, thermal_stress_out);
            p_thermal_law->CalculateValue(values.values, MECHANICAL_STRESS_VECTOR, mechanical_stress_out);

            Vector zero(strain.size());
            noalias(zero) = ZeroVector(strain.size());
            ExpectVectorComponentsNear(thermal_strain_out, zero,
                                       label + " THERMAL_STRAIN_VECTOR");
            ExpectVectorComponentsNear(thermal_stress_out, zero,
                                       label + " THERMAL_STRESS_VECTOR");
            ExpectVectorComponentsNear(mechanical_stress_out, thermal_stress,
                                       label + " MECHANICAL_STRESS_VECTOR == total stress");
        }
    }
}


//************************************************************************************
// CalculateValue: non-Dam variables must reach the CLA/SMA base implementation
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalFamily_CalculateValueBaseFallback, KratosDamFastSuite)
{
    // Variables that are not Dam-specific outputs are served by the CLA/SMA
    // base-class CalculateValue implementations (ElasticIsotropic3D): the
    // adapters must delegate to ElasticIsotropic3D::CalculateValue instead of
    // returning the caller's rValue untouched. Sentinels detect the shadowing
    // for both the Vector and the Matrix overload independently.
    //
    // The state is a genuinely thermal one (scenario 1: uniform nodal
    // T = T_ref + 25 = 45 with T_ref = 20, so Delta_T = +25 for every law).
    // CAUCHY_STRESS_VECTOR therefore only matches the expected value if the
    // delegation reaches the Dam thermal response through the virtual chain
    // CalculateMaterialResponsePK2 -> SubstractThermalStrain; a purely
    // mechanical C*epsilon answer would omit that term.
    for (const std::string& r_law_name : TleFamilyLawNames()) {
        Model model;
        ModelPart& r_model_part = CreateFamilyModelPart(
            model, "BaseFallback" + r_law_name, TleFamilyElementName(r_law_name),
            r_law_name, TleIs3D(r_law_name) ? 3 : 2, 0);
        PrescribeFamilyScenario(r_model_part, 1, TleIs3D(r_law_name) ? 3 : 2);
        auto p_element = r_model_part.pGetElement(1);
        const auto& r_geometry = p_element->GetGeometry();
        const Properties& r_properties = p_element->GetProperties();
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        const Vector strain = TleFamilyTestStrain(r_law_name);
        ConstitutiveLaw::Pointer p_law = r_properties.GetValue(CONSTITUTIVE_LAW);

        // Reference thermal response: element-provided strain plus exactly one
        // in-place thermal subtraction (CalculateMaterialResponsePK2 subtracts
        // the thermal strain on the Parameters strain vector on every call, so
        // every assertion below runs on its own pristine Parameters).
        FamilyLawParametersBundle expected_values(r_geometry, r_properties, r_pi, strain);
        p_law->CalculateMaterialResponse(expected_values.values, ConstitutiveLaw::StressMeasure_PK2);
        const Vector expected_stress = expected_values.values.GetStressVector();
        const Matrix expected_matrix = expected_values.values.GetConstitutiveMatrix();

        // Guard: the expected stress must carry a genuine thermal contribution,
        // otherwise the CAUCHY_STRESS_VECTOR comparison below would be
        // satisfied by a purely mechanical law and prove nothing about the
        // dispatch chain.
        const Vector mechanical_reference = prod(expected_matrix, strain);
        KRATOS_EXPECT_GT(
            norm_2(expected_stress - mechanical_reference),
            0.1 * norm_2(mechanical_reference));

        // Fresh Parameters for the fallback assertions: the base recomputes the
        // response itself, so the element-provided strain is still intact.
        FamilyLawParametersBundle values(r_geometry, r_properties, r_pi, strain);

        // Vector overload, strain variable: base law returns the strain state.
        Vector strain_out(1);
        strain_out[0] = 1.0e30;
        p_law->CalculateValue(values.values, STRAIN, strain_out);
        ExpectVectorComponentsNear(
            strain_out, strain,
            r_law_name + " CalculateValue(STRAIN) base fallback");

        // Matrix overload: base law assembles the constitutive matrix.
        Matrix constitutive_out(1, 1);
        constitutive_out(0, 0) = 1.0e30;
        p_law->CalculateValue(values.values, CONSTITUTIVE_MATRIX, constitutive_out);
        ExpectMatrixComponentsNear(
            constitutive_out, expected_matrix,
            r_law_name + " CalculateValue(CONSTITUTIVE_MATRIX) base fallback");

        // Vector overload, stress variable: base law recomputes the thermal
        // response through the virtual chain. Evaluated last because it
        // subtracts the thermal strain in place on this Parameters strain
        // vector, and compared against the independently obtained reference
        // instead of the Parameters stress buffer it just wrote.
        Vector stress_out(1);
        stress_out[0] = 1.0e30;
        p_law->CalculateValue(values.values, CAUCHY_STRESS_VECTOR, stress_out);
        ExpectVectorComponentsNear(
            stress_out, expected_stress,
            r_law_name + " CalculateValue(CAUCHY_STRESS_VECTOR) base fallback");
    }
}


} // namespace Testing
} // namespace Kratos
