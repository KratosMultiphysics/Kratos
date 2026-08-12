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

// Phase 3D.2 characterization: the new
// DamNodalCauchyStressExtrapolationProcess must reproduce exactly the nodal
// Cauchy-stress recovery currently embedded in
// SmallDisplacementThermoMechanicElement::FinalizeSolutionStep, while operating
// on the standard StructuralMechanicsApplication small-displacement elements.
//
// The comparison uses the RAW historical accumulators:
//   NODAL_CAUCHY_STRESS_TENSOR (sum of element_measure * extrapolated_stress)
//   NODAL_AREA                 (sum of element_measure)
// matching the state produced by the legacy element before the Dam smoothing
// scheme normalizes NODAL_CAUCHY_STRESS_TENSOR by NODAL_AREA.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <limits>
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
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_processes/dam_nodal_cauchy_stress_extrapolation_process.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances (same philosophy as the previous characterization).
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;
constexpr double machine_precision_allowance = 1.0e-15;

/// Material data shared by all model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_thickness = 0.15;

/// Creates the constitutive law for the given configuration.
ConstitutiveLaw::Pointer CreateThermalLaw(const std::string& rLawName)
{
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        return ThermalLinearElastic2DPlaneStrain().Clone();
    } else {
        return ThermalLinearElastic3DLaw().Clone();
    }
}

/// Creates a model part with the given element and the thermal law, with the
/// nodal accumulator variables initialized (dimension x dimension) to zero.
ModelPart& CreateProcessModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rLawName,
    const std::size_t rDimension)
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
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

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
    if (rDimension == 2) {
        (*p_prop)[THICKNESS] = test_thickness;
    }
    p_prop->SetValue(CONSTITUTIVE_LAW, CreateThermalLaw(rLawName));

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
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(rDimension, rDimension);
        noalias(zero_stress_tensor) = ZeroMatrix(rDimension, rDimension);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Prescribes a thermo-mechanical state with a non-uniform Gauss-point stress.
void PrescribeVaryingState(ModelPart& rModelPart, const std::size_t rDimension)
{
    KRATOS_TRY;

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = 1.0e-3 * r_x0[0] + 3.0e-4 * r_x0[1] + 1.0e-4 * r_x0[0] * r_x0[1]
                           + 2.0e-4 * r_x0[0] * r_x0[2];
        r_displacement[1] = 3.0e-4 * r_x0[0] - 1.0e-3 * r_x0[1] + 2.0e-4 * r_x0[0] * r_x0[0]
                           + 1.0e-4 * r_x0[1] * r_x0[2];
        r_displacement[2] = 5.0e-4 * r_x0[2] + 2.0e-4 * r_x0[0] + 1.0e-4 * r_x0[1]
                           + 1.0e-4 * r_x0[0] * r_x0[1];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + 5.0 * static_cast<double>(node_index);
        ++node_index;
    }

    KRATOS_CATCH("");
}

/// Prescribes a thermo-mechanical state whose yz/xz shear strains are
/// theoretically zero (the field used in the earlier diagnostic). It is used only
/// to characterize the upstream element-level roundoff in those components.
void PrescribeZeroShearState(ModelPart& rModelPart, const std::size_t rDimension)
{
    KRATOS_TRY;

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = 1.0e-3 * r_x0[0] + 3.0e-4 * r_x0[1] + 1.0e-4 * r_x0[0] * r_x0[1];
        r_displacement[1] = 3.0e-4 * r_x0[0] - 1.0e-3 * r_x0[1] + 2.0e-4 * r_x0[0] * r_x0[0];
        r_displacement[2] = 5.0e-4 * r_x0[2];
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + 5.0 * static_cast<double>(node_index);
        ++node_index;
    }

    KRATOS_CATCH("");
}

/// Initializes an element through the standard lifecycle.
void InitializeProcessElement(ModelPart& rModelPart)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    KRATOS_EXPECT_EQ(p_element->Check(r_process_info), 0);
    p_element->Initialize(r_process_info);
    p_element->InitializeSolutionStep(r_process_info);
    p_element->InitializeNonLinearIteration(r_process_info);

    KRATOS_CATCH("");
}

/// Runs the legacy element finalization (accumulates the raw nodal values).

/// Runs the new Dam process (accumulation only; the caller must have
/// initialized the nodal accumulators, as the legacy smoothing scheme does).
void RunCandidateProcess(ModelPart& rModelPart)
{
    DamNodalCauchyStressExtrapolationProcess process(rModelPart);
    process.ExecuteFinalizeSolutionStep();
}

/// Replicates the legacy scheme reset: zeroes NODAL_AREA and
/// NODAL_CAUCHY_STRESS_TENSOR on all nodes.
void ResetNodalAccumulators(ModelPart& rModelPart, const std::size_t rDimension)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        noalias(r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) = ZeroMatrix(rDimension, rDimension);
    }
}

/// Replicates the legacy scheme normalization: divides
/// NODAL_CAUCHY_STRESS_TENSOR by NODAL_AREA when the area is non-negligible.

/// Compares raw historical NODAL_CAUCHY_STRESS_TENSOR and NODAL_AREA.

/// Minimal test-only element that returns a prescribed CAUCHY_STRESS_TENSOR
/// field at its integration points. It is NOT registered in KratosComponents and
/// introduces no production dependency; it is used to isolate the process-level
/// extrapolation from the element-level strain computation.
class PrescribedStressElement : public Element
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrescribedStressElement);
    using BaseType = Element;

    PrescribedStressElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry) {}

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == CAUCHY_STRESS_TENSOR) {
            rOutput = mPrescribedStress;
        } else {
            KRATOS_ERROR << "PrescribedStressElement: unsupported variable " << rVariable << std::endl;
        }
    }

    void SetPrescribedStress(std::vector<Matrix> rStress) { mPrescribedStress = std::move(rStress); }

private:
    std::vector<Matrix> mPrescribedStress;
};

/// Creates a model part with the prescribed-stress element and the accumulator
/// nodal variables initialized to zero. Returns the element pointer through
/// rOutElement.
ModelPart& CreatePrescribedGPModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rPrototypeElementName,
    const std::size_t rDimension,
    PrescribedStressElement*& rOutElement)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;

    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    const Element& r_prototype = KratosComponents<Element>::Get(rPrototypeElementName);
    const auto& r_prototype_geometry = r_prototype.GetGeometry();
    Matrix local_coordinates;
    r_prototype_geometry.PointsLocalCoordinates(local_coordinates);

    const double geometry_scale = 2.5;
    array_1d<double, 3> geometry_offset;
    geometry_offset[0] = 0.75;
    geometry_offset[1] = 1.25;
    geometry_offset[2] = (rDimension == 3) ? 0.5 : 0.0;

    const std::size_t number_of_nodes = r_prototype_geometry.PointsNumber();
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const double x = geometry_scale * local_coordinates(i, 0) + geometry_offset[0];
        const double y = geometry_scale * local_coordinates(i, 1) + geometry_offset[1];
        const double z = (rDimension == 3) ? geometry_scale * local_coordinates(i, 2) + geometry_offset[2] : 0.0;
        r_model_part.CreateNewNode(i + 1, x, y, z);
    }
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(rDimension, rDimension);
        noalias(zero_stress_tensor) = ZeroMatrix(rDimension, rDimension);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    Geometry<Node>::PointsArrayType points;
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        points.push_back(r_model_part.pGetNode(i + 1));
    }
    Geometry<Node>::Pointer p_geometry = r_prototype_geometry.Create(points);
    auto p_element = Kratos::make_intrusive<PrescribedStressElement>(1, p_geometry);
    r_model_part.AddElement(p_element);
    rOutElement = p_element.get();

    return r_model_part;

    KRATOS_CATCH("");
}

/// Computes the expected legacy nodal raw accumulator for a prescribed GP
/// stress field, following exactly the legacy mathematical path.
Matrix ComputeExpectedNodalRawStress(
    const Geometry<Node>& rGeometry,
    const std::vector<Matrix>& rGaussStress,
    const std::size_t rDimension,
    const std::size_t rNodeIndex)
{
    const std::size_t number_of_gps = rGaussStress.size();
    const std::size_t voigt_size = (rDimension == 2) ? 3 : 6;
    Matrix stress_container(number_of_gps, voigt_size);
    noalias(stress_container) = ZeroMatrix(number_of_gps, voigt_size);
    for (std::size_t gp = 0; gp < number_of_gps; ++gp) {
        stress_container(gp, 0) = rGaussStress[gp](0, 0);
        stress_container(gp, 1) = rGaussStress[gp](1, 1);
        if (rDimension == 2) {
            stress_container(gp, 2) = rGaussStress[gp](0, 1);
        } else {
            stress_container(gp, 2) = rGaussStress[gp](2, 2);
            stress_container(gp, 3) = rGaussStress[gp](0, 1);
            stress_container(gp, 4) = rGaussStress[gp](1, 2);
            stress_container(gp, 5) = rGaussStress[gp](0, 2);
        }
    }

    Matrix aux_nodal_stress;
    if (rDimension == 2) {
        BoundedMatrix<double, 4, 4> q4;
        PoroElementUtilities::Calculate2DExtrapolationMatrix(q4);
        aux_nodal_stress = prod(Matrix(q4), stress_container);
    } else {
        BoundedMatrix<double, 8, 8> h8;
        PoroElementUtilities::Calculate3DExtrapolationMatrix(h8);
        aux_nodal_stress = prod(Matrix(h8), stress_container);
    }

    Vector nodal_voigt(voigt_size);
    for (std::size_t comp = 0; comp < voigt_size; ++comp) {
        nodal_voigt(comp) = aux_nodal_stress(rNodeIndex, comp);
    }
    const Matrix nodal_tensor = MathUtils<double>::StressVectorToTensor(nodal_voigt);
    return rGeometry.Area() * nodal_tensor;
}

/// Verifies the process result against the independently computed analytic
/// legacy formula: NODAL_CAUCHY_STRESS_TENSOR = sum over incident elements of
/// (element_measure * extrapolated_stress), NODAL_AREA = sum of measures.
/// After Phase 3D.4B the process is the single production implementation and
/// there is no element-level reference to compare against.
void VerifyProcessAgainstAnalyticRaw(ModelPart& rModelPart, const std::string& rLabel)
{
    const std::size_t dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
    std::map<ModelPart::IndexType, Matrix> expected_stress;
    std::map<ModelPart::IndexType, double> expected_area;
    for (auto& r_element : rModelPart.Elements()) {
        if (!r_element.IsActive()) continue;
        std::vector<Matrix> gauss_stress;
        r_element.CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, gauss_stress, rModelPart.GetProcessInfo());
        const auto& r_geom = r_element.GetGeometry();
        const std::size_t n_nodes = r_geom.size();
        const std::size_t n_gps = gauss_stress.size();
        const std::size_t voigt = (dim == 2) ? 3 : 6;
        Matrix container(n_gps, voigt);
        noalias(container) = ZeroMatrix(n_gps, voigt);
        for (std::size_t gp = 0; gp < n_gps; ++gp) {
            container(gp, 0) = gauss_stress[gp](0, 0);
            container(gp, 1) = gauss_stress[gp](1, 1);
            if (dim == 2) container(gp, 2) = gauss_stress[gp](0, 1);
            else {
                container(gp, 2) = gauss_stress[gp](2, 2);
                container(gp, 3) = gauss_stress[gp](0, 1);
                container(gp, 4) = gauss_stress[gp](1, 2);
                container(gp, 5) = gauss_stress[gp](0, 2);
            }
        }
        Matrix aux;
        if (n_gps == 1) {
            aux = ZeroMatrix(n_nodes, voigt);
            for (std::size_t i = 0; i < n_nodes; ++i)
                for (std::size_t c = 0; c < voigt; ++c) aux(i, c) = container(0, c);
        } else if (dim == 2) {
            BoundedMatrix<double, 4, 4> q4;
            PoroElementUtilities::Calculate2DExtrapolationMatrix(q4);
            aux = prod(Matrix(q4), container);
        } else {
            BoundedMatrix<double, 8, 8> h8;
            PoroElementUtilities::Calculate3DExtrapolationMatrix(h8);
            aux = prod(Matrix(h8), container);
        }
        const double area = r_geom.Area();
        for (std::size_t i = 0; i < n_nodes; ++i) {
            Vector nv(voigt);
            for (std::size_t c = 0; c < voigt; ++c) nv(c) = aux(i, c);
            Matrix contrib = area * MathUtils<double>::StressVectorToTensor(nv);
            const ModelPart::IndexType nid = r_geom[i].Id();
            if (expected_stress.find(nid) == expected_stress.end()) {
                expected_stress[nid] = ZeroMatrix(dim, dim);
                expected_area[nid] = 0.0;
            }
            expected_stress[nid] += contrib;
            expected_area[nid] += area;
        }
    }
    for (auto& r_node : rModelPart.Nodes()) {
        const ModelPart::IndexType nid = r_node.Id();
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(NODAL_AREA), expected_area[nid],
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(expected_area[nid])));
        const Matrix& stored = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < dim; ++i) {
            for (std::size_t j = 0; j < dim; ++j) {
                KRATOS_EXPECT_NEAR(stored(i, j), expected_stress[nid](i, j),
                                   std::max(comparison_absolute_tolerance,
                                            comparison_relative_tolerance * std::abs(expected_stress[nid](i, j))));
            }
        }
    }
}

} // namespace

//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_SingleElement_Triangle2D3, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "CandidateT3", "SmallDisplacementElement2D3N",
        "ThermalLinearElastic2DPlaneStrain", 2);
    PrescribeVaryingState(r_candidate, 2);
    InitializeProcessElement(r_candidate);
    RunCandidateProcess(r_candidate);
    VerifyProcessAgainstAnalyticRaw(r_candidate, "Triangle2D3");
}

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_SingleElement_Quadrilateral2D4, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "CandidateQ4", "SmallDisplacementElement2D4N",
        "ThermalLinearElastic2DPlaneStrain", 2);
    PrescribeVaryingState(r_candidate, 2);
    InitializeProcessElement(r_candidate);
    RunCandidateProcess(r_candidate);
    VerifyProcessAgainstAnalyticRaw(r_candidate, "Quadrilateral2D4");
}

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_SingleElement_Tetrahedra3D4, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "CandidateT4", "SmallDisplacementElement3D4N",
        "ThermalLinearElastic3DLaw", 3);
    PrescribeVaryingState(r_candidate, 3);
    InitializeProcessElement(r_candidate);
    RunCandidateProcess(r_candidate);
    VerifyProcessAgainstAnalyticRaw(r_candidate, "Tetrahedra3D4");
}

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_SingleElement_Hexahedra3D8, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "CandidateH8", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3);
    PrescribeVaryingState(r_candidate, 3);
    InitializeProcessElement(r_candidate);

    // Verify that the GP stress field is non-uniform, all six components are
    // non-degenerate, and the legacy and candidate integration-point stresses
    // agree within the established mechanical-equivalence tolerance.
    {
        std::vector<Vector> candidate_sv;
        r_candidate.pGetElement(1)->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, candidate_sv, r_candidate.GetProcessInfo());
        double min_abs_component = std::numeric_limits<double>::max();
        double max_abs_component = 0.0;
        for (std::size_t gp = 0; gp < candidate_sv.size(); ++gp) {
            KRATOS_EXPECT_EQ(candidate_sv[gp].size(), 6);
            for (std::size_t c = 0; c < 6; ++c) {
                min_abs_component = std::min(min_abs_component, std::abs(candidate_sv[gp](c)));
                max_abs_component = std::max(max_abs_component, std::abs(candidate_sv[gp](c)));
            }
        }
        // Non-uniform GP stress field (max-min of the largest component).
        double max_sxx = -std::numeric_limits<double>::max();
        double min_sxx = std::numeric_limits<double>::max();
        for (std::size_t gp = 0; gp < candidate_sv.size(); ++gp) {
            max_sxx = std::max(max_sxx, candidate_sv[gp](0));
            min_sxx = std::min(min_sxx, candidate_sv[gp](0));
        }
        KRATOS_EXPECT_TRUE((max_sxx - min_sxx) > 1.0e-3 * max_abs_component);
        KRATOS_EXPECT_TRUE(min_abs_component > 1.0e-4 * max_abs_component);
    }

    RunCandidateProcess(r_candidate);
    VerifyProcessAgainstAnalyticRaw(r_candidate, "Hexahedra3D8");
}

//************************************************************************************
// Shared-node unequal-area equivalence
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_SharedNode_UnequalAreas, KratosDamFastSuite)
{
    // Two connected quads sharing nodes 2 and 3, areas 1 and 2.
    // Quad 1: 1(0,0), 2(1,0), 3(1,1), 4(0,1)  - area 1
    // Quad 2: 2(1,0), 5(3,0), 6(3,1), 3(1,1)  - area 2
    const double coords[6][2] = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {3.0, 0.0}, {3.0, 1.0}};

    auto build_model = [&coords](Model& r_model, const std::string& r_element_name) -> ModelPart& {
        ModelPart& r_model_part = r_model.CreateModelPart("SharedNodes" + r_element_name, 2);
        ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        r_pi[DOMAIN_SIZE] = 2;
        r_pi[SPACE_DIMENSION] = 2;
        r_pi[IS_RESTARTED] = false;
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
        r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
        r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
        for (std::size_t i = 0; i < 6; ++i)
            r_model_part.CreateNewNode(i + 1, coords[i][0], coords[i][1], 0.0);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
            r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
            r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            Matrix zero_stress_tensor(2, 2);
            noalias(zero_stress_tensor) = ZeroMatrix(2, 2);
            r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
        }
        auto p_prop = r_model_part.CreateNewProperties(1);
        (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
        (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
        (*p_prop)[DENSITY] = test_density;
        (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
        (*p_prop)[THICKNESS] = test_thickness;
        p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic2DPlaneStrain().Clone());
        r_model_part.CreateNewElement(r_element_name, 1, {1, 2, 3, 4}, p_prop);
        r_model_part.CreateNewElement(r_element_name, 2, {2, 5, 6, 3}, p_prop);
        for (auto& r_node : r_model_part.Nodes()) {
            const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
            array_1d<double, 3>& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            r_disp[0] = 1.0e-3 * r_x0[0] + 3.0e-4 * r_x0[1] + 1.0e-4 * r_x0[0] * r_x0[1];
            r_disp[1] = 3.0e-4 * r_x0[0] - 1.0e-3 * r_x0[1] + 2.0e-4 * r_x0[0] * r_x0[0];
            r_node.X() = r_x0[0] + r_disp[0];
            r_node.Y() = r_x0[1] + r_disp[1];
            r_node.FastGetSolutionStepValue(TEMPERATURE) =
                test_reference_temperature + 3.0 * r_x0[0] + 2.0 * r_x0[1];
        }
        for (auto& r_element : r_model_part.Elements()) {
            KRATOS_EXPECT_EQ(r_element.Check(r_pi), 0);
            r_element.Initialize(r_pi);
            r_element.InitializeSolutionStep(r_pi);
            r_element.InitializeNonLinearIteration(r_pi);
        }
        return r_model_part;
    };

    Model model;
    ModelPart& r_candidate = build_model(model, "SmallDisplacementElement2D4N");
    RunCandidateProcess(r_candidate);
    VerifyProcessAgainstAnalyticRaw(r_candidate, "SharedNode2Quads");
}

//************************************************************************************
// Multi-step equivalence
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_MultiStep, KratosDamFastSuite)
{
    // Two consecutive steps with different thermo-mechanical states. After
    // Phase 3D.4B the process is the single production implementation, so each
    // step is verified against the independently computed analytic formula:
    //   1. reset NODAL_AREA / NODAL_CAUCHY_STRESS_TENSOR (as the smoothing scheme);
    //   2. accumulate through the process;
    //   3. verify the raw accumulators and NODAL_AREA against the analytic sum;
    //   4. normalize and verify the final smoothed nodal stress is non-trivial.
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "CandidateMS", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3);
    InitializeProcessElement(r_candidate);

    for (std::size_t step = 0; step < 2; ++step) {
        PrescribeVaryingState(r_candidate, 3);

        // Stage 1: reset as the smoothing scheme does.
        ResetNodalAccumulators(r_candidate, 3);

        // Stage 2: accumulate through the single production implementation.
        RunCandidateProcess(r_candidate);

        // Stage 3: verify against the analytic formula.
        VerifyProcessAgainstAnalyticRaw(r_candidate, "MultiStep step " + std::to_string(step + 1));
    }
}

//************************************************************************************
// Standard-element independence (only StructuralMechanics elements)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_StandardElementIndependence, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "StandardOnly", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3);
    PrescribeVaryingState(r_candidate, 3);
    InitializeProcessElement(r_candidate);
    RunCandidateProcess(r_candidate);

    // The process must have produced a non-trivial raw accumulated stress.
    bool non_trivial = false;
    for (auto& r_node : r_candidate.Nodes()) {
        if (r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(0, 0) != 0.0)
            non_trivial = true;
        KRATOS_EXPECT_TRUE(r_node.FastGetSolutionStepValue(NODAL_AREA) > 0.0);
    }
    KRATOS_EXPECT_TRUE(non_trivial);
}

//************************************************************************************
// No equilibrium side effects
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_NoEquilibriumSideEffects, KratosDamFastSuite)
{
    for (const char* p_element_name : {"SmallDisplacementElement2D4N",
                                       "SmallDisplacementElement3D8N"}) {
        const std::string r_element_name(p_element_name);
        const std::size_t dimension = (r_element_name.find("2D") != std::string::npos) ? 2 : 3;
        const std::string r_law_name =
            (dimension == 2) ? "ThermalLinearElastic2DPlaneStrain" : "ThermalLinearElastic3DLaw";

        Model model;
        ModelPart& r_model_part = CreateProcessModelPart(
            model, "SideEffect" + r_element_name, r_element_name, r_law_name, dimension);
        PrescribeVaryingState(r_model_part, dimension);
        InitializeProcessElement(r_model_part);
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

        Matrix lhs_before, lhs_after;
        Vector rhs_before, rhs_after;
        std::vector<Vector> cauchy_before, cauchy_after, pk2_before, pk2_after, strain_before, strain_after;
        p_element->CalculateLocalSystem(lhs_before, rhs_before, r_pi);
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_before, r_pi);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_before, r_pi);
        p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, strain_before, r_pi);

        RunCandidateProcess(r_model_part);

        p_element->CalculateLocalSystem(lhs_after, rhs_after, r_pi);
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_after, r_pi);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_after, r_pi);
        p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, strain_after, r_pi);

        for (std::size_t i = 0; i < lhs_before.size1(); ++i) {
            for (std::size_t j = 0; j < lhs_before.size2(); ++j) {
                KRATOS_EXPECT_NEAR(lhs_after(i, j), lhs_before(i, j), comparison_absolute_tolerance);
            }
        }
        for (std::size_t i = 0; i < rhs_before.size(); ++i) {
            KRATOS_EXPECT_NEAR(rhs_after(i), rhs_before(i), comparison_absolute_tolerance);
        }
        for (std::size_t gp = 0; gp < cauchy_before.size(); ++gp) {
            for (std::size_t c = 0; c < cauchy_before[gp].size(); ++c) {
                KRATOS_EXPECT_NEAR(cauchy_after[gp](c), cauchy_before[gp](c), comparison_absolute_tolerance);
                KRATOS_EXPECT_NEAR(pk2_after[gp](c), pk2_before[gp](c), comparison_absolute_tolerance);
                KRATOS_EXPECT_NEAR(strain_after[gp](c), strain_before[gp](c), comparison_absolute_tolerance);
            }
        }
    }
}

//************************************************************************************
// Prescribed-GP process-equivalence (identical GP inputs)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_PrescribedGP_Hexahedra3D8, KratosDamFastSuite)
{
    // Isolates the PROCESS extrapolation from the element strain computation:
    // both paths receive exactly the same prescribed GP stress field, so the
    // process must reproduce the legacy mathematical path at machine precision.
    Model model;
    PrescribedStressElement* p_element = nullptr;
    ModelPart& r_model_part = CreatePrescribedGPModelPart(
        model, "PrescribedH8", "SmallDisplacementSolidElement3D8N", 3, p_element);

    // Deliberately non-uniform 8-GP field with every stress component non-zero.
    std::vector<Matrix> prescribed(8);
    for (std::size_t gp = 0; gp < 8; ++gp) {
        Matrix t(3, 3);
        noalias(t) = ZeroMatrix(3, 3);
        t(0, 0) = 1000.0 + 100.0 * static_cast<double>(gp);
        t(1, 1) = 2000.0 + 200.0 * static_cast<double>(gp);
        t(2, 2) = 3000.0 + 300.0 * static_cast<double>(gp);
        t(0, 1) = t(1, 0) = 400.0 + 40.0 * static_cast<double>(gp);
        t(0, 2) = t(2, 0) = 500.0 + 50.0 * static_cast<double>(gp);
        t(1, 2) = t(2, 1) = 600.0 + 60.0 * static_cast<double>(gp);
        prescribed[gp] = t;
    }
    p_element->SetPrescribedStress(prescribed);

    DamNodalCauchyStressExtrapolationProcess process(r_model_part);
    process.ExecuteFinalizeSolutionStep();

    const auto& r_geometry = p_element->GetGeometry();
    const double expected_area = r_geometry.Area();
    for (auto& r_node : r_model_part.Nodes()) {
        const std::size_t node_index = r_node.Id() - 1;
        const Matrix expected_stress = ComputeExpectedNodalRawStress(r_geometry, prescribed, 3, node_index);
        const Matrix& r_actual = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                KRATOS_EXPECT_NEAR(r_actual(i, j), expected_stress(i, j),
                                   std::max(1.0e-12, 1.0e-10 * std::abs(expected_stress(i, j))));
            }
        }
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(NODAL_AREA), expected_area, 1.0e-12);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_PrescribedGP_Quadrilateral2D4, KratosDamFastSuite)
{
    Model model;
    PrescribedStressElement* p_element = nullptr;
    ModelPart& r_model_part = CreatePrescribedGPModelPart(
        model, "PrescribedQ4", "SmallDisplacementSolidElement2D4N", 2, p_element);

    std::vector<Matrix> prescribed(4);
    for (std::size_t gp = 0; gp < 4; ++gp) {
        Matrix t(2, 2);
        noalias(t) = ZeroMatrix(2, 2);
        t(0, 0) = 1000.0 + 100.0 * static_cast<double>(gp);
        t(1, 1) = 2000.0 + 200.0 * static_cast<double>(gp);
        t(0, 1) = t(1, 0) = 300.0 + 30.0 * static_cast<double>(gp);
        prescribed[gp] = t;
    }
    p_element->SetPrescribedStress(prescribed);

    DamNodalCauchyStressExtrapolationProcess process(r_model_part);
    process.ExecuteFinalizeSolutionStep();

    const auto& r_geometry = p_element->GetGeometry();
    const double expected_area = r_geometry.Area();
    for (auto& r_node : r_model_part.Nodes()) {
        const std::size_t node_index = r_node.Id() - 1;
        const Matrix expected_stress = ComputeExpectedNodalRawStress(r_geometry, prescribed, 2, node_index);
        const Matrix& r_actual = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                KRATOS_EXPECT_NEAR(r_actual(i, j), expected_stress(i, j),
                                   std::max(1.0e-12, 1.0e-10 * std::abs(expected_stress(i, j))));
            }
        }
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(NODAL_AREA), expected_area, 1.0e-12);
    }
}

//************************************************************************************
// Zero-shear roundoff characterization (documented, non-failing)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamNodalStressExtrapolation_ZeroShearRoundoffCharacterization, KratosDamFastSuite)
{
    // After Phase 3D.4B the process is the single production implementation and
    // there is no element-level reference. For a field with theoretically zero
    // yz/xz shears the process still reproduces the analytic legacy formula to
    // machine precision from the element's own integration-point stresses (the
    // process consumes the same GP stresses used to build the expected values).
    Model model;
    ModelPart& r_candidate = CreateProcessModelPart(
        model, "CandidateZS", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3);
    PrescribeZeroShearState(r_candidate, 3);
    InitializeProcessElement(r_candidate);

    ResetNodalAccumulators(r_candidate, 3);
    RunCandidateProcess(r_candidate);

    // The process result is verified against the analytic formula (scale-aware
    // tolerance); the zero-shear GP residues remain at the element arithmetic
    // level and do not produce pathological nodal values.
    double max_nodal_deviation = 0.0;
    double max_nodal_scale = 0.0;
    // Build the analytic expected map.
    std::map<ModelPart::IndexType, Matrix> expected_stress;
    std::map<ModelPart::IndexType, double> expected_area;
    for (auto& r_element : r_candidate.Elements()) {
        std::vector<Matrix> gauss_stress;
        r_element.CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, gauss_stress, r_candidate.GetProcessInfo());
        const auto& r_geom = r_element.GetGeometry();
        const std::size_t n_nodes = r_geom.size();
        const std::size_t n_gps = gauss_stress.size();
        Matrix container(n_gps, 6);
        noalias(container) = ZeroMatrix(n_gps, 6);
        for (std::size_t gp = 0; gp < n_gps; ++gp) {
            container(gp, 0) = gauss_stress[gp](0, 0);
            container(gp, 1) = gauss_stress[gp](1, 1);
            container(gp, 2) = gauss_stress[gp](2, 2);
            container(gp, 3) = gauss_stress[gp](0, 1);
            container(gp, 4) = gauss_stress[gp](1, 2);
            container(gp, 5) = gauss_stress[gp](0, 2);
        }
        Matrix aux;
        if (n_gps == 1) {
            aux = ZeroMatrix(n_nodes, 6);
            for (std::size_t i = 0; i < n_nodes; ++i)
                for (std::size_t c = 0; c < 6; ++c) aux(i, c) = container(0, c);
        } else {
            BoundedMatrix<double, 8, 8> h8;
            PoroElementUtilities::Calculate3DExtrapolationMatrix(h8);
            aux = prod(Matrix(h8), container);
        }
        const double area = r_geom.Area();
        for (std::size_t i = 0; i < n_nodes; ++i) {
            Vector nv(6);
            for (std::size_t c = 0; c < 6; ++c) nv(c) = aux(i, c);
            Matrix contrib = area * MathUtils<double>::StressVectorToTensor(nv);
            const ModelPart::IndexType nid = r_geom[i].Id();
            if (expected_stress.find(nid) == expected_stress.end()) {
                expected_stress[nid] = ZeroMatrix(3, 3);
                expected_area[nid] = 0.0;
            }
            expected_stress[nid] += contrib;
            expected_area[nid] += area;
        }
    }
    for (auto& r_node : r_candidate.Nodes()) {
        const ModelPart::IndexType nid = r_node.Id();
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(NODAL_AREA), expected_area[nid],
                           std::max(comparison_absolute_tolerance,
                                    comparison_relative_tolerance * std::abs(expected_area[nid])));
        const Matrix& stored = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                max_nodal_deviation = std::max(max_nodal_deviation, std::abs(stored(i, j) - expected_stress[nid](i, j)));
                max_nodal_scale = std::max(max_nodal_scale, std::abs(stored(i, j)));
            }
        }
    }
    std::cout << "[CHARACTERIZATION] zero-shear roundoff: max nodal raw deviation = "
              << max_nodal_deviation << " (scale " << max_nodal_scale << ")" << std::endl;
    KRATOS_EXPECT_TRUE(max_nodal_deviation < 1.0e-6 * std::max(1.0, max_nodal_scale));
}

} // namespace Testing
} // namespace Kratos
