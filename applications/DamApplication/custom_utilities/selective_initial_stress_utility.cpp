//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication contributors
//

// System includes
#include <cmath>
#include <vector>

// Project includes
#include "includes/checks.h"
#include "includes/initial_state.h"
#include "includes/variables.h"
#include "custom_elements/solid_elements/small_displacement.h"
#include "custom_constitutive/linear_plane_strain.h"

// Application includes
#include "custom_utilities/selective_initial_stress_utility.h"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "dam_application_variables.h"

namespace Kratos
{
namespace
{
constexpr double Tolerance = 1.e-12;
constexpr std::size_t StressSize = 3;

void CheckZeroDisplacement(const Node& rNode)
{
    KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(DISPLACEMENT))
        << "SelectiveInitialStressUtility: node " << rNode.Id()
        << " does not have DISPLACEMENT in its solution step data." << std::endl;

    for (std::size_t step = 0; step < rNode.GetBufferSize(); ++step) {
        KRATOS_ERROR_IF(norm_2(rNode.FastGetSolutionStepValue(DISPLACEMENT, step)) > Tolerance)
            << "SelectiveInitialStressUtility: destination node " << rNode.Id()
            << " has a nonzero displacement in buffer step " << step
            << ". Selective initial stresses require the original zero displacement datum."
            << std::endl;
    }
}

void CheckMaterialProperty(
    const Properties& rSourceProperties,
    const Properties& rDestinationProperties,
    const Variable<double>& rVariable)
{
    KRATOS_ERROR_IF_NOT(rSourceProperties.Has(rVariable))
        << "SelectiveInitialStressUtility: source properties do not define "
        << rVariable.Name() << "." << std::endl;
    KRATOS_ERROR_IF_NOT(rDestinationProperties.Has(rVariable))
        << "SelectiveInitialStressUtility: destination properties do not define "
        << rVariable.Name() << "." << std::endl;
    KRATOS_ERROR_IF(rSourceProperties[rVariable] != rDestinationProperties[rVariable])
        << "SelectiveInitialStressUtility: source and destination material mismatch for "
        << rVariable.Name() << " (" << rSourceProperties[rVariable]
        << " vs " << rDestinationProperties[rVariable] << ")." << std::endl;
}

void CheckThermalDatum(const Node& rNode)
{
    const double temperature = rNode.FastGetSolutionStepValue(TEMPERATURE);
    const double reference_temperature = rNode.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE);
    const double scale = 1. + std::abs(reference_temperature);
    KRATOS_ERROR_IF(std::abs(temperature - reference_temperature) > Tolerance * scale)
        << "SelectiveInitialStressUtility: thermoelastic initial-stress import requires "
        << "T == T_ref at the moment of import, but node " << rNode.Id()
        << " has T = " << temperature << " and T_ref = " << reference_temperature
        << ". Prescribe a matching reference temperature before importing stresses."
        << std::endl;
}

struct PointData
{
    IndexType ElementId;
    std::vector<Vector> Stresses;
    std::vector<ConstitutiveLaw::Pointer> Laws;
};
} // namespace

void SelectiveInitialStressUtility::AssignInitialStresses(
    ModelPart& rSourceModelPart,
    ModelPart& rDestinationModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(&rSourceModelPart == &rDestinationModelPart)
        << "SelectiveInitialStressUtility: source and destination model parts must be independent."
        << std::endl;
    KRATOS_ERROR_IF(rSourceModelPart.NumberOfNodes() == 0 || rDestinationModelPart.NumberOfNodes() == 0)
        << "SelectiveInitialStressUtility: source and destination model parts must contain nodes."
        << std::endl;
    KRATOS_ERROR_IF(rSourceModelPart.NumberOfElements() != rDestinationModelPart.NumberOfElements())
        << "SelectiveInitialStressUtility: element count mismatch (" << rSourceModelPart.NumberOfElements()
        << " vs " << rDestinationModelPart.NumberOfElements() << ")." << std::endl;
    KRATOS_ERROR_IF(rSourceModelPart.NumberOfNodes() != rDestinationModelPart.NumberOfNodes())
        << "SelectiveInitialStressUtility: node count mismatch (" << rSourceModelPart.NumberOfNodes()
        << " vs " << rDestinationModelPart.NumberOfNodes() << ")." << std::endl;

    std::vector<PointData> point_data;
    point_data.reserve(rDestinationModelPart.NumberOfElements());

    // Validate every element and integration point before touching any law.
    for (auto& r_destination_element : rDestinationModelPart.Elements()) {
        const IndexType element_id = r_destination_element.Id();
        KRATOS_ERROR_IF_NOT(rSourceModelPart.HasElement(element_id))
            << "SelectiveInitialStressUtility: source model part does not contain element "
            << element_id << "." << std::endl;
        auto& r_source_element = rSourceModelPart.GetElement(element_id);

        KRATOS_ERROR_IF(typeid(r_source_element) != typeid(SmallDisplacement) ||
                        typeid(r_destination_element) != typeid(SmallDisplacement))
            << "SelectiveInitialStressUtility: only SmallDisplacement elements are supported "
            << "for element " << element_id << "." << std::endl;

        auto& r_source_geometry = r_source_element.GetGeometry();
        auto& r_destination_geometry = r_destination_element.GetGeometry();
        KRATOS_ERROR_IF(r_destination_geometry.GetGeometryType() !=
                        GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 ||
                        r_source_geometry.GetGeometryType() !=
                        GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
            << "SelectiveInitialStressUtility: only quadrilateral 2D4 geometries are supported "
            << "for element " << element_id << "." << std::endl;

        const auto source_integration_method = r_source_element.GetIntegrationMethod();
        const auto destination_integration_method = r_destination_element.GetIntegrationMethod();
        KRATOS_ERROR_IF(source_integration_method != destination_integration_method)
            << "SelectiveInitialStressUtility: source and destination integration rules differ "
            << "for element " << element_id << "." << std::endl;
        KRATOS_ERROR_IF(destination_integration_method != GeometryData::IntegrationMethod::GI_GAUSS_2)
            << "SelectiveInitialStressUtility: only the GI_GAUSS_2 integration rule is supported "
            << "for element " << element_id << "." << std::endl;
        const std::size_t num_points =
            r_destination_geometry.IntegrationPointsNumber(destination_integration_method);
        KRATOS_ERROR_IF(num_points != 4)
            << "SelectiveInitialStressUtility: expected four integration points for element "
            << element_id << ", found " << num_points << "." << std::endl;

        KRATOS_ERROR_IF(r_source_geometry.size() != 4 || r_destination_geometry.size() != 4)
            << "SelectiveInitialStressUtility: expected four nodes per element "
            << element_id << "." << std::endl;
        for (std::size_t i = 0; i < 4; ++i) {
            const auto& r_source_node = r_source_geometry[i];
            const auto& r_destination_node = r_destination_geometry[i];
            KRATOS_ERROR_IF(r_source_node.Id() != r_destination_node.Id())
                << "SelectiveInitialStressUtility: node ordering mismatch at position " << i
                << " of element " << element_id << " (" << r_source_node.Id()
                << " vs " << r_destination_node.Id() << ")." << std::endl;
            KRATOS_ERROR_IF(norm_2(r_destination_node.GetInitialPosition() -
                                  r_source_node.GetInitialPosition()) > Tolerance)
                << "SelectiveInitialStressUtility: initial coordinates differ for node "
                << r_destination_node.Id() << "." << std::endl;
            KRATOS_ERROR_IF(norm_2(r_destination_node.Coordinates() -
                                  r_destination_node.GetInitialPosition()) > Tolerance)
                << "SelectiveInitialStressUtility: destination node " << r_destination_node.Id()
                << " has moved from its initial position." << std::endl;
            KRATOS_ERROR_IF(norm_2(r_source_node.Coordinates() -
                                  r_source_node.GetInitialPosition()) > Tolerance)
                << "SelectiveInitialStressUtility: source node " << r_source_node.Id()
                << " has moved from its initial position. Only fixed-mesh sources are supported."
                << std::endl;
            CheckZeroDisplacement(r_destination_node);
        }

        auto& r_source_properties = r_source_element.GetProperties();
        auto& r_destination_properties = r_destination_element.GetProperties();
        for (const auto* p_variable : {&YOUNG_MODULUS, &POISSON_RATIO, &DENSITY, &THICKNESS}) {
            CheckMaterialProperty(r_source_properties, r_destination_properties, *p_variable);
        }

        std::vector<ConstitutiveLaw::Pointer> source_laws;
        std::vector<ConstitutiveLaw::Pointer> destination_laws;
        r_source_element.CalculateOnIntegrationPoints(
            CONSTITUTIVE_LAW, source_laws, rSourceModelPart.GetProcessInfo());
        r_destination_element.CalculateOnIntegrationPoints(
            CONSTITUTIVE_LAW, destination_laws, rDestinationModelPart.GetProcessInfo());
        KRATOS_ERROR_IF(source_laws.size() != num_points || destination_laws.size() != num_points)
            << "SelectiveInitialStressUtility: uninitialized constitutive laws or integration-point "
            << "mismatch for element " << element_id << "." << std::endl;

        bool has_thermal_law = false;
        bool has_linear_law = false;
        for (std::size_t i = 0; i < num_points; ++i) {
            KRATOS_ERROR_IF(!source_laws[i] || !destination_laws[i])
                << "SelectiveInitialStressUtility: null constitutive law at element "
                << element_id << ", integration point " << i << "." << std::endl;
            const auto& r_source_law = *source_laws[i];
            const auto& r_destination_law = *destination_laws[i];
            const bool is_thermal =
                typeid(r_destination_law) == typeid(ThermalLinearElastic2DPlaneStrain);
            has_thermal_law = has_thermal_law || is_thermal;
            has_linear_law = has_linear_law || typeid(r_destination_law) == typeid(LinearPlaneStrain);
            KRATOS_ERROR_IF(typeid(r_destination_law) != typeid(LinearPlaneStrain) && !is_thermal)
                << "SelectiveInitialStressUtility: unsupported constitutive law at element "
                << element_id << ", integration point " << i << ". Supported laws are "
                << "LinearElasticPlaneStrain2DLaw and ThermalLinearElastic2DPlaneStrain."
                << std::endl;
            KRATOS_ERROR_IF(typeid(r_source_law) != typeid(r_destination_law))
                << "SelectiveInitialStressUtility: source and destination constitutive laws "
                << "differ at element " << element_id << ", integration point " << i << "."
                << std::endl;
            KRATOS_ERROR_IF(r_destination_law.GetStrainSize() != StressSize)
                << "SelectiveInitialStressUtility: unsupported strain size "
                << r_destination_law.GetStrainSize() << " at element " << element_id
                << ", integration point " << i << "." << std::endl;
            KRATOS_ERROR_IF(r_destination_law.HasInitialState())
                << "SelectiveInitialStressUtility: destination element " << element_id
                << ", integration point " << i
                << " already has an initial state. The initial state is applied only once."
                << std::endl;
        }

        KRATOS_ERROR_IF(has_thermal_law && has_linear_law)
            << "SelectiveInitialStressUtility: mixed thermal and non-thermal constitutive laws "
            << "are not supported for element " << element_id << "." << std::endl;

        if (has_thermal_law) {
            CheckMaterialProperty(r_source_properties, r_destination_properties, THERMAL_EXPANSION);
            for (const auto& r_node : r_destination_geometry) {
                CheckThermalDatum(r_node);
            }
        }

        std::vector<Vector> stresses;
        r_source_element.CalculateOnIntegrationPoints(
            CAUCHY_STRESS_VECTOR, stresses, rSourceModelPart.GetProcessInfo());
        KRATOS_ERROR_IF(stresses.size() != num_points)
            << "SelectiveInitialStressUtility: expected " << num_points
            << " stress vectors for element " << element_id << ", found "
            << stresses.size() << "." << std::endl;
        for (std::size_t i = 0; i < num_points; ++i) {
            KRATOS_ERROR_IF(stresses[i].size() != StressSize)
                << "SelectiveInitialStressUtility: unsupported Cauchy stress vector size "
                << stresses[i].size() << " at element " << element_id
                << ", integration point " << i << ". Expected " << StressSize << " components."
                << std::endl;
            for (const double component : stresses[i]) {
                KRATOS_ERROR_IF(!std::isfinite(component))
                    << "SelectiveInitialStressUtility: non-finite stress component at element "
                    << element_id << ", integration point " << i << "." << std::endl;
            }
        }

        point_data.push_back({element_id, std::move(stresses), std::move(destination_laws)});
    }

    // All checks passed. Assign the selected stresses as initial states.
    for (const auto& r_data : point_data) {
        for (std::size_t i = 0; i < r_data.Laws.size(); ++i) {
            r_data.Laws[i]->SetInitialState(Kratos::make_intrusive<InitialState>(
                ZeroVector(StressSize), r_data.Stresses[i], IdentityMatrix(2)));
        }
    }

    KRATOS_CATCH("SelectiveInitialStressUtility::AssignInitialStresses")
}

} // namespace Kratos
