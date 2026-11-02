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
#include <typeinfo>
#include <vector>

// Project includes
#include "includes/checks.h"
#include "includes/initial_state.h"
#include "includes/variables.h"
#include "custom_elements/solid_elements/small_displacement.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/elastic_isotropic_3d.h"

// Application includes
#include "custom_utilities/selective_initial_stress_utility.h"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{
namespace
{
constexpr double Tolerance = 1.e-12;

bool IsSupportedGeometry(const GeometryData::KratosGeometryType GeometryType)
{
    switch (GeometryType) {
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
        case GeometryData::KratosGeometryType::Kratos_Triangle2D6:
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8:
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9:
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10:
        case GeometryData::KratosGeometryType::Kratos_Prism3D6:
        case GeometryData::KratosGeometryType::Kratos_Prism3D15:
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27:
            return true;
        default:
            return false;
    }
}

bool IsSupportedIntegrationMethod(const GeometryData::IntegrationMethod IntegrationMethod)
{
    switch (IntegrationMethod) {
        case GeometryData::IntegrationMethod::GI_GAUSS_1:
        case GeometryData::IntegrationMethod::GI_GAUSS_2:
        case GeometryData::IntegrationMethod::GI_GAUSS_3:
        case GeometryData::IntegrationMethod::GI_GAUSS_4:
        case GeometryData::IntegrationMethod::GI_GAUSS_5:
            return true;
        default:
            return false;
    }
}

bool IsThermalLaw(const ConstitutiveLaw& rLaw)
{
    return typeid(rLaw) == typeid(ThermalLinearElastic2DPlaneStrain) ||
           typeid(rLaw) == typeid(ThermalLinearElastic2DPlaneStress) ||
           typeid(rLaw) == typeid(ThermalLinearElastic3DLaw);
}

bool IsSupportedLaw(const ConstitutiveLaw& rLaw)
{
    return typeid(rLaw) == typeid(LinearPlaneStrain) ||
           typeid(rLaw) == typeid(LinearPlaneStress) ||
           typeid(rLaw) == typeid(ElasticIsotropic3D) ||
           IsThermalLaw(rLaw);
}

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
            << "SelectiveInitialStressUtility: only SmallDisplacement solid elements are supported "
            << "for element " << element_id << "." << std::endl;

        auto& r_source_geometry = r_source_element.GetGeometry();
        auto& r_destination_geometry = r_destination_element.GetGeometry();
        const auto geometry_type = r_destination_geometry.GetGeometryType();
        KRATOS_ERROR_IF(r_source_geometry.GetGeometryType() != geometry_type)
            << "SelectiveInitialStressUtility: source and destination geometry types differ "
            << "for element " << element_id << "." << std::endl;
        KRATOS_ERROR_IF(!IsSupportedGeometry(geometry_type))
            << "SelectiveInitialStressUtility: unsupported geometry '"
            << r_destination_geometry.Name() << "' for element " << element_id
            << ". Supported geometries are T3, T6, Q4, Q8, Q9, T4, T10, P6, P15, H8, H20 and H27."
            << std::endl;

        const auto source_integration_method = r_source_element.GetIntegrationMethod();
        const auto destination_integration_method = r_destination_element.GetIntegrationMethod();
        KRATOS_ERROR_IF(source_integration_method != destination_integration_method)
            << "SelectiveInitialStressUtility: source and destination integration rules differ "
            << "for element " << element_id << "." << std::endl;
        KRATOS_ERROR_IF(!IsSupportedIntegrationMethod(destination_integration_method))
            << "SelectiveInitialStressUtility: unsupported integration method "
            << static_cast<int>(destination_integration_method) << " for element "
            << element_id << ". Only Gaussian integration rules GI_GAUSS_1 to GI_GAUSS_5 "
            << "are supported." << std::endl;
        const std::size_t num_points =
            r_destination_geometry.IntegrationPointsNumber(destination_integration_method);
        KRATOS_ERROR_IF(num_points == 0)
            << "SelectiveInitialStressUtility: the integration rule of element " << element_id
            << " has no integration points." << std::endl;

        const auto& r_source_points = r_source_geometry.IntegrationPoints(source_integration_method);
        const auto& r_destination_points =
            r_destination_geometry.IntegrationPoints(destination_integration_method);
        KRATOS_ERROR_IF(r_source_points.size() != num_points || r_destination_points.size() != num_points)
            << "SelectiveInitialStressUtility: integration-point count mismatch for element "
            << element_id << "." << std::endl;
        for (std::size_t i = 0; i < num_points; ++i) {
            KRATOS_ERROR_IF(norm_2(r_source_points[i].Coordinates() -
                                  r_destination_points[i].Coordinates()) > Tolerance)
                << "SelectiveInitialStressUtility: integration-point ordering differs at point "
                << i << " of element " << element_id << "." << std::endl;
        }

        for (std::size_t i = 0; i < r_destination_geometry.size(); ++i) {
            const auto& r_source_node = r_source_geometry[i];
            const auto& r_destination_node = r_destination_geometry[i];
            KRATOS_ERROR_IF(r_source_node.Id() != r_destination_node.Id())
                << "SelectiveInitialStressUtility: node ordering mismatch at position " << i
                << " of element " << element_id << " (" << r_source_node.Id()
                << " vs " << r_destination_node.Id() << ")." << std::endl;
            KRATOS_ERROR_IF(norm_2(r_destination_node.GetInitialPosition() -
                                  r_source_node.GetInitialPosition()) > Tolerance)
                << "SelectiveInitialStressUtility: reference coordinates differ for node "
                << r_destination_node.Id() << "." << std::endl;
            KRATOS_ERROR_IF(norm_2(r_destination_node.Coordinates() -
                                  r_destination_node.GetInitialPosition()) > Tolerance)
                << "SelectiveInitialStressUtility: destination node " << r_destination_node.Id()
                << " has moved from its initial position." << std::endl;
            CheckZeroDisplacement(r_destination_node);
        }

        auto& r_source_properties = r_source_element.GetProperties();
        auto& r_destination_properties = r_destination_element.GetProperties();
        for (const auto* p_variable : {&YOUNG_MODULUS, &POISSON_RATIO, &DENSITY}) {
            CheckMaterialProperty(r_source_properties, r_destination_properties, *p_variable);
        }
        if (r_destination_geometry.WorkingSpaceDimension() == 2) {
            CheckMaterialProperty(r_source_properties, r_destination_properties, THICKNESS);
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

        std::size_t strain_size = 0;
        bool has_thermal_law = false;
        bool has_nonthermal_law = false;
        for (std::size_t i = 0; i < num_points; ++i) {
            KRATOS_ERROR_IF(!source_laws[i] || !destination_laws[i])
                << "SelectiveInitialStressUtility: null constitutive law at element "
                << element_id << ", integration point " << i << "." << std::endl;
            const auto& r_source_law = *source_laws[i];
            auto& r_destination_law = *destination_laws[i];
            KRATOS_ERROR_IF(!IsSupportedLaw(r_destination_law))
                << "SelectiveInitialStressUtility: unsupported constitutive law '"
                << r_destination_law.Info() << "' at element " << element_id
                << ", integration point " << i << ". Supported laws are "
                << "LinearElasticPlaneStrain2DLaw, LinearElasticPlaneStress2DLaw, "
                << "LinearElastic3DLaw, ThermalLinearElastic2DPlaneStrain, "
                << "ThermalLinearElastic2DPlaneStress and ThermalLinearElastic3DLaw."
                << std::endl;
            KRATOS_ERROR_IF(typeid(r_source_law) != typeid(r_destination_law))
                << "SelectiveInitialStressUtility: source and destination constitutive laws "
                << "differ at element " << element_id << ", integration point " << i << "."
                << std::endl;
            const std::size_t law_strain_size = r_destination_law.GetStrainSize();
            if (i == 0) {
                strain_size = law_strain_size;
            } else {
                KRATOS_ERROR_IF(law_strain_size != strain_size)
                    << "SelectiveInitialStressUtility: inconsistent strain sizes within element "
                    << element_id << "." << std::endl;
            }
            KRATOS_ERROR_IF(r_destination_law.WorkingSpaceDimension() !=
                            r_destination_geometry.WorkingSpaceDimension())
                << "SelectiveInitialStressUtility: constitutive law and geometry dimensions "
                << "differ at element " << element_id << ", integration point " << i << "."
                << std::endl;
            KRATOS_ERROR_IF(r_destination_law.HasInitialState())
                << "SelectiveInitialStressUtility: destination element " << element_id
                << ", integration point " << i
                << " already has an initial state. The initial state is applied only once."
                << std::endl;
            const bool is_thermal = IsThermalLaw(r_destination_law);
            has_thermal_law = has_thermal_law || is_thermal;
            has_nonthermal_law = has_nonthermal_law || !is_thermal;
        }

        KRATOS_ERROR_IF(has_thermal_law && has_nonthermal_law)
            << "SelectiveInitialStressUtility: mixed thermal and non-thermal constitutive laws "
            << "are not supported for element " << element_id << "." << std::endl;
        const std::size_t expected_strain_size =
            r_destination_geometry.WorkingSpaceDimension() == 2 ? 3 : 6;
        KRATOS_ERROR_IF(strain_size != expected_strain_size)
            << "SelectiveInitialStressUtility: unsupported strain size " << strain_size
            << " for element " << element_id << "; expected " << expected_strain_size
            << " for a " << r_destination_geometry.WorkingSpaceDimension()
            << "D geometry." << std::endl;

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
            KRATOS_ERROR_IF(stresses[i].size() != strain_size)
                << "SelectiveInitialStressUtility: unsupported Cauchy stress vector size "
                << stresses[i].size() << " at element " << element_id
                << ", integration point " << i << ". Expected " << strain_size
                << " components." << std::endl;
            for (const double component : stresses[i]) {
                KRATOS_ERROR_IF(!std::isfinite(component))
                    << "SelectiveInitialStressUtility: non-finite stress component at element "
                    << element_id << ", integration point " << i << "." << std::endl;
            }
        }

        point_data.push_back({std::move(stresses), std::move(destination_laws)});
    }

    // All checks passed. Assign the selected stresses as initial states.
    for (const auto& r_data : point_data) {
        for (std::size_t i = 0; i < r_data.Laws.size(); ++i) {
            auto& r_law = *r_data.Laws[i];
            r_law.SetInitialState(Kratos::make_intrusive<InitialState>(
                ZeroVector(r_data.Stresses[i].size()), r_data.Stresses[i],
                IdentityMatrix(r_law.WorkingSpaceDimension())));
        }
    }

    KRATOS_CATCH("SelectiveInitialStressUtility::AssignInitialStresses")
}

} // namespace Kratos
