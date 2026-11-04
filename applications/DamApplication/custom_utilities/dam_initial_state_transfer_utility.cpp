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
#include <map>
#include <string>
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
#include "custom_utilities/dam_initial_state_transfer_utility.h"
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
        << "DamInitialStateTransferUtility: node " << rNode.Id()
        << " does not have DISPLACEMENT in its solution step data." << std::endl;

    for (std::size_t step = 0; step < rNode.GetBufferSize(); ++step) {
        KRATOS_ERROR_IF(norm_2(rNode.FastGetSolutionStepValue(DISPLACEMENT, step)) > Tolerance)
            << "DamInitialStateTransferUtility: destination node " << rNode.Id()
            << " has a nonzero displacement in buffer step " << step
            << ". Selective initial stresses require the original zero displacement datum."
            << std::endl;
    }
}

void CheckThermalDatum(const Node& rNode)
{
    const double temperature = rNode.FastGetSolutionStepValue(TEMPERATURE);
    const double reference_temperature = rNode.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE);
    const double scale = 1. + std::abs(reference_temperature);
    KRATOS_ERROR_IF(std::abs(temperature - reference_temperature) > Tolerance * scale)
        << "DamInitialStateTransferUtility: thermoelastic initial-stress import requires "
        << "T == T_ref at the moment of import, but node " << rNode.Id()
        << " has T = " << temperature << " and T_ref = " << reference_temperature
        << ". Prescribe a matching reference temperature before importing stresses."
        << std::endl;
}

struct PointData
{
    std::vector<Vector> Stresses;
    std::vector<Vector> InitialStrains;
    std::vector<ConstitutiveLaw::Pointer> Laws;
};

void AssignPointData(const std::vector<PointData>& rPointData)
{
    for (const auto& r_data : rPointData) {
        for (std::size_t i = 0; i < r_data.Laws.size(); ++i) {
            auto& r_law = *r_data.Laws[i];
            Vector initial_strain(r_data.Stresses[i].size(), 0.0);
            if (r_data.InitialStrains.size() == r_data.Stresses.size()) {
                initial_strain = r_data.InitialStrains[i];
            }
            r_law.SetInitialState(Kratos::make_intrusive<InitialState>(
                initial_strain, r_data.Stresses[i],
                IdentityMatrix(r_law.WorkingSpaceDimension())));
        }
    }
}

Vector ComputeThermalStrain(
    const Element& rElement,
    ConstitutiveLaw& rLaw,
    const ProcessInfo& rProcessInfo,
    const std::size_t IntegrationPoint)
{
    const auto& r_geometry = rElement.GetGeometry();
    ConstitutiveLaw::Parameters values(
        r_geometry, rElement.GetProperties(), rProcessInfo);

    const Matrix& r_shape_functions =
        r_geometry.ShapeFunctionsValues(rElement.GetIntegrationMethod());
    Vector shape_functions(r_geometry.size(), 0.0);
    for (std::size_t i = 0; i < r_geometry.size(); ++i) {
        shape_functions[i] = r_shape_functions(IntegrationPoint, i);
    }
    values.SetShapeFunctionsValues(shape_functions);

    Vector strain(rLaw.GetStrainSize(), 0.0);
    values.SetStrainVector(strain);

    Vector thermal_strain;
    rLaw.CalculateValue(values, THERMAL_STRAIN_VECTOR, thermal_strain);
    KRATOS_ERROR_IF(thermal_strain.size() != rLaw.GetStrainSize())
        << "DamInitialStateTransferUtility: thermal-strain size "
        << thermal_strain.size() << " does not match strain size "
        << rLaw.GetStrainSize() << " at element " << rElement.Id()
        << ", integration point " << IntegrationPoint << "." << std::endl;
    return thermal_strain;
}
} // namespace

// -----------------------------------------------------------------------------
// Lightweight initial-state file support
// -----------------------------------------------------------------------------

std::map<IndexType, std::vector<Vector>> DamInitialStateTransferUtility::GetFinalCauchyStresses(
    ModelPart& rModelPart)
{
    KRATOS_TRY
    std::map<IndexType, std::vector<Vector>> result;

    for (auto& r_element : rModelPart.Elements()) {
        KRATOS_ERROR_IF(typeid(r_element) != typeid(SmallDisplacement))
            << "DamInitialStateTransferUtility: only SmallDisplacement solid elements are supported "
            << "for element " << r_element.Id() << "." << std::endl;
        auto& r_geometry = r_element.GetGeometry();
        KRATOS_ERROR_IF(!IsSupportedGeometry(r_geometry.GetGeometryType()))
            << "DamInitialStateTransferUtility: unsupported geometry for element "
            << r_element.Id() << "." << std::endl;

        std::vector<ConstitutiveLaw::Pointer> laws;
        std::vector<Vector> stresses;
        r_element.CalculateOnIntegrationPoints(
            CONSTITUTIVE_LAW, laws, rModelPart.GetProcessInfo());
        r_element.CalculateOnIntegrationPoints(
            CAUCHY_STRESS_VECTOR, stresses, rModelPart.GetProcessInfo());
        KRATOS_ERROR_IF(laws.size() != stresses.size())
            << "DamInitialStateTransferUtility: constitutive-law and stress integration-point "
            << "counts differ for element " << r_element.Id() << "." << std::endl;
        for (std::size_t i = 0; i < stresses.size(); ++i) {
            KRATOS_ERROR_IF(!laws[i] || !IsSupportedLaw(*laws[i]))
                << "DamInitialStateTransferUtility: unsupported constitutive law at element "
                << r_element.Id() << ", integration point " << i << "." << std::endl;
            KRATOS_ERROR_IF(stresses[i].size() != laws[i]->GetStrainSize())
                << "DamInitialStateTransferUtility: stress-vector size "
                << stresses[i].size() << " does not match strain size "
                << laws[i]->GetStrainSize() << " at element " << r_element.Id()
                << ", integration point " << i << "." << std::endl;
            for (const double component : stresses[i]) {
                KRATOS_ERROR_IF(!std::isfinite(component))
                    << "DamInitialStateTransferUtility: non-finite Cauchy stress at element "
                    << r_element.Id() << ", integration point " << i << "." << std::endl;
            }
        }
        result.emplace(r_element.Id(), std::move(stresses));
    }

    return result;
    KRATOS_CATCH("DamInitialStateTransferUtility::GetFinalCauchyStresses")
}

namespace
{
std::vector<PointData> BuildStoredPointData(
    ModelPart& rDestinationModelPart,
    const std::map<IndexType, std::vector<Vector>>& rStresses,
    const bool CompensateInitialThermalStrain)
{
    KRATOS_ERROR_IF(rStresses.size() != rDestinationModelPart.NumberOfElements())
        << "DamInitialStateTransferUtility: imported stress map contains "
        << rStresses.size() << " elements but the destination contains "
        << rDestinationModelPart.NumberOfElements() << "." << std::endl;

    std::vector<PointData> point_data;
    point_data.reserve(rDestinationModelPart.NumberOfElements());
    for (auto& r_element : rDestinationModelPart.Elements()) {
        const IndexType element_id = r_element.Id();
        const auto it_stress = rStresses.find(element_id);
        KRATOS_ERROR_IF(it_stress == rStresses.end())
            << "DamInitialStateTransferUtility: imported stresses do not contain element "
            << element_id << "." << std::endl;
        KRATOS_ERROR_IF(typeid(r_element) != typeid(SmallDisplacement))
            << "DamInitialStateTransferUtility: only SmallDisplacement solid elements are supported "
            << "for element " << element_id << "." << std::endl;
        auto& r_geometry = r_element.GetGeometry();
        KRATOS_ERROR_IF(!IsSupportedGeometry(r_geometry.GetGeometryType()))
            << "DamInitialStateTransferUtility: unsupported geometry for element "
            << element_id << "." << std::endl;

        for (const auto& r_node : r_geometry) {
            CheckZeroDisplacement(r_node);
        }

        std::vector<ConstitutiveLaw::Pointer> laws;
        r_element.CalculateOnIntegrationPoints(
            CONSTITUTIVE_LAW, laws, rDestinationModelPart.GetProcessInfo());
        const auto& r_stresses = it_stress->second;
        KRATOS_ERROR_IF(laws.size() != r_stresses.size())
            << "DamInitialStateTransferUtility: imported stress count "
            << r_stresses.size() << " does not match the integration-point count "
            << laws.size() << " of element " << element_id << "." << std::endl;

        bool has_thermal_law = false;
        bool has_nonthermal_law = false;
        for (std::size_t i = 0; i < laws.size(); ++i) {
            KRATOS_ERROR_IF(!laws[i] || !IsSupportedLaw(*laws[i]))
                << "DamInitialStateTransferUtility: unsupported constitutive law at element "
                << element_id << ", integration point " << i << "." << std::endl;
            KRATOS_ERROR_IF(r_stresses[i].size() != laws[i]->GetStrainSize())
                << "DamInitialStateTransferUtility: imported stress size "
                << r_stresses[i].size() << " does not match strain size "
                << laws[i]->GetStrainSize() << " at element " << element_id
                << ", integration point " << i << "." << std::endl;
            KRATOS_ERROR_IF(laws[i]->HasInitialState())
                << "DamInitialStateTransferUtility: destination element " << element_id
                << ", integration point " << i
                << " already has an initial state. The initial state is applied only once."
                << std::endl;
            for (const double component : r_stresses[i]) {
                KRATOS_ERROR_IF(!std::isfinite(component))
                    << "DamInitialStateTransferUtility: non-finite imported stress at element "
                    << element_id << ", integration point " << i << "." << std::endl;
            }
            const bool is_thermal = IsThermalLaw(*laws[i]);
            has_thermal_law = has_thermal_law || is_thermal;
            has_nonthermal_law = has_nonthermal_law || !is_thermal;
        }
        KRATOS_ERROR_IF(has_thermal_law && has_nonthermal_law)
            << "DamInitialStateTransferUtility: mixed thermal and non-thermal constitutive laws "
            << "are not supported for element " << element_id << "." << std::endl;

        std::vector<Vector> initial_strains;
        if (has_thermal_law && CompensateInitialThermalStrain) {
            initial_strains.reserve(laws.size());
            for (std::size_t i = 0; i < laws.size(); ++i) {
                Vector initial_strain = ComputeThermalStrain(
                    r_element, *laws[i], rDestinationModelPart.GetProcessInfo(), i);
                initial_strain *= -1.0;
                initial_strains.push_back(std::move(initial_strain));
            }
        }
        if (has_thermal_law && !CompensateInitialThermalStrain) {
            for (const auto& r_node : r_geometry) {
                CheckThermalDatum(r_node);
            }
        }
        point_data.push_back(
            {r_stresses, std::move(initial_strains), std::move(laws)});
    }
    return point_data;
}
} // namespace

void DamInitialStateTransferUtility::ValidateInitialStresses(
    ModelPart& rDestinationModelPart,
    const std::map<IndexType, std::vector<Vector>>& rStresses,
    const bool CompensateInitialThermalStrain)
{
    KRATOS_TRY
    BuildStoredPointData(
        rDestinationModelPart, rStresses, CompensateInitialThermalStrain);
    KRATOS_CATCH("DamInitialStateTransferUtility::ValidateInitialStresses")
}

void DamInitialStateTransferUtility::AssignInitialStresses(
    ModelPart& rDestinationModelPart,
    const std::map<IndexType, std::vector<Vector>>& rStresses,
    const bool CompensateInitialThermalStrain)
{
    KRATOS_TRY
    AssignPointData(BuildStoredPointData(
        rDestinationModelPart, rStresses, CompensateInitialThermalStrain));
    KRATOS_CATCH("DamInitialStateTransferUtility::AssignInitialStresses")
}

std::string DamInitialStateTransferUtility::GetIntegrationMethodName(const Element& rElement)
{
    switch (rElement.GetIntegrationMethod()) {
        case GeometryData::IntegrationMethod::GI_GAUSS_1: return "GI_GAUSS_1";
        case GeometryData::IntegrationMethod::GI_GAUSS_2: return "GI_GAUSS_2";
        case GeometryData::IntegrationMethod::GI_GAUSS_3: return "GI_GAUSS_3";
        case GeometryData::IntegrationMethod::GI_GAUSS_4: return "GI_GAUSS_4";
        case GeometryData::IntegrationMethod::GI_GAUSS_5: return "GI_GAUSS_5";
        default: return "UNSUPPORTED";
    }
}

std::vector<Vector> DamInitialStateTransferUtility::GetIntegrationPointCoordinates(
    const Element& rElement)
{
    const auto& r_points = rElement.GetGeometry().IntegrationPoints(
        rElement.GetIntegrationMethod());
    std::vector<Vector> result;
    result.reserve(r_points.size());
    for (const auto& r_point : r_points) {
        Vector coordinates(3);
        for (std::size_t i = 0; i < 3; ++i) {
            coordinates[i] = r_point.Coordinates()[i];
        }
        result.push_back(std::move(coordinates));
    }
    return result;
}

} // namespace Kratos
