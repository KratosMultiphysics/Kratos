// Test-only bridge: existing C++ interfaces missing from the Python API.
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "includes/initial_state.h"
#include "custom_elements/solid_elements/small_displacement.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "dam_application_variables.h"

namespace Kratos {
namespace {
void assign_stresses(ModelPart& rSource, ModelPart& rTarget)
{
    KRATOS_ERROR_IF(&rSource == &rTarget) << "Source and destination must be independent";
    KRATOS_ERROR_IF(rSource.NumberOfElements() != rTarget.NumberOfElements() ||
                    rSource.NumberOfNodes() != rTarget.NumberOfNodes()) << "Mesh mismatch";
    // Validate everything before assigning any point. Only the two reviewed laws.
    for (auto& r_target : rTarget.Elements()) {
        KRATOS_ERROR_IF_NOT(rSource.HasElement(r_target.Id())) << "Element ID mismatch";
        auto& r_source = rSource.GetElement(r_target.Id());
        KRATOS_ERROR_IF(typeid(r_source) != typeid(SmallDisplacement) ||
                        typeid(r_target) != typeid(SmallDisplacement)) << "Unsupported element";
        const auto& r_geometry = r_target.GetGeometry();
        KRATOS_ERROR_IF(r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 ||
            r_geometry.GetGeometryType() != r_source.GetGeometry().GetGeometryType()) << "Unsupported geometry";
        KRATOS_ERROR_IF(r_source.GetIntegrationMethod() != r_target.GetIntegrationMethod() ||
            r_target.GetIntegrationMethod() != GeometryData::IntegrationMethod::GI_GAUSS_2) << "Integration rule mismatch";
        for (std::size_t i = 0; i < r_geometry.size(); ++i) {
            const auto& r_node = r_geometry[i];
            const auto& r_old = r_source.GetGeometry()[i];
            KRATOS_ERROR_IF(r_node.Id() != r_old.Id() ||
                norm_2(r_node.GetInitialPosition() - r_old.GetInitialPosition()) > 1.e-12 ||
                norm_2(r_node.Coordinates() - r_node.GetInitialPosition()) > 1.e-12 ||
                norm_2(r_old.Coordinates() - r_old.GetInitialPosition()) > 1.e-12) << "Mesh mismatch";
            KRATOS_ERROR_IF(norm_2(r_node.FastGetSolutionStepValue(DISPLACEMENT)) != 0.) << "Destination must have zero displacement";
        }
        std::vector<ConstitutiveLaw::Pointer> source_laws, target_laws;
        r_source.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, source_laws, rSource.GetProcessInfo());
        r_target.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, target_laws, rTarget.GetProcessInfo());
        KRATOS_ERROR_IF(source_laws.size() != 4 || target_laws.size() != 4) << "Uninitialized laws or integration rule mismatch";
        for (std::size_t i = 0; i < 4; ++i) {
            KRATOS_ERROR_IF(!source_laws[i] || !target_laws[i]) << "Uninitialized law";
            const auto& r_law = *target_laws[i];
            const bool thermal = typeid(r_law) == typeid(ThermalLinearElastic2DPlaneStrain);
            KRATOS_ERROR_IF(typeid(r_law) != typeid(LinearPlaneStrain) && !thermal) << "Unsupported constitutive law";
            KRATOS_ERROR_IF(typeid(*source_laws[i]) != typeid(r_law)) << "Constitutive law mismatch";
            KRATOS_ERROR_IF(r_law.HasInitialState()) << "Destination already has initial state";
            if (thermal) {
                KRATOS_ERROR_IF(r_source.GetProperties()[THERMAL_EXPANSION] != r_target.GetProperties()[THERMAL_EXPANSION]) << "Material mismatch";
                for (const auto& r_node : r_geometry) {
                    KRATOS_ERROR_IF(r_node.FastGetSolutionStepValue(TEMPERATURE) !=
                        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE)) << "Stress-only probe requires initial T = T_ref";
                }
            }
        }
        for (const auto* p_variable : {&YOUNG_MODULUS, &POISSON_RATIO, &DENSITY, &THICKNESS}) {
            KRATOS_ERROR_IF(r_source.GetProperties()[*p_variable] != r_target.GetProperties()[*p_variable]) << "Material mismatch";
        }
    }
    for (auto& r_target : rTarget.Elements()) {
        std::vector<Vector> stresses;
        rSource.GetElement(r_target.Id()).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stresses, rSource.GetProcessInfo());
        std::vector<ConstitutiveLaw::Pointer> laws;
        r_target.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, rTarget.GetProcessInfo());
        for (std::size_t i = 0; i < laws.size(); ++i) {
            laws[i]->SetInitialState(Kratos::make_intrusive<InitialState>(
                ZeroVector(3), stresses[i], IdentityMatrix(2)));
        }
    }
}
}
PYBIND11_MODULE(dam_selective_state_probe, m)
{
    m.def("assign_stresses", &assign_stresses,
        "Test-only same-mesh Q4 stress assignment through constitutive InitialState.");
}
}
