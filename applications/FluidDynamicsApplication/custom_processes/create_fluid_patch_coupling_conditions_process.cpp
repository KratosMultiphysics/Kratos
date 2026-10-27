// License: BSD License (kratos/license.txt)

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

// Project includes
#include "custom_processes/create_fluid_patch_coupling_conditions_process.h"
#include "custom_conditions/fluid_patch_coupling_condition.h"

namespace Kratos
{
namespace
{
struct BoundaryEdge
{
    Node::Pointer mpFirst;
    Node::Pointer mpSecond;
    Element::Pointer mpParent;
    std::size_t mMultiplicity = 0;
};

std::vector<BoundaryEdge> extract_boundary_edges(ModelPart& rPatch)
{
    std::map<std::array<ModelPart::IndexType, 2>, BoundaryEdge> edge_map;
    for (auto& r_element : rPatch.Elements()) {
        if (r_element.IsDefined(ACTIVE) && r_element.IsNot(ACTIVE)) {
            continue;
        }
        auto& r_geometry = r_element.GetGeometry();
        KRATOS_ERROR_IF(r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Triangle2D3 &&
                        r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
            << "FEM fluid coupling requires Triangle2D3 (P1) or Quadrilateral2D4 (Q1) elements." << std::endl;
        KRATOS_ERROR_IF(r_geometry.Area() <= 0.0) << "Parent elements must have positive area." << std::endl;
        for (const auto& r_node : r_geometry) {
            KRATOS_ERROR_IF(std::abs(r_node.Z()) > 1e-10)
                << "FEM fluid coupling currently requires meshes in the z=0 plane." << std::endl;
        }
        for (auto& r_edge : r_geometry.GenerateEdges()) {
            std::array<ModelPart::IndexType, 2> key{r_edge[0].Id(), r_edge[1].Id()};
            std::sort(key.begin(), key.end());
            auto& r_data = edge_map[key];
            ++r_data.mMultiplicity;
            KRATOS_ERROR_IF(r_data.mMultiplicity > 2) << "Non-manifold patch edge." << std::endl;
            r_data.mpFirst = r_edge.pGetPoint(0);
            r_data.mpSecond = r_edge.pGetPoint(1);
            r_data.mpParent = rPatch.pGetElement(r_element.Id());
        }
    }
    std::vector<BoundaryEdge> result;
    for (const auto& [key, r_edge] : edge_map) {
        if (r_edge.mMultiplicity == 1) {
            result.push_back(r_edge);
        }
    }
    return result;
}
} // unnamed namespace

CreateFluidPatchCouplingConditionsProcess::CreateFluidPatchCouplingConditionsProcess(
    Model& rModel, Parameters Settings)
    : mrModel(rModel), mSettings(Settings.Clone())
{
    mSettings.ValidateAndAssignDefaults(GetDefaultParameters());
}

const Parameters CreateFluidPatchCouplingConditionsProcess::GetDefaultParameters() const
{
    return Parameters(R"({
        "first_model_part_name": "",
        "second_model_part_name": "",
        "interface_model_part_name": "",
        "properties_id": 0,
        "tolerance": 1e-10
    })");
}

void CreateFluidPatchCouplingConditionsProcess::Execute()
{
    KRATOS_TRY
    const double tolerance = mSettings["tolerance"].GetDouble();
    KRATOS_ERROR_IF(!std::isfinite(tolerance) || tolerance <= 0.0)
        << "Coupling tolerance must be finite and positive." << std::endl;
    auto& r_first = mrModel.GetModelPart(mSettings["first_model_part_name"].GetString());
    auto& r_second = mrModel.GetModelPart(mSettings["second_model_part_name"].GetString());
    auto& r_root = r_first.GetRootModelPart();
    KRATOS_ERROR_IF(&r_root != &r_second.GetRootModelPart())
        << "Coupled patches must belong to the same root ModelPart." << std::endl;
    const auto properties_id = mSettings["properties_id"].GetInt();
    KRATOS_ERROR_IF(properties_id < 0 || !r_root.HasProperties(properties_id))
        << "Coupling Properties must already exist in the root ModelPart." << std::endl;
    auto p_properties = r_root.pGetProperties(properties_id);
    KRATOS_ERROR_IF_NOT(p_properties->Has(PENALTY_COEFFICIENT))
        << "Coupling Properties require PENALTY_COEFFICIENT." << std::endl;

    std::set<ModelPart::IndexType> first_node_ids;
    for (const auto& r_node : r_first.Nodes()) {
        first_node_ids.insert(r_node.Id());
    }
    for (const auto& r_node : r_second.Nodes()) {
        KRATOS_ERROR_IF(first_node_ids.count(r_node.Id()) != 0)
            << "Coupled patches must have independent nodes; do not merge their interface nodes." << std::endl;
    }

    const auto first_edges = extract_boundary_edges(r_first);
    const auto second_edges = extract_boundary_edges(r_second);
    struct CouplingSegment
    {
        Element::Pointer mpFirst;
        Element::Pointer mpSecond;
        Matrix mCoordinates;
    };
    std::vector<CouplingSegment> segments;
    for (const auto& r_a : first_edges) {
        const double ax = r_a.mpFirst->X();
        const double ay = r_a.mpFirst->Y();
        const double dx = r_a.mpSecond->X() - ax;
        const double dy = r_a.mpSecond->Y() - ay;
        const double length = std::hypot(dx, dy);
        const double tx = dx / length;
        const double ty = dy / length;
        for (const auto& r_b : second_edges) {
            // Collinearity is tested in physical distance, not in a local
            // parameter. Point contacts have zero measure and are skipped.
            const double distance_first = tx * (r_b.mpFirst->Y() - ay) - ty * (r_b.mpFirst->X() - ax);
            const double distance_second = tx * (r_b.mpSecond->Y() - ay) - ty * (r_b.mpSecond->X() - ax);
            if (std::abs(distance_first) > tolerance || std::abs(distance_second) > tolerance) {
                continue;
            }
            const double b0 = tx * (r_b.mpFirst->X() - ax) + ty * (r_b.mpFirst->Y() - ay);
            const double b1 = tx * (r_b.mpSecond->X() - ax) + ty * (r_b.mpSecond->Y() - ay);
            const double start = std::max(0.0, std::min(b0, b1));
            const double end = std::min(length, std::max(b0, b1));
            if (end - start <= tolerance) {
                continue;
            }
            const auto center_a = r_a.mpParent->GetGeometry().Center();
            const auto center_b = r_b.mpParent->GetGeometry().Center();
            KRATOS_ERROR_IF(ty * (center_b.X() - center_a.X()) - tx * (center_b.Y() - center_a.Y()) <= 0.0)
                << "Interface parents must occupy opposite sides of their common edge." << std::endl;
            KRATOS_ERROR_IF(r_a.mpParent->GetGeometry().GetGeometryType() !=
                            r_b.mpParent->GetGeometry().GetGeometryType())
                << "Coupled parents must have the same type (P1/P1 or Q1/Q1)." << std::endl;
            Matrix coordinates = ZeroMatrix(2, 3);
            coordinates(0, 0) = ax + start * tx;
            coordinates(0, 1) = ay + start * ty;
            coordinates(1, 0) = ax + end * tx;
            coordinates(1, 1) = ay + end * ty;
            segments.push_back({r_a.mpParent, r_b.mpParent, coordinates});
        }
    }
    KRATOS_ERROR_IF(segments.empty())
        << "No common positive-length boundary segments were found." << std::endl;

    const auto interface_name = mSettings["interface_model_part_name"].GetString();
    KRATOS_ERROR_IF(interface_name.empty()) << "interface_model_part_name must be provided." << std::endl;
    auto& r_interface = mrModel.HasModelPart(interface_name)
        ? mrModel.GetModelPart(interface_name) : mrModel.CreateModelPart(interface_name);
    KRATOS_ERROR_IF(&r_interface.GetRootModelPart() != &r_root || &r_interface == &r_root)
        << "The interface must be a submodelpart of the patches' root." << std::endl;
    KRATOS_ERROR_IF(r_interface.NumberOfConditions() != 0 || r_interface.NumberOfElements() != 0)
        << "The destination interface must have no conditions or elements; execute the creator only once." << std::endl;

    ModelPart::IndexType next_id = 1;
    for (const auto& r_condition : r_root.Conditions()) {
        next_id = std::max(next_id, r_condition.Id() + 1);
    }
    for (const auto& r_segment : segments) {
        std::vector<ModelPart::IndexType> node_ids;
        for (const auto* p_parent : {r_segment.mpFirst.get(), r_segment.mpSecond.get()}) {
            for (const auto& r_node : p_parent->GetGeometry()) {
                node_ids.push_back(r_node.Id());
            }
        }
        const std::string condition_name = node_ids.size() == 6
            ? "FluidPatchCouplingCondition2D6N" : "FluidPatchCouplingCondition2D8N";
        auto p_condition = r_interface.CreateNewCondition(
            condition_name, next_id++, node_ids, p_properties);
        GlobalPointersVector<Element> parents;
        parents.resize(2);
        // The intrusive-pointer overload takes a non-const reference. Use
        // the addresses of the parent elements, which remain model-owned.
        parents(0) = GlobalPointer<Element>(r_segment.mpFirst.get());
        parents(1) = GlobalPointer<Element>(r_segment.mpSecond.get());
        p_condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
        p_condition->SetValue(FEM_COUPLING_INTERFACE_COORDINATES, r_segment.mCoordinates);
        r_interface.AddNodes(node_ids);
        p_condition->Set(ACTIVE, true);
    }
    KRATOS_CATCH("")
}

} // namespace Kratos
