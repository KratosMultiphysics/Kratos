// BSD License, Kratos default license: kratos/license.txt
#include <algorithm>
#include <set>
#include "custom_utilities/scalar_gap_geometry_utility.h"
#include "convection_diffusion_application_variables.h"
#include "utilities/multipatch_gap_utilities.h"

namespace Kratos {
void ScalarGapGeometryUtility::Create(ModelPart& rPhase, ModelPart& rComputing, const Matrix& rContour,
                                      Properties::Pointer pConditionProperties, double Tolerance)
{
    KRATOS_TRY
    auto& r_root = rPhase.GetRootModelPart();
    KRATOS_ERROR_IF(&r_root != &rComputing.GetRootModelPart() || &r_root == &rComputing)
        << "Scalar gap model parts must share a root and computing must be a submodelpart." << std::endl;
    KRATOS_ERROR_IF(rComputing.NumberOfElements() || rComputing.NumberOfConditions() || rComputing.NumberOfSubModelParts())
        << "Scalar gap computing model part must be empty." << std::endl;
    const auto data = FemMultipatchGapUtilities::BuildLevelSetPhase(rPhase, rContour, false, Tolerance);
    auto& r_volume = rComputing.CreateSubModelPart("GapVolumes");
    auto& r_interface = rComputing.CreateSubModelPart("PhysicalInterfaces");
    auto& r_connectors = rComputing.CreateSubModelPart("GapInterfaces");
    ModelPart::IndexType element_id = 1, condition_id = 1;
    for (const auto& r_element : r_root.Elements()) element_id = std::max(element_id, r_element.Id() + 1);
    for (const auto& r_condition : r_root.Conditions()) condition_id = std::max(condition_id, r_condition.Id() + 1);
    for (auto& r_element : rPhase.Elements()) rComputing.AddElement(rPhase.pGetElement(r_element.Id()));
    for (const auto& r_cell : data.Cells) {
        auto p_element = r_cell.pParent->Create(element_id++, r_cell.pQuadrature, r_cell.pParent->pGetProperties());
        r_volume.AddElement(p_element);
    }
    const auto create = [&](const FemMultipatchGapUtilities::Segment& rSegment, ModelPart& rPart, bool Connector) {
        std::vector<ModelPart::IndexType> ids;
        std::set<ModelPart::IndexType> inserted;
        GlobalPointersVector<Element> parents;
        for (auto p_parent : {rSegment.pFirst, rSegment.pSecond}) {
            if (!p_parent) continue;
            parents.push_back(GlobalPointer<Element>(p_parent.get()));
            for (const auto& r_node : p_parent->GetGeometry())
                if (inserted.insert(r_node.Id()).second) ids.push_back(r_node.Id());
        }
        auto p_condition = rPart.CreateNewCondition(Connector ? "LaplacianGapConnectorCondition2D" :
                                                    "LaplacianGapDirichletCondition2D", condition_id++, ids,
                                                    pConditionProperties);
        p_condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
        Matrix coordinates = ZeroMatrix(2, 3);
        for (std::size_t d = 0; d < 3; ++d) {
            coordinates(0, d) = rSegment.Start[d];
            coordinates(1, d) = rSegment.End[d];
        }
        p_condition->SetValue(Connector ? SURROGATE_BOUNDARY_FACE_COORDINATES : GAP_INTERFACE_COORDINATES,
                              coordinates);
    };
    for (const auto& r_segment : data.Physical) create(r_segment, r_interface, false);
    for (const auto& r_segment : data.Connectors) create(r_segment, r_connectors, true);
    std::set<ModelPart::IndexType> node_ids;
    for (const auto& r_element : rComputing.Elements())
        for (const auto& r_node : r_element.GetGeometry()) node_ids.insert(r_node.Id());
    for (const auto& r_condition : rComputing.Conditions())
        for (const auto& r_node : r_condition.GetGeometry()) node_ids.insert(r_node.Id());
    for (const auto id : node_ids) rComputing.AddNode(r_root.pGetNode(id));
    KRATOS_CATCH("")
}
}
