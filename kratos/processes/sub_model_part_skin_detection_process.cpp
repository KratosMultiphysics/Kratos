//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// Project includes
#include "sub_model_part_skin_detection_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

template<SizeType TDim>
void SubModelPartSkinDetectionProcess<TDim>::SelectIfAllNodesOnSubModelPart::Prepare(ModelPart& rMainModelPart) const
{
    VariableUtils().SetFlag(SubModelPartSkinDetectionProcess::NODE_SELECTED, true, rMainModelPart.GetSubModelPart(mName).Nodes());
}

template<SizeType TDim>
bool SubModelPartSkinDetectionProcess<TDim>::SelectIfAllNodesOnSubModelPart::IsSelected(const Geometry<Node>::PointsArrayType& rNodes) const
{
    bool select = true;
    for (auto i_node = rNodes.begin(); i_node != rNodes.end(); ++i_node) {
        select &= i_node->Is(SubModelPartSkinDetectionProcess::NODE_SELECTED);
    }
    return select;
}

template<SizeType TDim>
void SubModelPartSkinDetectionProcess<TDim>::SelectIfOneNodeNotOnSubModelPart::Prepare(ModelPart& rMainModelPart) const
{
    VariableUtils().SetFlag(SubModelPartSkinDetectionProcess::NODE_SELECTED, true, rMainModelPart.Nodes());
    for (const auto& r_name : mNames) {
        VariableUtils().SetFlag(SubModelPartSkinDetectionProcess::NODE_SELECTED, false, rMainModelPart.GetSubModelPart(r_name).Nodes());
    }
}

template<SizeType TDim>
bool SubModelPartSkinDetectionProcess<TDim>::SelectIfOneNodeNotOnSubModelPart::IsSelected(const Geometry<Node>::PointsArrayType& rNodes) const
{
    for (const auto& r_node : rNodes) {
        if (r_node.Is(SubModelPartSkinDetectionProcess::NODE_SELECTED)) {
            return true;
        }
    }
    return false;
}

template<SizeType TDim>
SubModelPartSkinDetectionProcess<TDim>::SubModelPartSkinDetectionProcess(
    ModelPart& rModelPart, Parameters Settings)
    : SkinDetectionProcess<TDim>(rModelPart, Settings, this->GetDefaultParameters())
{
    KRATOS_TRY;

    Parameters settings = this->GetSettings();
    if (settings["selection_criteria"].GetString() == "nodes_on_sub_model_part")
    {
        KRATOS_ERROR_IF_NOT(settings["selection_settings"].Has("sub_model_part_name"))
        << "When using \"selection_criteria\" == \"nodes_on_sub_model_part\","
        << " SubModelPartSkinDetectionProcess requires the name of the target SubModelPart,"
        << " given as the \"sub_model_part_name\" string argument." << std::endl;
        mpFaceSelector = Kratos::make_shared<SelectIfAllNodesOnSubModelPart>(
            settings["selection_settings"]["sub_model_part_name"].GetString()
        );
   } else if (settings["selection_criteria"].GetString() == "node_not_on_sub_model_part") {
        KRATOS_ERROR_IF_NOT(settings["selection_settings"].Has("sub_model_part_names"))
        << "When using \"selection_criteria\" == \"node_not_on_sub_model_part\","
        << " SubModelPartSkinDetectionProcess requires the name of the target SubModelParts,"
        << " given as the \"sub_model_part_names\" string array argument." << std::endl;
        mpFaceSelector = Kratos::make_shared<SelectIfOneNodeNotOnSubModelPart>(
            settings["selection_settings"]["sub_model_part_names"].GetStringArray()
        );
    } else {
        KRATOS_ERROR << "Unsupported \"selection_criteria\" \"" << settings["selection_criteria"].GetString() << "\"." << std::endl;
    }

    KRATOS_CATCH("");
}

template<SizeType TDim>
void SubModelPartSkinDetectionProcess<TDim>::Execute()
{
    mpFaceSelector->Prepare(this->GetModelPart());
    SkinDetectionProcess<TDim>::Execute();
}

template<SizeType TDim>
void SubModelPartSkinDetectionProcess<TDim>::CreateConditions(
    ModelPart& rMainModelPart,
    ModelPart& rSkinModelPart,
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap,
    std::unordered_set<IndexType>& rNodesInTheSkin,
    const std::string& rConditionName) const
{
    IndexType condition_id = rMainModelPart.GetRootModelPart().Conditions().size();
    const auto& r_process_info = rMainModelPart.GetProcessInfo();

    // Create the auxiliar conditions
    for (auto& map : rInverseFaceMap) {
        condition_id += 1;

        const VectorIndexType& nodes_face = map.second;

        Properties::Pointer p_prop;
        const IndexType property_id = rPropertiesFaceMap[map.first];
         if (rMainModelPart.RecursivelyHasProperties(property_id)) {
             p_prop = rMainModelPart.pGetProperties(property_id);
         } else {
             p_prop = rMainModelPart.CreateNewProperties(property_id);
         }

        const std::string complete_name = rConditionName + std::to_string(TDim) + "D" + std::to_string(nodes_face.size()) + "N"; // If the condition doesn't follow this structure...sorry, we then need to modify this...

        Geometry<Node>::PointsArrayType condition_nodes;
        for (unsigned int i = 0; i < nodes_face.size(); i++)
        {
            condition_nodes.push_back(rMainModelPart.pGetNode(nodes_face[i]));
        }

        if (mpFaceSelector->IsSelected(condition_nodes))
        {
            auto p_cond = rMainModelPart.CreateNewCondition(complete_name, condition_id, condition_nodes, p_prop);
            rSkinModelPart.AddCondition(p_cond);
            p_cond->Set(INTERFACE, true);
            p_cond->Initialize(r_process_info);

            for (auto& index : nodes_face)
                rNodesInTheSkin.insert(index);
        }
    }
}

template<SizeType TDim>
const Parameters SubModelPartSkinDetectionProcess<TDim>::GetDefaultParameters() const
{
    Parameters defaults = SkinDetectionProcess<TDim>::GetDefaultParameters();
    defaults.AddEmptyValue("selection_criteria");
    defaults["selection_criteria"].SetString("");
    defaults.AddValue("selection_settings", Parameters("{}"));
    const Parameters const_defaults(defaults);
    return const_defaults;
}

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template<SizeType TDim>
const Kratos::Flags SubModelPartSkinDetectionProcess<TDim>::NODE_SELECTED(Kratos::Flags::Create(0));


template class SubModelPartSkinDetectionProcess<2>;
template class SubModelPartSkinDetectionProcess<3>;

}  // namespace Kratos.
