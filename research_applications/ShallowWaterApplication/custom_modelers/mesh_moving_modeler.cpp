//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "mesh_moving_modeler.h"
#include "includes/model_part_io.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "processes/sub_model_part_skin_detection_process.h"


namespace Kratos
{

KRATOS_CREATE_LOCAL_FLAG(MeshMovingModeler, TO_COPY, 1);

MeshMovingModeler::MeshMovingModeler() : Modeler() {}

MeshMovingModeler::MeshMovingModeler(Model& rModel, Parameters ModelerParameters)
    : Modeler(rModel, ModelerParameters)
    , mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

void MeshMovingModeler::SetupGeometryModel()
{
    const auto input_file_name = mParameters["input_file_name"].GetString();
    const auto fixed_model_part_name = mParameters["fixed_model_part_name"].GetString();
    ModelPart& fixed_model_part = mpModel->GetModelPart(fixed_model_part_name);
    Flags io_options = IO::READ;
    if (mParameters["skip_timer"].GetBool())
        io_options = IO::SKIP_TIMER | io_options;
    if (mParameters["ignore_variables_not_in_solution_step_data"].GetBool())
        io_options = IO::IGNORE_VARIABLES_ERROR | io_options;
    ModelPartIO(input_file_name, io_options).ReadModelPart(fixed_model_part);

    // Both model parts share the same ProcessInfo
    const auto moving_model_part_name = mParameters["moving_model_part_name"].GetString();
    ModelPart& moving_model_part = mpModel->GetModelPart(moving_model_part_name);
    fixed_model_part.SetProcessInfo(moving_model_part.pGetProcessInfo());
}

void MeshMovingModeler::PrepareGeometryModel()
{}

void MeshMovingModeler::SetupModelPart()
{
    const auto fixed_model_part_name = mParameters["fixed_model_part_name"].GetString();
    const auto moving_model_part_name = mParameters["moving_model_part_name"].GetString();
    ModelPart& fixed_model_part = mpModel->GetModelPart(fixed_model_part_name);
    ModelPart& moving_model_part = mpModel->GetModelPart(moving_model_part_name);
    const double relative_dry_height = mParameters["relative_dry_height"].GetDouble();
    ShallowWaterUtilities().FlagWetElements(fixed_model_part, TO_COPY, relative_dry_height);
    ShallowWaterUtilities().ExtrapolateElementalFlagToNodes(fixed_model_part, TO_COPY);

    // Nodes
    for (auto& r_node : fixed_model_part.Nodes())
    {
        if (r_node.Is(TO_COPY)) {
            auto p_new_node = moving_model_part.CreateNewNode(r_node.Id(), r_node);
            p_new_node->FastGetSolutionStepValue(HEIGHT) = r_node.FastGetSolutionStepValue(HEIGHT);
        }
    }

    // Properties
    for (auto& r_prop : fixed_model_part.rProperties()) {
        moving_model_part.AddProperties(Kratos::make_shared<Properties>(r_prop));
    }

    // Elements: we need to clone the geometry with the new nodes
    ModelPart::ElementsContainerType aux_array_with_element_pointers;
    for (auto& r_elem : fixed_model_part.Elements())
    {
        if (r_elem.Is(TO_COPY)) {
            NodesArrayType new_elem_nodes;
            for (auto& r_node : r_elem.GetGeometry())
            {
                new_elem_nodes.push_back(moving_model_part.pGetNode(r_node.Id()));
            }
            auto new_elem = r_elem.Clone(r_elem.Id(), new_elem_nodes);
            auto new_property = moving_model_part.pGetProperties(r_elem.GetProperties().Id());
            new_elem->SetProperties(new_property);
            aux_array_with_element_pointers.push_back(new_elem);
        }
    }
    moving_model_part.AddElements(aux_array_with_element_pointers.begin(), aux_array_with_element_pointers.end());

    // Conditions: we need to clone the geometry with the new nodes
    ModelPart::ConditionsContainerType aux_array_with_condition_pointers;
    for (auto& r_cond : fixed_model_part.Conditions())
    {
        if (r_cond.Is(TO_COPY)) {
            NodesArrayType new_cond_nodes;
            for (auto& r_node : r_cond.GetGeometry())
            {
                new_cond_nodes.push_back(moving_model_part.pGetNode(r_node.Id()));
            }
            auto new_cond = r_cond.Clone(r_cond.Id(), new_cond_nodes);
            auto new_property = moving_model_part.pGetProperties(r_cond.GetProperties().Id());
            new_cond->SetProperties(new_property);
            aux_array_with_condition_pointers.push_back(new_cond);
        }
    }
    moving_model_part.AddConditions(aux_array_with_condition_pointers.begin(), aux_array_with_condition_pointers.end());

    // SubModelParts
    for (const auto& sub_model_part : fixed_model_part.SubModelParts())
    {
        auto& dest_sub_model_part = moving_model_part.CreateSubModelPart(sub_model_part.Name());

        std::vector<std::size_t> nodes_ids;
        std::vector<std::size_t> elements_ids;
        std::vector<std::size_t> conditions_ids;

        for (auto& r_node : sub_model_part.Nodes()) {
            if (r_node.Is(TO_COPY)) {
                nodes_ids.push_back(r_node.Id());
            }
        }
        for (auto& r_elem : sub_model_part.Elements()) {
            if (r_elem.Is(TO_COPY)) {
                elements_ids.push_back(r_elem.Id());
            }
        }
        for (auto& r_cond : sub_model_part.Conditions()) {
            if (r_cond.Is(TO_COPY)) {
                conditions_ids.push_back(r_cond.Id());
            }
        }

        dest_sub_model_part.AddNodes(nodes_ids);
        dest_sub_model_part.AddElements(elements_ids);
        dest_sub_model_part.AddConditions(conditions_ids);
    }

    // Generate the conditions for the shoreline
    auto skin_settings = Parameters();
    skin_settings.AddValue("name_auxiliar_model_part", mParameters["interface_sub_model_part_name"]);
    skin_settings.AddString("selection_criteria", "node_not_on_sub_model_part");
    skin_settings.AddEmptyValue("selection_settings");
    skin_settings["selection_settings"].AddValue("sub_model_part_names", mParameters["solid_boundary_sub_model_part_names"]);
    SubModelPartSkinDetectionProcess<2>(moving_model_part, skin_settings).Execute();

    // Set a serial numbering for the two model parts
    ShallowWaterUtilities().OffsetIds(moving_model_part.Nodes(), fixed_model_part.Nodes().size());
    ShallowWaterUtilities().OffsetIds(moving_model_part.Elements(), fixed_model_part.Elements().size());
    ShallowWaterUtilities().OffsetIds(moving_model_part.Conditions(), fixed_model_part.Conditions().size());
    ShallowWaterUtilities().OffsetIds(moving_model_part.rProperties(), fixed_model_part.rProperties().size());
}

const Parameters MeshMovingModeler::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "input_file_name"                            : "",
        "fixed_model_part_name"                      : "eulerian",
        "moving_model_part_name"                     : "lagrangian",
        "interface_sub_model_part_name"              : "shoreline",
        "solid_boundary_sub_model_part_names"        : [],
        "skip_timer"                                 : true,
        "ignore_variables_not_in_solution_step_data" : false,
        "relative_dry_height"                        : 0.1
    })");
    return default_parameters;
}

}  // namespace Kratos.
