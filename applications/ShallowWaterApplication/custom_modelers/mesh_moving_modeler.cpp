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
    auto input_file_name = mParameters["input_file_name"].GetString();
    auto fixed_model_part_name = mParameters["fixed_model_part_name"].GetString();
    ModelPart& fixed_model_part = mpModel->GetModelPart(fixed_model_part_name);
    Flags io_options = IO::READ;
    if (mParameters["skip_timer"].GetBool())
        io_options = IO::SKIP_TIMER | io_options;
    if (mParameters["ignore_variables_not_in_solution_step_data"].GetBool())
        io_options = IO::IGNORE_VARIABLES_ERROR | io_options;
    ModelPartIO(input_file_name, io_options).ReadModelPart(fixed_model_part);
}

void MeshMovingModeler::PrepareGeometryModel()
{}

void MeshMovingModeler::SetupModelPart()
{
    auto fixed_model_part_name = mParameters["fixed_model_part_name"].GetString();
    auto moving_model_part_name = mParameters["moving_model_part_name"].GetString();
    ModelPart& fixed_model_part = mpModel->GetModelPart(fixed_model_part_name);
    ModelPart& moving_model_part = mpModel->GetModelPart(moving_model_part_name);
    const double relative_dry_height = mParameters["relative_dry_height"].GetDouble();
    ShallowWaterUtilities().IdentifyWetDomain(fixed_model_part, TO_COPY, relative_dry_height);

    // Both model parts share the same ProcessInfo
    moving_model_part.SetProcessInfo(fixed_model_part.pGetProcessInfo());

    // Nodes
    for (auto& r_node : fixed_model_part.Nodes())
    {
        if (r_node.Is(TO_COPY)) {
            auto p_new_node = moving_model_part.CreateNewNode(r_node.Id(), r_node);

            p_new_node->FastGetSolutionStepValue(HEIGHT) = r_node.FastGetSolutionStepValue(HEIGHT);
        }
    }

    // Properties
    for (auto& r_prop : fixed_model_part.rProperties())
    {
        auto new_prop = Kratos::make_shared<Properties>(r_prop);
        moving_model_part.AddProperties(new_prop);
    }

    // Elements: we need to clone the geometry with the new nodes
    std::size_t last_elem_id = 0;
    ModelPart::ElementsContainerType aux_array_with_element_pointers;
    for (auto& r_elem : fixed_model_part.Elements())
    {
        if (r_elem.Is(TO_COPY)) {
            NodesArrayType new_elem_nodes;
            for (auto& r_node : r_elem.GetGeometry())
            {
                new_elem_nodes.push_back(moving_model_part.pGetNode(r_node.Id()));
            }
            auto new_elem = r_elem.Clone(++last_elem_id, new_elem_nodes);
            new_elem->SetProperties(moving_model_part.pGetProperties(r_elem.GetProperties().Id()));
            aux_array_with_element_pointers.push_back(new_elem);
        }
    }
    moving_model_part.AddElements(aux_array_with_element_pointers.begin(), aux_array_with_element_pointers.end());

    // Conditions: we need to clone the geometry with the new nodes
    std::size_t last_cond_id = 0;
    ModelPart::ConditionsContainerType aux_array_with_condition_pointers;
    for (auto& r_cond : fixed_model_part.Conditions())
    {
        if (r_cond.Is(TO_COPY)) {
            NodesArrayType new_cond_nodes;
            for (auto& r_node : r_cond.GetGeometry())
            {
                new_cond_nodes.push_back(moving_model_part.pGetNode(r_node.Id()));
            }
            auto new_cond = r_cond.Clone(++last_cond_id, new_cond_nodes);
            new_cond->SetProperties(moving_model_part.pGetProperties(r_cond.GetProperties().Id()));
            aux_array_with_condition_pointers.push_back(new_cond);
        }
    }
    moving_model_part.AddConditions(aux_array_with_condition_pointers.begin(), aux_array_with_condition_pointers.end());
}

const Parameters MeshMovingModeler::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "input_file_name"                            : "",
        "fixed_model_part_name"                      : "eulerian",
        "moving_model_part_name"                     : "lagrangian",
        "skip_timer"                                 : true,
        "ignore_variables_not_in_solution_step_data" : false,
        "relative_dry_height"                        : 0.1
    })");
    return default_parameters;
}

}  // namespace Kratos.
