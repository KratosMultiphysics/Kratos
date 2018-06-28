//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "custom_processes/multi_scale_refining_process.h"
#include "utilities/sub_model_parts_list_utility.h"
#include "custom_utilities/uniform_refine_utility.h"

namespace Kratos
{

MultiScaleRefiningProcess::MultiScaleRefiningProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters)
    : mrRootModelPart(rThisModelPart)
    , mParameters(ThisParameters)
{
    Parameters DefaultParameters = Parameters(R"(
    {
        "model_part_name"              : "MainModelPart",
        "own_model_part_name"          : "own",
        "refined_model_part_name"      : "refined",
        "echo_level"                   : 0,
        "number_of_divisions_at_level" : 2,
        "refining_boundary_condition"  : "Condition2D3N"
    }
    )");

    mParameters.ValidateAndAssignDefaults(DefaultParameters);

    mEchoLevel = mParameters["echo_level"].GetInt();
    mDivisions = mParameters["number_of_divisions_at_level"].GetInt();
    mConditionName = mParameters["refining_boundary_condition"].GetString();

    std::string own_name = mParameters["own_name"].GetString();
    std::string refined_name = mParameters["refined_name"].GetString();

    // Get the model part hierarchy
    StringVectorType sub_model_parts_names;
    if (mrRootModelPart.HasSubModelPart(own_name))
        sub_model_parts_names = RecursiveGetSubModelPartNames(mrRootModelPart.GetSubModelPart(own_name));
    else
        sub_model_parts_names = RecursiveGetSubModelPartNames(mrRootModelPart);

    // Clone the model part at the own level
    InitializeOwnModelPart(own_name, sub_model_parts_names);

    // Initialize the refined model part
    InitializeRefinedModelPart(refined_name, own_name, sub_model_parts_names);
}


void MultiScaleRefiningProcess::InitializeModelPart()
{
    // ModelPart& root_mp = mrModel.GetModelPart(mModelPartsNames[0]);
    
    // StringVectorType sub_model_parts_names = root_mp.GetSubModelPartNames();
    StringVectorType sub_model_parts_names = this->RecursiveGetSubModelPartNames(mrRootModelPart);

    // for (unsigned int i = 1; i <= mMaxMultiScaleLevels; i++)
    // {
    //     std::string name = mModelPartsNames[0] + "-level_" + std::to_string(i);
    //     KRATOS_INFO_IF(this->Info(), mEchoLevel>1) << "Creating Model Part " << name << std::endl;
    //     mModelPartsNames.push_back(name);

    //     ModelPart::Pointer new_model_part = Kratos::make_shared<ModelPart>( name );
    //     mrModel.AddModelPart(new_model_part);
    //     KRATOS_WATCH(mrModel)
    //     // TODO: using the new Model interface:
    //     // mrModel.CreateModelPart(name);

    //     for (auto sub_full_name : sub_model_parts_names)
    //     {
    //         ModelPart::Pointer sub_model_part = new_model_part;
    //         KRATOS_WATCH(sub_full_name)
    //         std::istringstream iss(sub_full_name);
    //         std::string token;
    //         while (std::getline(iss, token, '.'))
    //         {
    //             if (sub_model_part->HasSubModelPart(token))
    //                 sub_model_part = sub_model_part->pGetSubModelPart(token);
    //             else
    //                 sub_model_part = sub_model_part->CreateSubModelPart(token);
    //         }
    //         KRATOS_WATCH(*new_model_part)
    //     }
    // }

}


MultiScaleRefiningProcess::StringVectorType MultiScaleRefiningProcess::RecursiveGetSubModelPartNames(
    ModelPart& rThisModelPart,
    std::string Prefix
    )
{
    StringVectorType names = rThisModelPart.GetSubModelPartNames();
    if (!Prefix.empty())
        Prefix += ".";
    
    for (auto& name : names)
    {
        ModelPart& sub_model_part = rThisModelPart.GetSubModelPart(name);
        auto sub_names = this->RecursiveGetSubModelPartNames(sub_model_part, Prefix + name);
        name.insert(0, Prefix);
        for (auto sub_name : sub_names)
            names.push_back(sub_name);
    }

    return names;
}


ModelPart& MultiScaleRefiningProcess::RecursiveGetSubModelPart(ModelPart& rThisModelPart, std::string FullName)
{
    std::istringstream iss(FullName);
    std::string token;
    if (std::getline(iss, token, '.'))
    {
        ModelPart& aux_model_part = rThisModelPart.GetSubModelPart(token);
        return RecursiveGetSubModelPart(aux_model_part, iss.str());
    }
    return rThisModelPart;
}


void MultiScaleRefiningProcess::InitializeOwnModelPart(
    const std::string& rOwnName,
    const StringVectorType& rNames
    )
{
    // Get the own model part
    if (mrRootModelPart.HasSubModelPart(rOwnName))
        mpOwnModelPart = mrRootModelPart.pGetSubModelPart(rOwnName);
    else
    {
        mpOwnModelPart = mrRootModelPart.CreateSubModelPart(rOwnName);
        // Copy the nodes, elements and conditions
        AddNodesToSubModelPart(mrRootModelPart, mpOwnModelPart);
        AddElementsToSubModelPart(mrRootModelPart, mpOwnModelPart);
        AddConditionsToSubModelPart(mrRootModelPart, mpOwnModelPart);
    }

    // Copy the hierarchy to the own model part
    for (auto full_name : rNames)
    {
        ModelPart::Pointer aux_model_part = mpOwnModelPart;
        std::istringstream iss(full_name);
        std::string token;
        while (std::getline(iss, token, '.'))
        {
            if (aux_model_part->HasSubModelPart(token))
                aux_model_part = aux_model_part->pGetSubModelPart(token);
            else
            {
                aux_model_part = aux_model_part->CreateSubModelPart(token);
            }
        }
        // Copy all the nodes, elements and conditions
        ModelPart& origin_model_part = RecursiveGetSubModelPart(mrRootModelPart, full_name);
        AddNodesToSubModelPart(origin_model_part, aux_model_part);
        AddElementsToSubModelPart(origin_model_part, aux_model_part);
        AddConditionsToSubModelPart(origin_model_part, aux_model_part);
    }
}


void MultiScaleRefiningProcess::InitializeRefinedModelPart(
    const std::string& rRefinedName,
    const std::string& rOwnName,
    const StringVectorType& rNames
    )
{
    // Create the refined sub model part
    ModelPart::Pointer refined_model_part = mrRootModelPart.CreateSubModelPart(rRefinedName);
    mpRefinedModelPart = refined_model_part->pGetSubModelPart(rOwnName);

    // Copy the hierarchy to the refined model part
    for (auto full_name : rNames)
    {
        ModelPart::Pointer aux_model_part = mpRefinedModelPart;
        std::istringstream iss(full_name);
        std::string token;
        while (std::getline(iss, token, '.'))
        {
            if (aux_model_part->HasSubModelPart(token))
                aux_model_part = aux_model_part->pGetSubModelPart(token);
            else
            {
                aux_model_part = aux_model_part->CreateSubModelPart(token);
            }
        }
    }
}


void MultiScaleRefiningProcess::AddNodesToSubModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
{
    const int nnodes = static_cast<int>(rOriginModelPart.Nodes().size());
    IndexVectorType origin_nodes(nnodes);
    ModelPart::NodesContainerType::iterator node_begin = rOriginModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i < nnodes; i++)
    {
        auto node = node_begin + i;
        origin_nodes[i] = node->Id();
    }
    pDestinationModelPart->AddNodes(origin_nodes);
}


void MultiScaleRefiningProcess::AddElementsToSubModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
{
    const int nelems = static_cast<int>(rOriginModelPart.Elements().size());
    IndexVectorType origin_elems(nelems);
    ModelPart::ElementsContainerType::iterator elem_begin = rOriginModelPart.ElementsBegin();

    #pragma omp parallel for
    for (int i = 0; i < nelems; i++)
    {
        auto elem = elem_begin + i;
        origin_elems[i] = elem->Id();
    }
    pDestinationModelPart->AddElements(origin_elems);
}


void MultiScaleRefiningProcess::AddConditionsToSubModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
{
    const int nconds = static_cast<int>(rOriginModelPart.Conditions().size());
    IndexVectorType origin_conds(nconds);
    ModelPart::ConditionsContainerType::iterator cond_begin = rOriginModelPart.ConditionsBegin();

    #pragma omp parallel for
    for (int i = 0; i < nconds; i++)
    {
        auto cond = cond_begin + i;
        origin_conds[i] = cond->Id();
    }
    pDestinationModelPart->AddConditions(origin_conds);
}

} // namespace Kratos