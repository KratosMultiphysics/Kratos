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

    mOwnName = mParameters["own_model_part_name"].GetString();
    mRefinedName = mParameters["refined_model_part_name"].GetString();

    std::string own_name = mParameters["own_model_part_name"].GetString();
    std::string refined_name = mParameters["refined_model_part_name"].GetString();

    // Get the model part hierarchy
    StringVectorType sub_model_parts_names;
    if (mrRootModelPart.HasSubModelPart(own_name))
        sub_model_parts_names = mrRootModelPart.GetSubModelPart(own_name).GetSubModelPartNames();
        // sub_model_parts_names = RecursiveGetSubModelPartNames(mrRootModelPart.GetSubModelPart(own_name));
    else
        sub_model_parts_names = mrRootModelPart.GetSubModelPartNames();
        // sub_model_parts_names = RecursiveGetSubModelPartNames(mrRootModelPart);

    // Clone the model part at the own level
    InitializeOwnModelPart(sub_model_parts_names);

    // Initialize the refined model part
    InitializeRefinedModelPart(sub_model_parts_names);
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
    //         ModelPart::Pointer p_sub_model_part = new_model_part;
    //         KRATOS_WATCH(sub_full_name)
    //         std::istringstream iss(sub_full_name);
    //         std::string token;
    //         while (std::getline(iss, token, '.'))
    //         {
    //             if (p_sub_model_part->HasSubModelPart(token))
    //                 p_sub_model_part = p_sub_model_part->pGetSubModelPart(token);
    //             else
    //                 p_sub_model_part = p_sub_model_part->CreateSubModelPart(token);
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


void MultiScaleRefiningProcess::InitializeOwnModelPart(const StringVectorType& rNames)
{
    // Get the own model part
    if (mrRootModelPart.HasSubModelPart(mOwnName))
        mpOwnModelPart = mrRootModelPart.pGetSubModelPart(mOwnName);
    else
    {
        mpOwnModelPart = mrRootModelPart.CreateSubModelPart(mOwnName);
        
        // Copy all the tables and properties
        AddAllTablesToModelPart(mrRootModelPart, mpOwnModelPart);
        AddAllPropertiesToModelPart(mrRootModelPart, mpOwnModelPart);

        // Copy all the nodes, elements and conditions
        AddAllNodesToModelPart(mrRootModelPart, mpOwnModelPart);
        AddAllElementsToModelPart(mrRootModelPart, mpOwnModelPart);
        AddAllConditionsToModelPart(mrRootModelPart, mpOwnModelPart);
    
        // Copy the hierarchy from the root model part to the own model part
        for (auto name : rNames)
        {
            ModelPart::Pointer p_sub_model_part;
            if (mpOwnModelPart->HasSubModelPart(name))
                p_sub_model_part = mpOwnModelPart->pGetSubModelPart(name);
            else
                p_sub_model_part = mpOwnModelPart->CreateSubModelPart(name);

            ModelPart& origin_model_part = mrRootModelPart.GetSubModelPart(name);

            // Copy all the tables and properties
            AddAllTablesToModelPart(origin_model_part, p_sub_model_part);
            AddAllPropertiesToModelPart(origin_model_part, p_sub_model_part);

            // Copy all the nodes, elements and conditions
            AddAllNodesToModelPart(origin_model_part, p_sub_model_part);
            AddAllElementsToModelPart(origin_model_part, p_sub_model_part);
            AddAllConditionsToModelPart(origin_model_part, p_sub_model_part);
        }
    }
}


// TODO: remove this method
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

        // Copy all the tables and properties
        AddAllTablesToModelPart(mrRootModelPart, mpOwnModelPart);
        AddAllPropertiesToModelPart(mrRootModelPart, mpOwnModelPart);

        // Copy the nodes, elements and conditions
        AddAllNodesToModelPart(mrRootModelPart, mpOwnModelPart);
        AddAllElementsToModelPart(mrRootModelPart, mpOwnModelPart);
        AddAllConditionsToModelPart(mrRootModelPart, mpOwnModelPart);
    }

    // Copy the hierarchy from the root model part to the own model part
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
                aux_model_part = aux_model_part->CreateSubModelPart(token);
        }

        ModelPart& origin_model_part = RecursiveGetSubModelPart(mrRootModelPart, full_name);

        // Copy all the tables and properties
        AddAllTablesToModelPart(origin_model_part, aux_model_part);
        AddAllPropertiesToModelPart(origin_model_part, aux_model_part);

        // Copy the nodes, elements and conditions
        AddAllNodesToModelPart(origin_model_part, aux_model_part);
        AddAllElementsToModelPart(origin_model_part, aux_model_part);
        AddAllConditionsToModelPart(origin_model_part, aux_model_part);
    }
}


void MultiScaleRefiningProcess::InitializeRefinedModelPart(const StringVectorType& rNames)
{
    // Create the refined sub model part
    KRATOS_ERROR_IF(mrRootModelPart.HasSubModelPart(mRefinedName)) << "MultiScaleRefiningProcess: a refined model part with name : " << mRefinedName << " is already present in the model part : " << mrRootModelPart.Name() << std::endl;
    ModelPart::Pointer refined_model_part = mrRootModelPart.CreateSubModelPart(mRefinedName);
    mpRefinedModelPart = refined_model_part->CreateSubModelPart(mOwnName);

    // Copy all the tables and properties
    AddAllTablesToModelPart(mrRootModelPart, mpRefinedModelPart);
    AddAllPropertiesToModelPart(mrRootModelPart, mpRefinedModelPart);

    // Copy the hierarchy to the refined model part
    for (auto name : rNames)
    {
        ModelPart::Pointer p_sub_model_part = mpRefinedModelPart->CreateSubModelPart(name);

        ModelPart& origin_model_part = mrRootModelPart.GetSubModelPart(name);

        // Copy all the tables and properties
        AddAllTablesToModelPart(origin_model_part, p_sub_model_part);
        AddAllPropertiesToModelPart(origin_model_part, p_sub_model_part);
    
        // Note: we don't add the nodes, elements and conditions
        // This operation is the refining process itself
    }
}


// TODO: remove this method
void MultiScaleRefiningProcess::InitializeRefinedModelPart(
    const std::string& rRefinedName,
    const std::string& rOwnName,
    const StringVectorType& rNames
    )
{
    // Create the refined sub model part
    ModelPart::Pointer refined_model_part = mrRootModelPart.CreateSubModelPart(rRefinedName);
    mpRefinedModelPart = refined_model_part->CreateSubModelPart(rOwnName);

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


void MultiScaleRefiningProcess::AddAllPropertiesToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
{
    const IndexType nprop = rOriginModelPart.NumberOfProperties();
    ModelPart::PropertiesContainerType::iterator prop_begin = rOriginModelPart.PropertiesBegin();

    for (IndexType i = 0; i < nprop; i++)
    {
        auto prop = prop_begin + i;
        pDestinationModelPart->AddProperties(*prop.base());
    }
}


void MultiScaleRefiningProcess::AddAllTablesToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
{
    const IndexType ntables = rOriginModelPart.NumberOfTables();
    ModelPart::TablesContainerType::iterator table_begin = rOriginModelPart.TablesBegin();

    for (IndexType i = 0; i < ntables; i++)
    {
        auto table = table_begin + i;
        pDestinationModelPart->AddTable(table.base()->first, table.base()->second);
    }
}

void MultiScaleRefiningProcess::AddAllNodesToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
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


void MultiScaleRefiningProcess::AddAllElementsToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
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


void MultiScaleRefiningProcess::AddAllConditionsToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart)
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