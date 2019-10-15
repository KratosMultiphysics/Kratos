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
#include "id_renumbering_process.h"


namespace Kratos
{

IdRenumberingProcess::IdRenumberingProcess(Model& rThisModel) : mrModel(rThisModel)
{
    mModelPartNames = GetRootModelPartNames();
}

IdRenumberingProcess::IdRenumberingProcess(Model& rThisModel, const StringVectorType& rModelPartNames)
 : mrModel(rThisModel)
{
    mModelPartNames = rModelPartNames;
}

void IdRenumberingProcess::RenumberNodes()
{
    // The new absolute unique_id
    IndexType unique_id = 0;
    mOriginNodesIdsMap.clear();

    // Loop the model parts
    for (const auto& name : mModelPartNames)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (IndexType i = 0; i < model_part.NumberOfNodes(); ++i)
        {
            auto it_node = model_part.NodesBegin() + i;
            mOriginNodesIdsMap[++unique_id] = it_node->Id();
            it_node->SetId(unique_id);
        }
    }
}

void IdRenumberingProcess::RenumberElements()
{
    // The new absolute unique_id
    IndexType unique_id = 0;
    mOriginElementsIdsMap.clear();

    // Loop the model parts
    for (const auto& name : mModelPartNames)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (IndexType i = 0; i < model_part.NumberOfElements(); ++i)
        {
            auto it_elem = model_part.ElementsBegin() + i;
            mOriginElementsIdsMap[++unique_id] = it_elem->Id();
            it_elem->SetId(unique_id);
        }
    }
}

void IdRenumberingProcess::RenumberConditions()
{
    // The new absolute unique_id
    IndexType unique_id = 0;
    mOriginConditionsIdsMap.clear();

    // Loop the model parts
    for (const auto& name : mModelPartNames)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (IndexType i = 0; i < model_part.NumberOfConditions(); ++i)
        {
            auto it_cond = model_part.ConditionsBegin() + i;
            mOriginConditionsIdsMap[++unique_id] = it_cond->Id();
            it_cond->SetId(unique_id);
        }
    }
}

void IdRenumberingProcess::RestoreNodes()
{
    if(mOriginNodesIdsMap.size() != 0)
    {
        // Loop the model parts
        for (const auto& name : mModelPartNames)
        {
            ModelPart& model_part = mrModel.GetModelPart(name);
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(model_part.NumberOfNodes()); ++i)
            {
                auto it_node = model_part.NodesBegin() + i;
                it_node->SetId(mOriginNodesIdsMap[it_node->Id()]);
            }
        }
    }
    // Avoiding the possibility to restore the ids twice
    mOriginNodesIdsMap.clear();
}

void IdRenumberingProcess::RestoreElements()
{
    if (mOriginElementsIdsMap.size() != 0)
    {
        // Loop the model parts
        for (const auto& name : mModelPartNames)
        {
            ModelPart& model_part = mrModel.GetModelPart(name);
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(model_part.NumberOfElements()); ++i)
            {
                auto it_elem = model_part.ElementsBegin() + i;
                it_elem->SetId(mOriginElementsIdsMap[it_elem->Id()]);
            }
        }
    }
    // Avoiding the possibility to restore the ids twice
    mOriginElementsIdsMap.clear();
}

void IdRenumberingProcess::RestoreConditions()
{
    if (mOriginConditionsIdsMap.size() != 0)
    {
        // Loop the model parts
        for (const auto& name : mModelPartNames)
        {
            ModelPart& model_part = mrModel.GetModelPart(name);
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(model_part.NumberOfConditions()); ++i)
            {
                auto it_cond = model_part.ConditionsBegin() + i;
                it_cond->SetId(mOriginConditionsIdsMap[it_cond->Id()]);
            }
        }
    }
    // Avoiding the possibility to restore the ids twice
    mOriginConditionsIdsMap.clear();
}

StringVectorType IdRenumberingProcess::GetRootModelPartNames() const
{
    StringVectorType root_names;
    StringVectorType all_names = mrModel.GetModelPartNames();
    for (const auto name : all_names)
    {
        // Get the root name of the model part.
        std::string root_name = mrModel.GetModelPart(name).GetRootModelPart().Name();

        // Add the root name if it is not present in the vector.
        if (std::find(root_names.begin(), root_names.end(), root_name) == root_names.end())
        {
            root_names.push_back(root_name);
        }
    }
    return root_names;
}

}  // namespace Kratos.
