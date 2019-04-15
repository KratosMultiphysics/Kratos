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

IdRenumberingProcess::IdRenumberingProcess(Model& rThisModel) : mrModel(rThisModel){}


void IdRenumberingProcess::RenumberNodes()
{
    // The new absolute id
    IndexType id = 0;
    mOriginNodesIds.clear();

    // Loop all the model parts
    StringVectorType model_part_names = mrModel.GetModelPartNames();
    for (auto& name : model_part_names)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfNodes()); ++i)
        {
            ModelPart::NodeIterator it_node = model_part.NodesBegin() + i;
            mOriginNodesIds[++id] = it_node->Id();
            it_node->SetId(id);
        }
    }
}

void IdRenumberingProcess::RenumberElements()
{
    // The new absolute id
    IndexType id = 0;
    mOriginElementsIds.clear();

    // Loop all the model parts
    StringVectorType model_part_names = mrModel.GetModelPartNames();
    for (auto& name : model_part_names)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfElements()); ++i)
        {
            ModelPart::ElementIterator it_elem = model_part.ElementsBegin() + i;
            mOriginElementsIds[++id] = it_elem->Id();
            it_elem->SetId(id);
        }
    }
}

void IdRenumberingProcess::RenumberConditions()
{
    // The new absolute id
    IndexType id = 0;
    mOriginConditionsIds.clear();

    // Loop all the model parts
    StringVectorType model_part_names = mrModel.GetModelPartNames();
    for (auto& name : model_part_names)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfConditions()); ++i)
        {
            ModelPart::ConditionIterator it_cond = model_part.ConditionsBegin() + i;
            mOriginConditionsIds[++id] = it_cond->Id();
            it_cond->SetId(id);
        }
    }
}

void IdRenumberingProcess::RestoreNodes()
{
    // Loop all the model parts
    StringVectorType model_part_names = mrModel.GetModelPartNames();
    for (auto& name : model_part_names)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfNodes()); ++i)
        {
            ModelPart::NodeIterator it_node = model_part.NodesBegin() + i;
            it_node->SetId(mOriginNodesIds[it_node->Id()]);
        }
    }
}

void IdRenumberingProcess::RestoreElements()
{
    // Loop all the model parts
    StringVectorType model_part_names = mrModel.GetModelPartNames();
    for (auto& name : model_part_names)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfElements()); ++i)
        {
            ModelPart::ElementIterator it_elem = model_part.ElementsBegin() + i;
            it_elem->SetId(mOriginElementsIds[it_elem->Id()]);
        }
    }
}

void IdRenumberingProcess::RestoreConditions()
{
    // Loop all the model parts
    StringVectorType model_part_names = mrModel.GetModelPartNames();
    for (auto& name : model_part_names)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfConditions()); ++i)
        {
            ModelPart::ConditionIterator it_cond = model_part.ConditionsBegin() + i;
            it_cond->SetId(mOriginConditionsIds[it_cond->Id()]);
        }
    }
}

}  // namespace Kratos.
