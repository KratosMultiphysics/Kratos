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
#include "custom_processes/multi_scale_refining_process.h"
#include "geometries/point.h"
#include "processes/fast_transfer_between_model_parts_process.h"
#include "utilities/sub_model_parts_list_utility.h"
#include "custom_utilities/meshing_flags.h"

namespace Kratos
{

MultiScaleRefiningProcess::MultiScaleRefiningProcess(
    ModelPart& rThisCoarseModelPart,
    ModelPart& rThisRefinedModelPart,
    ModelPart& rThisVisualizationModelPart,
    Parameters ThisParameters)
    : mrCoarseModelPart(rThisCoarseModelPart)
    , mrRefinedModelPart(rThisRefinedModelPart)
    , mrVisualizationModelPart(rThisVisualizationModelPart)
    , mParameters(ThisParameters)
    , mUniformRefinement(mrRefinedModelPart)
{
    Parameters DefaultParameters = Parameters(R"(
    {
        "number_of_divisions_at_subscale"     : 2,
        "echo_level"                          : 0,
        "subscale_interface_base_name"        : "refined_interface",
        "subscale_boundary_condition"         : "Condition2D2N"
    }
    )");

    mParameters.ValidateAndAssignDefaults(DefaultParameters);

    mDivisionsAtSubscale = mParameters["number_of_divisions_at_subscale"].GetInt();
    mEchoLevel = mParameters["echo_level"].GetInt();

    std::string interface_base_name = mParameters["subscale_interface_base_name"].GetString();
    mRefinedInterfaceName = interface_base_name + "_" + std::to_string(mrCoarseModelPart.GetValue(SUBSCALE_INDEX) + 1);
    mInterfaceConditionName = mParameters["subscale_boundary_condition"].GetString();

    if (mEchoLevel > 1) KRATOS_WATCH(mParameters);

    Check();

    mStepDataSize = mrCoarseModelPart.GetNodalSolutionStepDataSize();

    // Initialize the coarse model part
    InitializeCoarseModelPart();

    // Get the model part hierarchy
    StringVectorType sub_model_parts_names;
    sub_model_parts_names = mrCoarseModelPart.GetSubModelPartNames();
    // sub_model_parts_names = RecursiveGetSubModelPartNames(mrCoarseModelPart);

    // Initialize the refined model part
    InitializeRefinedModelPart(sub_model_parts_names);

    // Copy all the entities to the visualization model part
    InitializeVisualizationModelPart(sub_model_parts_names);
}


void MultiScaleRefiningProcess::Check()
{
    KRATOS_TRY

    KRATOS_CHECK(KratosComponents<Condition>::Has(mInterfaceConditionName));

    KRATOS_CHECK_NOT_EQUAL(mDivisionsAtSubscale, 0);

    KRATOS_CATCH("")
}


void MultiScaleRefiningProcess::ExecuteRefinement()
{
    // Initialize the maps
    IndexIndexMapType node_tag, elem_tag, cond_tag;
    SubModelPartsListUtility model_part_collection(mrCoarseModelPart);
    model_part_collection.ComputeSubModelPartsList(node_tag, cond_tag, elem_tag, mCollections);

    // Get the Id's
    IndexType node_id;
    IndexType elem_id;
    IndexType cond_id;
    GetLastId(node_id, elem_id, cond_id);

    // Clone the nodes and set the nodal flags
    CloneNodesToRefine(node_id);

    // Set the elements and conditions flags
    MarkElementsFromNodalFlag();
    MarkConditionsFromNodalFlag();

    // Create the auxiliary entities
    CreateElementsToRefine(elem_id, elem_tag);
    CreateConditionsToRefine(cond_id, cond_tag);

    // Check and prepare the interface
    IdentifyCurrentInterface();

    // Execute the refinement
    int divisions = mrRefinedModelPart.GetValue(SUBSCALE_INDEX) * mDivisionsAtSubscale;
    mUniformRefinement.SetCustomIds(node_id, elem_id, cond_id);
    mUniformRefinement.Refine(divisions);
    mUniformRefinement.GetLastCreatedIds(node_id, elem_id, cond_id);

    // Update the refined interface after creation
    UpdateRefinedInterface();

    // Update the visualization model part
    UpdateVisualizationAfterRefinement();

    // Reset the flags
    FinalizeRefinement();
}


void MultiScaleRefiningProcess::ExecuteCoarsening()
{
    IdentifyParentNodesToErase();
    IdentifyElementsToErase();
    IdentifyConditionsToErase();
    IdentifyRefinedNodesToErase();

    mUniformRefinement.RemoveRefinedEntities(TO_ERASE);

    // Update the visualization model part
    IdentifyCurrentInterface();
    UpdateVisualizationAfterCoarsening();

    FinalizeCoarsening();
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


void MultiScaleRefiningProcess::InitializeNewModelPart(ModelPart& rReferenceModelPart, ModelPart& rNewModelPart)
{
    // Copy all the tables and properties
    AddAllTablesToModelPart(rReferenceModelPart, rNewModelPart);
    AddAllPropertiesToModelPart(rReferenceModelPart, rNewModelPart);

    // Get the model part hierarchy
    StringVectorType sub_model_parts_names;
    sub_model_parts_names = rReferenceModelPart.GetSubModelPartNames();
    // sub_model_parts_names = RecursiveGetSubModelPartNames(mrCoarseModelPart);

    // Copy the hierarchy to the refined model part
    for (auto name : sub_model_parts_names)
    {
        ModelPart& sub_model_part = rNewModelPart.CreateSubModelPart(name);

        ModelPart& origin_model_part = rReferenceModelPart.GetSubModelPart(name);

        // Copy all the tables and properties
        AddAllTablesToModelPart(origin_model_part, sub_model_part);
        AddAllPropertiesToModelPart(origin_model_part, sub_model_part);
    }
}


void MultiScaleRefiningProcess::InitializeCoarseModelPart()
{
    // Create a model part to store the interface boundary conditions
    if (mrCoarseModelPart.HasSubModelPart(mRefinedInterfaceName))
    {
        mrCoarseModelPart.RemoveNodesFromAllLevels();
        mrCoarseModelPart.RemoveElementsFromAllLevels();
        mrCoarseModelPart.RemoveConditionsFromAllLevels();
    }
    else
        mrCoarseModelPart.CreateSubModelPart(mRefinedInterfaceName);
}


void MultiScaleRefiningProcess::InitializeRefinedModelPart(const StringVectorType& rNames)
{
    // Increase the refinement level
    int subscale_index = mrCoarseModelPart.GetValue(SUBSCALE_INDEX);
    mrRefinedModelPart.SetValue(SUBSCALE_INDEX, ++subscale_index);

    // Copy the variables
    AddVariablesToRefinedModelPart();

    // Create a model part to store the interface boundary conditions
    if (mrRefinedModelPart.HasSubModelPart(mRefinedInterfaceName))
    {
        mrRefinedModelPart.RemoveNodesFromAllLevels();
        mrRefinedModelPart.RemoveElementsFromAllLevels();
        mrRefinedModelPart.RemoveConditionsFromAllLevels();
    }
    else
        mrRefinedModelPart.CreateSubModelPart(mRefinedInterfaceName);
}


void MultiScaleRefiningProcess::InitializeVisualizationModelPart(const StringVectorType& rNames)
{
    // Create a model part to store the interface boundary conditions
    mrVisualizationModelPart.CreateSubModelPart(mRefinedInterfaceName);

    // Add the entities to the root model part
    FastTransferBetweenModelPartsProcess(mrVisualizationModelPart, mrCoarseModelPart)();

    // Add the entities to the submodel parts
    for (auto name : rNames)
    {
        ModelPart& destination = mrVisualizationModelPart.GetSubModelPart(name);
        ModelPart& origin = mrCoarseModelPart.GetSubModelPart(name);
        FastTransferBetweenModelPartsProcess(destination, origin)();
    }        
}


void MultiScaleRefiningProcess::UpdateVisualizationAfterRefinement()
{
    // Remove the refined elements and conditions to substitute them by the refined ones
    mrVisualizationModelPart.RemoveElementsFromAllLevels(MeshingFlags::REFINED);
    mrVisualizationModelPart.RemoveConditionsFromAllLevels(MeshingFlags::REFINED);

    // Remove the refined nodes which are not interface
    ModelPart::NodesContainerType::iterator nodes_begin = mrCoarseModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrCoarseModelPart.Nodes().size()); i++)
    {
        auto coarse_node = nodes_begin + i;
        if (coarse_node->Is((MeshingFlags::REFINED)) && (coarse_node->IsNot(INTERFACE)))
            coarse_node->Set(INSIDE, true);
        else
            coarse_node->Set(INSIDE, false);
    }
    mrVisualizationModelPart.RemoveNodesFromAllLevels(INSIDE);

    // Add the new entities which are refined
    FastTransferBetweenModelPartsProcess(mrVisualizationModelPart, mrRefinedModelPart,
        FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, NEW_ENTITY)();
}


void MultiScaleRefiningProcess::UpdateVisualizationAfterCoarsening()
{
    // Remove the coarsened entities
    mrVisualizationModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    mrVisualizationModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    mrVisualizationModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

    // Add the origin entities which are coarsened
    FastTransferBetweenModelPartsProcess(mrVisualizationModelPart, mrCoarseModelPart,
        FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, MeshingFlags::TO_COARSEN)();
    FastTransferBetweenModelPartsProcess(mrVisualizationModelPart, mrCoarseModelPart,
        FastTransferBetweenModelPartsProcess::EntityTransfered::NODES, INTERFACE)();
}


void MultiScaleRefiningProcess::AddVariablesToRefinedModelPart()
{
    auto variables_list = mrCoarseModelPart.GetNodalSolutionStepVariablesList();

    for (auto variable : variables_list)
    {
        // mrRefinedModelPart.AddNodalSolutionStepVariable(variable);
    }
}


void MultiScaleRefiningProcess::AddAllPropertiesToModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)
{
    const IndexType nprop = rOriginModelPart.NumberOfProperties();
    ModelPart::PropertiesContainerType::iterator prop_begin = rOriginModelPart.PropertiesBegin();

    for (IndexType i = 0; i < nprop; i++)
    {
        auto prop = prop_begin + i;
        rDestinationModelPart.AddProperties(*prop.base());
    }
}


void MultiScaleRefiningProcess::AddAllTablesToModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)
{
    const IndexType ntables = rOriginModelPart.NumberOfTables();
    ModelPart::TablesContainerType::iterator table_begin = rOriginModelPart.TablesBegin();

    for (IndexType i = 0; i < ntables; i++)
    {
        auto table = table_begin + i;
        rDestinationModelPart.AddTable(table.base()->first, table.base()->second);
    }
}


void MultiScaleRefiningProcess::MarkElementsFromNodalFlag()
{
    const int nelems = static_cast<int>(mrCoarseModelPart.Elements().size());
    ModelPart::ElementsContainerType::iterator elem_begin = mrCoarseModelPart.ElementsBegin();

    // We assume all the elements have the same number of nodes
    const IndexType number_of_nodes = elem_begin->GetGeometry().size();

    // We will refine the elements which:
    // 1. all the nodes are to refine
    // 2. at least, one node is marked as new entity
    // NEW_ENTITY flag is used to avoid duplication of refined entities on previous steps
    #pragma omp parallel for
    for (int i = 0; i < nelems; i++)
    {
        auto elem = elem_begin + i;
        bool to_refine = true;
        bool new_entity = false;
        for (IndexType node = 0; node < number_of_nodes; node++)
        {
            if (elem->GetGeometry()[node].IsNot(TO_REFINE))
                to_refine = false;
            
            if (elem->GetGeometry()[node].Is(NEW_ENTITY))
                new_entity = true;
        }
        elem->Set(TO_REFINE, (to_refine && new_entity));
    }
}


void MultiScaleRefiningProcess::MarkConditionsFromNodalFlag()
{
    const int nconds = static_cast<int>(mrCoarseModelPart.Conditions().size());
    ModelPart::ConditionsContainerType::iterator cond_begin = mrCoarseModelPart.ConditionsBegin();

    // We assume all the conditions have the same number of nodes
    const IndexType number_of_nodes = cond_begin->GetGeometry().size();

    // We will refine the conditions which:
    // 1. all the nodes are to refine
    // 2. at least, one node is marked as new entity
    // NEW_ENTITY flag is used to avoid duplication of refined entities on previous steps
    #pragma omp parallel for
    for (int i = 0; i < nconds; i++)
    {
        auto cond = cond_begin + i;
        bool to_refine = true;
        bool new_entity = false;
        for (IndexType node = 0; node < number_of_nodes; node++)
        {
            if (cond->GetGeometry()[node].IsNot(TO_REFINE))
                to_refine = false;
            
            if (cond->GetGeometry()[node].Is(NEW_ENTITY))
                new_entity = true;
        }
        cond->Set(TO_REFINE, (to_refine && new_entity));
    }
}


void MultiScaleRefiningProcess::CloneNodesToRefine(IndexType& rNodeId)
{
    const int nnodes = static_cast<int>(mrCoarseModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator nodes_begin = mrCoarseModelPart.NodesBegin();

    // Adding the nodes to the refined model part
    for (int i = 0; i < nnodes; i++)
    {
        auto coarse_node = nodes_begin + i;
        if (coarse_node->Is(TO_REFINE))
        {
            auto search = mCoarseToRefinedNodesMap.find(coarse_node->Id());
            if (search == mCoarseToRefinedNodesMap.end())
            {
                NodeType::Pointer new_node = mrRefinedModelPart.CreateNewNode(++rNodeId, *coarse_node);
                mCoarseToRefinedNodesMap[coarse_node->Id()] = new_node;
                mRefinedToCoarseNodesMap[rNodeId] = *coarse_node.base();
                new_node->Set(NEW_ENTITY, true);
                new_node->GetValue(FATHER_NODES).resize(0);
                new_node->GetValue(FATHER_NODES).push_back( NodeType::WeakPointer(*coarse_node.base()) );
                coarse_node->Set(NEW_ENTITY, true);
                coarse_node->Set(MeshingFlags::REFINED, true);
            }
        }
    }

    // Adding the nodes to the refined sub model parts
    StringVectorType sub_model_part_names = mrCoarseModelPart.GetSubModelPartNames();
    for (auto name : sub_model_part_names)
    {
        ModelPart& coarse_sub_model_part = mrCoarseModelPart.GetSubModelPart(name);
        ModelPart& refined_sub_model_part = mrRefinedModelPart.GetSubModelPart(name);

        const int nnodes = static_cast<int>(coarse_sub_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = coarse_sub_model_part.NodesBegin();

        for (int i = 0; i < nnodes; i++)
        {
            auto coarse_node = nodes_begin + i;
            if (coarse_node->Is(NEW_ENTITY))
                refined_sub_model_part.AddNode(mCoarseToRefinedNodesMap[coarse_node->Id()]);
        }
    }
}


void MultiScaleRefiningProcess::IdentifyParentNodesToErase()
{
    const int nnodes = static_cast<int>(mrCoarseModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator nodes_begin = mrCoarseModelPart.NodesBegin();

    // Identify the nodes to remove
    for (int i = 0; i < nnodes; i++)
    {
        auto coarse_node = nodes_begin + i;
        if (coarse_node->IsNot(TO_REFINE))
        {
            auto search = mCoarseToRefinedNodesMap.find(coarse_node->Id());
            if (search != mCoarseToRefinedNodesMap.end())
            {
                // We need to ensure the refined mesh does not has dependencies
                if (search->second->IsNot(MeshingFlags::REFINED))
                {
                    coarse_node->Set(MeshingFlags::TO_COARSEN, true);
                    coarse_node->Set(MeshingFlags::REFINED, false);
                    mCoarseToRefinedNodesMap.erase(search);
                    mRefinedToCoarseNodesMap.erase(search->second->Id());
                }
            }
        }
    }
}


void MultiScaleRefiningProcess::IdentifyElementsToErase()
{
    // Identify the parent elements to coarse
    const int nelems_coarse = static_cast<int>(mrCoarseModelPart.Elements().size());
    ModelPart::ElementsContainerType::iterator coarse_elem_begin = mrCoarseModelPart.ElementsBegin();

    // The number of nodes of the elements
    const IndexType element_nodes = coarse_elem_begin->GetGeometry().size();

    #pragma omp parallel for
    for (int i = 0; i < nelems_coarse; i++)
    {
        auto coarse_elem = coarse_elem_begin + i;
        bool to_coarse = false;
        for (IndexType inode = 0; inode < element_nodes; inode++)
        {
            if (coarse_elem->GetGeometry()[inode].Is(MeshingFlags::TO_COARSEN))
                to_coarse = true;
        }
        if (to_coarse)
        {
            coarse_elem->Set(MeshingFlags::TO_COARSEN, true);
            coarse_elem->Set(MeshingFlags::REFINED, false);
        }
    }

    // Identify the refined elements to remove
    const int nelems_ref = static_cast<int>(mrRefinedModelPart.Elements().size());
    ModelPart::ElementsContainerType::iterator ref_elem_begin = mrRefinedModelPart.ElementsBegin();

    #pragma omp parallel for
    for (int i = 0; i < nelems_ref; i++)
    {
        auto refined_elem = ref_elem_begin + i;
        if ((refined_elem->GetValue(FATHER_ELEMENT))->Is(MeshingFlags::TO_COARSEN))
            refined_elem->Set(TO_ERASE, true);
    }
}


void MultiScaleRefiningProcess::IdentifyConditionsToErase()
{
    // Identify the parent conditions  to coarse
    const int nconds_coarse = static_cast<int>(mrCoarseModelPart.Conditions().size());
    ModelPart::ConditionsContainerType::iterator coarse_cond_begin = mrCoarseModelPart.ConditionsBegin();

    // The number of nodes of the conditions
    const IndexType condition_nodes = coarse_cond_begin->GetGeometry().size();

    #pragma omp parallel for
    for (int i = 0; i < nconds_coarse; i++)
    {
        auto coarse_cond = coarse_cond_begin + i;
        bool to_coarse = false;
        for (IndexType inode = 0; inode < condition_nodes; inode++)
        {
            if (coarse_cond->GetGeometry()[inode].Is(MeshingFlags::TO_COARSEN))
                to_coarse = true;
        }
        if (to_coarse)
        {
            coarse_cond->Set(MeshingFlags::TO_COARSEN, true);
            coarse_cond->Set(MeshingFlags::REFINED, false);
        }
    }

    // Identify the refined conditions to remove
    const int nconds_ref = static_cast<int>(mrRefinedModelPart.Conditions().size());
    ModelPart::ConditionsContainerType::iterator ref_cond_begin = mrRefinedModelPart.ConditionsBegin();

    #pragma omp parallel for
    for (int i = 0; i < nconds_ref; i++)
    {
        auto refined_cond = ref_cond_begin + i;
        if ((refined_cond->GetValue(FATHER_CONDITION))->Is(MeshingFlags::TO_COARSEN))
            refined_cond->Set(TO_ERASE, true);
    }
}


void MultiScaleRefiningProcess::IdentifyRefinedNodesToErase()
{
    const IndexType nelems = mrRefinedModelPart.Elements().size();
    if (nelems != 0) // just avoiding segfault in case of an empty coarse model part
    {
        ModelPart::NodeIterator nodes_begin = mrRefinedModelPart.NodesBegin();
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrRefinedModelPart.Nodes().size()); i++)
        {
            auto node = nodes_begin + i;
            node->Set(TO_ERASE, true);
        }

        ModelPart::ElementIterator elements_begin = mrRefinedModelPart.ElementsBegin();
        const IndexType element_nodes = elements_begin->GetGeometry().size();

        for (IndexType i = 0; i < nelems; i++)
        {
            auto elem = elements_begin + i;
            if (elem->IsNot(TO_ERASE))
            {
                for (IndexType inode = 0; inode < element_nodes; inode++)
                {
                    (elem->GetGeometry()[inode]).Set(TO_ERASE, false);
                }
            }
        }
    }
}


void MultiScaleRefiningProcess::CreateElementsToRefine(IndexType& rElemId, IndexIndexMapType& rElemTag)
{
    const int nelems = static_cast<int>(mrCoarseModelPart.Elements().size());
    ModelPart::ElementsContainerType::iterator elements_begin = mrCoarseModelPart.ElementsBegin();

    // We assume all the elements have the same number of nodes
    const IndexType number_of_nodes = elements_begin->GetGeometry().size();

    // The map to add the elements to the sub model parts
    IndexVectorMapType tag_elems_map;

    // #pragma omp parallel for
    for (int i = 0; i < nelems; i++)
    {
        auto coarse_elem = elements_begin + i;
        if (coarse_elem->Is(TO_REFINE))
        {
            Geometry<NodeType>::PointsArrayType p_elem_nodes;
            for (IndexType node = 0; node < number_of_nodes; node++)
            {
                IndexType node_id = coarse_elem->GetGeometry()[node].Id();
                p_elem_nodes.push_back(mCoarseToRefinedNodesMap[node_id]);
            }

            Element::Pointer aux_elem = coarse_elem->Clone(++rElemId, p_elem_nodes);
            mrRefinedModelPart.AddElement(aux_elem);
            
            aux_elem->SetValue(FATHER_ELEMENT, *coarse_elem.base());
            aux_elem->Set(NEW_ENTITY, true);

            IndexType tag = rElemTag[coarse_elem->Id()];
            tag_elems_map[tag].push_back(rElemId);

            coarse_elem->Set(MeshingFlags::REFINED, true);
        }
    }

    // Loop the sub model parts and add the new elements to it
    for (auto& collection : mCollections)
    {
        const auto tag = collection.first;
        if (tag != 0)
        {
            for (auto name : collection.second)
            {
                ModelPart& sub_model_part = mrRefinedModelPart.GetSubModelPart(name);
                sub_model_part.AddElements(tag_elems_map[tag]);
            }
        }
    }
}


void MultiScaleRefiningProcess::CreateConditionsToRefine(IndexType& rCondId, IndexIndexMapType& rCondTag)
{
    const int nconds = static_cast<int>(mrCoarseModelPart.Conditions().size());
    ModelPart::ConditionsContainerType::iterator conditions_begin = mrCoarseModelPart.ConditionsBegin();

    // We assume all the conditions have the same number of nodes
    const IndexType number_of_nodes = conditions_begin->GetGeometry().size();

    // The map to add the conditions to the sub model parts
    IndexVectorMapType tag_conds_map;

    // #pragma omp parallel for
    for (int i = 0; i < nconds; i++)
    {
        auto coarse_cond = conditions_begin + i;
        if (coarse_cond->Is(TO_REFINE))
        {
            Geometry<NodeType>::PointsArrayType p_cond_nodes;
            for (IndexType node = 0; node < number_of_nodes; node++)
            {
                IndexType node_id = coarse_cond->GetGeometry()[node].Id();
                p_cond_nodes.push_back(mCoarseToRefinedNodesMap[node_id]);
            }

            Condition::Pointer aux_cond = coarse_cond->Clone(++rCondId, p_cond_nodes);
            mrRefinedModelPart.AddCondition(aux_cond);
            
            aux_cond->SetValue(FATHER_CONDITION, *coarse_cond.base());
            aux_cond->Set(NEW_ENTITY, true);

            IndexType tag = rCondTag[coarse_cond->Id()];
            tag_conds_map[tag].push_back(rCondId);

            coarse_cond->Set(MeshingFlags::REFINED, true);
        }
    }

    // Loop the sub model parts and add the new conditions to it
    for (auto& collection : mCollections)
    {
        const auto tag = collection.first;
        if (tag != 0)
        {
            for (auto name : collection.second)
            {
                ModelPart& sub_model_part = mrRefinedModelPart.GetSubModelPart(name);
                sub_model_part.AddConditions(tag_conds_map[tag]);
            }
        }
    }
}


void MultiScaleRefiningProcess::FinalizeRefinement()
{
    /// Coarse model part
    // Resetting the nodes flags
    ModelPart::NodeIterator coarse_nodes_begin = mrCoarseModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrCoarseModelPart.Nodes().size()); i++)
    {
        auto node = coarse_nodes_begin + i;
        node->Set(NEW_ENTITY, false);
    }

    /// Refined model part
    // Resetting the nodes flags
    ModelPart::NodeIterator refined_nodes_begin = mrRefinedModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrRefinedModelPart.Nodes().size()); i++)
    {
        auto node = refined_nodes_begin + i;
        node->Set(NEW_ENTITY, false);
    }

    // Resetting the elements flags
    ModelPart::ElementIterator refined_elem_begin = mrRefinedModelPart.ElementsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrRefinedModelPart.Elements().size()); i++)
    {
        auto elem = refined_elem_begin + i;
        elem->Set(NEW_ENTITY, false);
    }

    // Resetting the conditions flags
    ModelPart::ConditionIterator refined_cond_begin = mrRefinedModelPart.ConditionsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrRefinedModelPart.Conditions().size()); i++)
    {
        auto cond = refined_cond_begin + i;
        cond->Set(NEW_ENTITY, false);
    }
}

void MultiScaleRefiningProcess::FinalizeCoarsening()
{
    /// Coarse model part
    // Resetting the nodes flags
    ModelPart::NodeIterator nodes_begin = mrCoarseModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrCoarseModelPart.Nodes().size()); i++)
    {
        auto node = nodes_begin + i;
        node->Set(MeshingFlags::TO_COARSEN, false);
    }

    // Resetting the elements flags
    ModelPart::ElementIterator elements_begin = mrCoarseModelPart.ElementsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrCoarseModelPart.Elements().size()); i++)
    {
        auto elem = elements_begin + i;
        elem->Set(MeshingFlags::TO_COARSEN, false);
    }

    // Resetting the conditions flags
    ModelPart::ConditionIterator conditions_begin = mrCoarseModelPart.ConditionsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrCoarseModelPart.Conditions().size()); i++)
    {
        auto cond = conditions_begin + i;
        cond->Set(MeshingFlags::TO_COARSEN, false);
    }
}


void MultiScaleRefiningProcess::IdentifyCurrentInterface()
{
    // Resetting the flag
    int nnodes = static_cast<int>(mrCoarseModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator nodes_begin = mrCoarseModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i < nnodes; i++)
    {
        auto node = nodes_begin + i;
        node->Set(INTERFACE, false);
    }

    // The number of nodes of the elements
    ModelPart::ElementIterator elem_begin = mrCoarseModelPart.ElementsBegin();
    const IndexType element_nodes = elem_begin->GetGeometry().size();

    // Identify the current interface: 
    // Look for the elements which are not to refine and have some nodes to refine
    for (int i = 0; i < static_cast<int>(mrCoarseModelPart.Elements().size()); i++)
    {
        auto elem = elem_begin + i;
        if (elem->IsNot(MeshingFlags::REFINED))
        {
            // set the nodal flags
            for (IndexType node = 0; node < element_nodes; node++)
            {
                if (elem->GetGeometry()[node].Is(MeshingFlags::REFINED))
                    elem->GetGeometry()[node].Set(INTERFACE, true);
            }
        }
    }
}


void MultiScaleRefiningProcess::UpdateRefinedInterface()
{
    mRefinedInterfaceContainer.clear();

    ModelPart::NodeIterator refined_begin = mrRefinedModelPart.NodesBegin();
    for (int i = 0; i < static_cast<int>(mrRefinedModelPart.Nodes().size()); i++)
    {
        auto node = refined_begin + i;
        bool is_refined_interface = true;
        WeakPointerVector<NodeType>& father_nodes = node->GetValue(FATHER_NODES);
        for (auto father_node = father_nodes.begin(); father_node < father_nodes.end(); father_node++)
        {
            if (father_node->IsNot(INTERFACE))
                is_refined_interface = false;
        }
        if (is_refined_interface)
            mRefinedInterfaceContainer.push_back(NodeType::SharedPointer(*node.base()));
    }
}


void MultiScaleRefiningProcess::TransferLastStepToCoarseModelPart()
{
    ModelPart::NodeIterator refined_begin = mrRefinedModelPart.NodesBegin();
    for (int i = 0; i < static_cast<int>(mrRefinedModelPart.Nodes().size()); i++)
    {
        auto refined_node = refined_begin + i;
        if (refined_node->GetValue(FATHER_NODES).size() == 1)
        {
            double* dest_data = refined_node->GetValue(FATHER_NODES)[0].SolutionStepData().Data(0); // Current step only
            const double* src_data = refined_node->SolutionStepData().Data(0);

            for (IndexType j = 0; j < mStepDataSize; j++)
                dest_data[j] = src_data[j];
        }
    }
}


void MultiScaleRefiningProcess::GetLastId(
    IndexType& rNodesId,
    IndexType& rElemsId,
    IndexType& rCondsId)
{
    // Initialize the output
    rNodesId = 0;
    rElemsId = 0;
    rCondsId = 0;

    // Get the absolute root model part
    ModelPart& root_model_part = mrVisualizationModelPart.GetRootModelPart();

    // Get the maximum node id
    const IndexType nnodes = root_model_part.Nodes().size();
    ModelPart::NodesContainerType::iterator nodes_begin = root_model_part.NodesBegin();
    for (IndexType i = 0; i < nnodes; i++)
    {
        auto inode = nodes_begin + i;
        if (rNodesId < inode->Id())
            rNodesId = inode->Id();
    }

    // Get the maximum element id
    const IndexType nelems = root_model_part.Elements().size();
    ModelPart::ElementsContainerType::iterator elements_begin = root_model_part.ElementsBegin();
    for (IndexType i = 0; i < nelems; i++)
    {
        auto elem = elements_begin + i;
        if (rElemsId < elem->Id())
            rElemsId = elem->Id();
    }

    // Get the maximum condition id
    const IndexType nconds = root_model_part.Conditions().size();
    ModelPart::ConditionsContainerType::iterator conditions_begin = root_model_part.ConditionsBegin();
    for (IndexType i = 0; i < nconds; i++)
    {
        auto cond = conditions_begin + i;
        if (rCondsId < cond->Id())
            rCondsId = cond->Id();
    }
}

} // namespace Kratos