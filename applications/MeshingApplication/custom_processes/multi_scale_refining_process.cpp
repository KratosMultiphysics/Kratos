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
        "element_name"                 : "Element2D3N",
        "condition_name"               : "Condition2D3N",      
        "echo_level"                   : 0,
        "number_of_divisions_at_level" : 2,
        "refining_interface_model_part": "refining_interface",      
        "refining_boundary_condition"  : "Condition2D3N"
    }
    )");

    mParameters.ValidateAndAssignDefaults(DefaultParameters);

    mEchoLevel = mParameters["echo_level"].GetInt();
    mDivisionsAtLevel = mParameters["number_of_divisions_at_level"].GetInt();

    mOwnName = mParameters["own_model_part_name"].GetString();
    mRefinedName = mParameters["refined_model_part_name"].GetString();

    mElementName = mParameters["element_name"].GetString();
    mConditionName = mParameters["condition_name"].GetString();

    std::string own_name = mParameters["own_model_part_name"].GetString();
    std::string refined_name = mParameters["refined_model_part_name"].GetString();

    mInterfaceName = mParameters["refining_interface_model_part"].GetString();
    mInterfaceConditionName = mParameters["refining_boundary_condition"].GetString();

    // Get the model part hierarchy
    StringVectorType sub_model_parts_names;
    if (mrRootModelPart.HasSubModelPart(own_name))
        sub_model_parts_names = mrRootModelPart.GetSubModelPart(own_name).GetSubModelPartNames();
        // sub_model_parts_names = RecursiveGetSubModelPartNames(mrRootModelPart.GetSubModelPart(own_name));
    else
        sub_model_parts_names = mrRootModelPart.GetSubModelPartNames();
        // sub_model_parts_names = RecursiveGetSubModelPartNames(mrRootModelPart);

    Check();

    // Clone the model part at the own level
    InitializeOwnModelPart(sub_model_parts_names);

    // Initialize the refined model part
    InitializeRefinedModelPart(sub_model_parts_names);
}


void MultiScaleRefiningProcess::Check()
{
    KRATOS_TRY

    KRATOS_CHECK(KratosComponents<Element>::Has(mElementName));
    KRATOS_CHECK(KratosComponents<Condition>::Has(mConditionName));

    KRATOS_CHECK_NOT_EQUAL(mDivisionsAtLevel, 0);

    KRATOS_CATCH("")
}


void MultiScaleRefiningProcess::ExecuteRefinement()
{
    ModelPart& refined_model_part = *mpRefinedModelPart.get();

    IndexType node_id;
    IndexType elem_id;
    IndexType cond_id;

    GetLastId(node_id, elem_id, cond_id);

    CloneNodesToRefine(node_id);

    MarkElementsFromNodalFlag();
    MarkElementsFromNodalFlag();

    CreateElementsToRefine(elem_id);
    CreateConditionsToRefine(elem_id);
    
    int divisions = refined_model_part.GetValue(REFINEMENT_LEVEL) * mDivisionsAtLevel;
    auto uniform_refining = UniformRefineUtility<2>(refined_model_part, divisions);

    uniform_refining.Refine(node_id, elem_id, cond_id);
}


void MultiScaleRefiningProcess::ExecuteCoarsening()
{}


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

        // Create a model part to store the interface boundary conditions
        mpOwnModelPart->CreateSubModelPart(mInterfaceName);
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

    // Increase the refinement level
    mpRefinedModelPart->SetValue(REFINEMENT_LEVEL, mpOwnModelPart->GetValue(REFINEMENT_LEVEL)+1);

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

    // Create a model part to store the interface boundary conditions
    mpRefinedModelPart->CreateSubModelPart(mInterfaceName);
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


void MultiScaleRefiningProcess::MarkElementsFromNodalFlag()
{
    const int nelems = static_cast<int>(mpOwnModelPart->Elements().size());
    ModelPart::ElementsContainerType::iterator elem_begin = mpOwnModelPart->ElementsBegin();

    // We assume all the elements have the same number of nodes
    const IndexType number_of_nodes = elem_begin->GetGeometry().size();

    // We will refine the elements which all the nodes are to refine
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
    const int nconds = static_cast<int>(mpOwnModelPart->Conditions().size());
    ModelPart::ConditionsContainerType::iterator cond_begin = mpOwnModelPart->ConditionsBegin();

    const IndexType number_of_nodes = cond_begin->GetGeometry().size();

    // We will refine the conditions which all the nodes are to refine
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
    ModelPart& coarse_model_part = *mpOwnModelPart.get();
    ModelPart& refined_model_part = *mpRefinedModelPart.get();

    const int nnodes = static_cast<int>(coarse_model_part.Nodes().size());
    ModelPart::NodesContainerType::iterator nodes_begin = coarse_model_part.NodesBegin();

    // Adding the nodes to the refined model part
    for (int i = 0; i < nnodes; i++)
    {
        auto coarse_node = nodes_begin + i;
        auto search = mCoarseToRefinedNodesMap.find(coarse_node->Id());
        if ((coarse_node->Is(TO_REFINE)) && (search == mCoarseToRefinedNodesMap.end()))
        {
            NodeType::Pointer new_node = refined_model_part.CreateNewNode(++rNodeId, *coarse_node);
            mCoarseToRefinedNodesMap[coarse_node->Id()] = new_node;
            mRefinedToCoarseNodesMap[rNodeId] = *coarse_node.base();
            coarse_node->Set(NEW_ENTITY, true);
        }
        else if (!(coarse_node->Is(TO_REFINE)) && (search != mCoarseToRefinedNodesMap.end()))
        {
            coarse_node->Set(TO_ERASE, false);
            mCoarseToRefinedNodesMap.erase(search);
            mRefinedToCoarseNodesMap.erase(search->second->Id());
            /**
             * What happens if I need to remove a node?
             * and how do I remove all the divided nodes?
             * I should implement a TransferDataToRefinedNodes
             * Remove from the maps
             * Remove from the maps inside the utility
             * I think that the utility should remove the entities
             */
        }
    }

    // Adding the nodes to the refined sub model parts
    StringVectorType sub_model_part_names = coarse_model_part.GetSubModelPartNames();
    for (auto name : sub_model_part_names)
    {
        ModelPart& coarse_sub_model_part = coarse_model_part.GetSubModelPart(name);
        ModelPart& refined_sub_model_part = refined_model_part.GetSubModelPart(name);

        const int nnodes = static_cast<int>(coarse_sub_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = coarse_sub_model_part.NodesBegin();

        for (int i = 0; i < nnodes; i++)
        {
            auto coarse_node = nodes_begin + i;
            if (coarse_node->Is(NEW_ENTITY))
                refined_sub_model_part.AddNode(mCoarseToRefinedNodesMap[coarse_node->Id()]);
        }
    }

    // // Resetting the flag
    // #pragma omp parallel for
    // for (int i = 0; i < nnodes; i++)
    // {
    //     auto node = nodes_begin + i;
    //     node->Set(NEW_ENTITY, false);
    // }
}


void MultiScaleRefiningProcess::CreateElementsToRefine(IndexType& rElemId)
{
    ModelPart& coarse_model_part = *mpOwnModelPart.get();
    ModelPart& refined_model_part = *mpRefinedModelPart.get();

    const int nelems = static_cast<int>(coarse_model_part.Elements().size());
    ModelPart::ElementsContainerType::iterator elements_begin = coarse_model_part.ElementsBegin();

    // We assume all the elements have the same number of nodes
    const IndexType number_of_nodes = elements_begin->GetGeometry().size();

    #pragma omp parallel for
    for (int i = 0; i < nelems; i++)
    {
        auto elem = elements_begin + i;
        if (elem->Is(TO_REFINE))
        {
            Geometry<NodeType>::PointsArrayType p_elem_nodes;
            for (IndexType node = 0; node < number_of_nodes; node++)
            {
                IndexType node_id = elem->GetGeometry()[node].Id();
                p_elem_nodes.push_back(mCoarseToRefinedNodesMap[node_id]);
            }

            Element::Pointer aux_elem = refined_model_part.CreateNewElement(
                mElementName,
                ++rElemId,
                p_elem_nodes,
                elem->pGetProperties());
            
            aux_elem->SetValue(FATHER_ELEMENT, *elem.base());
        }
    }

    // Loop the sub model parts and add the new elements to it
}


void MultiScaleRefiningProcess::CreateConditionsToRefine(IndexType& rCondId)
{}


void MultiScaleRefiningProcess::IdentifyRefiningInterface()
{
    // 0. Reset the flags
    
    // 1. Identify the nodes which define the boundary
    const int nelems = static_cast<int>(mpOwnModelPart->Elements().size());
    ModelPart::ElementsContainerType::iterator elem_begin = mpOwnModelPart->ElementsBegin();

    // The number of nodes of the elements
    const IndexType element_nodes = elem_begin->GetGeometry().size();

    // Look for the elements which are not to refine and have some nodes to refine
    for (int i = 0; i < nelems; i++)
    {
        auto elem = elem_begin + i;
        if (elem->IsNot(TO_REFINE))
        {
            for (IndexType node = 0; node < element_nodes; node++)
            {
                if (elem->GetGeometry()[node].Is(TO_REFINE))
                    elem->GetGeometry()[node].Set(INTERFACE, true);
            }
        }
    }

    // TODO: here I need to tranfer the INTERFACE flag from the father nodes to the middle nodes
    /* do some stuff */

    // 2. Remove the old conditions
    ModelPart& interface = mpRefinedModelPart->GetSubModelPart(mInterfaceName);
    const int nconds = static_cast<int>(interface.Conditions().size());
    ModelPart::ConditionsContainerType::iterator cond_begin = interface.ConditionsBegin();

    // The number of nodes of the conditions
    const IndexType condition_nodes = cond_begin->GetGeometry().size();

    // Find the conditions which are not interface
    #pragma omp parallel for
    for (int i = 0; i < nconds; i++)
    {
        auto cond = cond_begin + i;
        bool to_erase = false;
        for (IndexType node = 0; node < condition_nodes; node++)
        {
            if (cond->GetGeometry()[node].IsNot(INTERFACE))
                to_erase = true;
        }
        cond->Set(TO_ERASE, to_erase);

        // We need to preserve the old interface
        // Warning: what should I do with OMP???????
        // NOTE: I am accessing the nodes twice
        if (!to_erase)
        {
            for (IndexType node = 0; node < condition_nodes; node++)
                cond->GetGeometry()[node].Set(OLD_ENTITY, true);
        }
    }

    // 3. And finally, create the new conditions where needed
    // The condition should inherit the refining level from the coarse element
    /**
     *  loop elements
     *      loop edges
     *          loop nodes
     *              if (all nodes are INTERFACE and are not OLD_ENTITY)
     *                  Create condition
     **/
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
    ModelPart& root_model_part = mpOwnModelPart->GetRootModelPart();

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

    // // Get the number of threads
    // const int num_threads = OpenMPUtils::GetNumThreads();

    // // Get the maximum node id
    // const int nnodes = static_cast<int>(root_model_part.Nodes().size());
    // ModelPart::NodesContainerType::iterator nodes_begin = root_model_part.NodesBegin();
    // std::vector<IndexType> nodes_id(num_threads);
    // #pragma omp parallel
    // {
    //     const int thread_id = OpenMPUtils::ThisThread();

    //     #pragma omp for
    //     for (int i = 0; i < nnodes; i++)
    //     {
    //         auto node = nodes_begin + i;
    //         if (nodes_id[thread_id] < node->Id())
    //             nodes_id[thread_id] = node->Id();
    //     }

    //     #pragma omp single
    //     {
    //         rNodesId = *std::max_element(nodes_id.begin(), nodes_id.end());
    //     }
    // }

    // // Get the maximum element id
    // const int nelems = static_cast<int>(root_model_part.Nodes().size());
    // ModelPart::ElementsContainerType::iterator elements_begin = root_model_part.ElementsBegin();
    // std::vector<IndexType> elems_id(num_threads);
    // #pragma omp parallel
    // {
    //     const int thread_id = OpenMPUtils::ThisThread();

    //     #pragma omp for
    //     for (int i = 0; i < nelems; i++)
    //     {   /* WARNING: HERE IS SOME SEGMENTATION FAULT */
    //         auto elem = elements_begin + i;
    //         if (elems_id[thread_id] < elem->Id())
    //             elems_id[thread_id] = elem->Id();
    //     }

    //     #pragma omp single
    //     {
    //         rElemsId = *(std::max_element(elems_id.begin(), elems_id.end()));
    //     }
    // }
}


} // namespace Kratos