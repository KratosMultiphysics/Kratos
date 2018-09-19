//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "processes/fast_transfer_between_model_parts_process.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

FastTransferBetweenModelPartsProcess::FastTransferBetweenModelPartsProcess(
    ModelPart& rDestinationModelPart,
    ModelPart& rOriginModelPart,
    const EntityTransfered Entity,
    const Flags Flag,
    const bool ReplicateEntities
    ) : Process(),
        mrDestinationModelPart(rDestinationModelPart),
        mrOriginModelPart(rOriginModelPart),
        mEntity(Entity),
        mFlag(Flag)
{
    KRATOS_TRY

    // If the entities are replicated or transfered
    if (ReplicateEntities) {
        this->Set(MODIFIED, true);
    } else {
        this->Set(MODIFIED, false);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::operator()()
{
    Execute();
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::Execute()
{
    KRATOS_TRY;

    if (this->IsNot(MODIFIED)) {
        // In case of not flag defined we transfer all the elements
        if (mFlag == Flags()) {
            TransferWithoutFlags();
        } else {
            TransferWithFlags();
        }
    } else {
        // In case of not flag defined we transfer all the elements
        if (mFlag == Flags()) {
            ReplicateWithoutFlags();
        } else {
            ReplicateWithFlags();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::TransferWithoutFlags()
{
    const SizeType num_nodes = mrOriginModelPart.Nodes().size();

    if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS))
        mrDestinationModelPart.AddNodes(mrOriginModelPart.NodesBegin(),mrOriginModelPart.NodesEnd());

    const SizeType num_elements = mrOriginModelPart.Elements().size();

    if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
        mrDestinationModelPart.AddElements(mrOriginModelPart.ElementsBegin(),mrOriginModelPart.ElementsEnd());

    const SizeType num_conditions = mrOriginModelPart.Conditions().size();

    if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
        mrDestinationModelPart.AddConditions(mrOriginModelPart.ConditionsBegin(),mrOriginModelPart.ConditionsEnd());

    const SizeType num_constraints = mrOriginModelPart.MasterSlaveConstraints().size();

    if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
        mrDestinationModelPart.AddMasterSlaveConstraints(mrOriginModelPart.MasterSlaveConstraintsBegin(),mrOriginModelPart.MasterSlaveConstraintsEnd());
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::TransferWithFlags()
{
    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<NodesArrayType> nodes_buffer(num_threads);
    std::vector<ElementsArrayType> elements_buffer(num_threads);
    std::vector<ConditionsArrayType> conditions_buffer(num_threads);
    std::vector<MasterSlaveConstraintArrayType> constraints_buffer(num_threads);

    // Auxiliar sizes
    const int num_nodes = static_cast<int>(mrOriginModelPart.Nodes().size());
    const int num_elements = static_cast<int>(mrOriginModelPart.Elements().size());
    const int num_conditions = static_cast<int>(mrOriginModelPart.Conditions().size());
    const int num_constraints = static_cast<int>(mrOriginModelPart.Conditions().size());

    #pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            #pragma omp for
            for(int i = 0; i < num_nodes; ++i) {
                auto it_node = mrOriginModelPart.NodesBegin() + i;
                if (it_node->Is(mFlag)) {
                    (nodes_buffer[thread_id]).insert((nodes_buffer[thread_id]).begin(), *(it_node.base()));
                }
            }
        }

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
            #pragma omp for
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = mrOriginModelPart.ElementsBegin() + i;
                if (it_elem->Is(mFlag)) {
                    (elements_buffer[thread_id]).insert((elements_buffer[thread_id]).begin(), *(it_elem.base()));
                }
            }
        }

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            #pragma omp for
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
                if (it_cond->Is(mFlag)) {
                    (conditions_buffer[thread_id]).insert((conditions_buffer[thread_id]).begin(), *(it_cond.base()));
                }
            }
        }

        if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS)) {
            #pragma omp for
            for(int i = 0; i < num_constraints; ++i) {
                auto it_const = mrOriginModelPart.MasterSlaveConstraintsBegin() + i;
                if (it_const->Is(mFlag)) {
                    (constraints_buffer[thread_id]).insert((constraints_buffer[thread_id]).begin(), *(it_const.base()));
                }
            }
        }

        // We transfer
        #pragma omp single
        {
            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                for( auto& node_buffer : nodes_buffer)
                    mrDestinationModelPart.AddNodes(node_buffer.begin(),node_buffer.end());
        }

        #pragma omp single
        {
            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                for( auto& element_buffer : elements_buffer)
                    mrDestinationModelPart.AddElements(element_buffer.begin(),element_buffer.end());
        }

        #pragma omp single
        {
            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                for( auto& condition_buffer : conditions_buffer)
                    mrDestinationModelPart.AddConditions(condition_buffer.begin(),condition_buffer.end());
        }

        #pragma omp single
        {
            if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
                for( auto& constraint_buffer : constraints_buffer)
                    mrDestinationModelPart.AddMasterSlaveConstraints(constraint_buffer.begin(),constraint_buffer.end());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::ReorderAllIds(ModelPart& rThisModelPart)
{
    NodesArrayType& nodes_array = rThisModelPart.Nodes();
    for(SizeType i = 0; i < nodes_array.size(); ++i)
        (nodes_array.begin() + i)->SetId(i + 1);

    ElementsArrayType& element_array = rThisModelPart.Elements();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->SetId(i + 1);

    ConditionsArrayType& condition_array = rThisModelPart.Conditions();
    for(SizeType i = 0; i < condition_array.size(); ++i)
        (condition_array.begin() + i)->SetId(i + 1);

    MasterSlaveConstraintArrayType& constraint_array = rThisModelPart.MasterSlaveConstraints();
    for(SizeType i = 0; i < constraint_array.size(); ++i)
        (constraint_array.begin() + i)->SetId(i + 1);
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::ReplicateWithoutFlags()
{
    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<NodesArrayType> nodes_buffer(num_threads);
    std::vector<ElementsArrayType> elements_buffer(num_threads);
    std::vector<ConditionsArrayType> conditions_buffer(num_threads);
    std::vector<MasterSlaveConstraintArrayType> constraints_buffer(num_threads);

    // Reordering ids (necessary for consistency)
    ModelPart& root_model_part = mrOriginModelPart.GetRootModelPart();
    ReorderAllIds(root_model_part);

    // Getting the auxiliar values
    const SizeType total_num_nodes = root_model_part.Nodes().size();
    const int num_nodes = static_cast<int>(mrOriginModelPart.Nodes().size());
    const SizeType total_num_elements = root_model_part.Elements().size();
    const int num_elements = static_cast<int>(mrOriginModelPart.Elements().size());
    const SizeType total_num_conditions = root_model_part.Conditions().size();
    const int num_conditions = static_cast<int>(mrOriginModelPart.Conditions().size());
    const SizeType total_num_constraints = root_model_part.MasterSlaveConstraints().size();
    const int num_constraints = static_cast<int>(mrOriginModelPart.MasterSlaveConstraints().size());

    #pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            #pragma omp for
            for(int i = 0; i < num_nodes; ++i) {
                auto it_node = mrOriginModelPart.NodesBegin() + i;
                NodeType::Pointer p_new_node = it_node->Clone();
                p_new_node->SetId(total_num_nodes + i + 1);
                (nodes_buffer[thread_id]).insert((nodes_buffer[thread_id]).begin(), p_new_node);
            }
        }

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
            #pragma omp for
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = mrOriginModelPart.ElementsBegin() + i;
                Element::Pointer p_new_elem = it_elem->Clone(total_num_elements + i + 1, it_elem->GetGeometry());
                (elements_buffer[thread_id]).insert((elements_buffer[thread_id]).begin(), p_new_elem);
            }
        }

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            #pragma omp for
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
                Condition::Pointer p_new_cond = it_cond->Clone(total_num_conditions + i + 1, it_cond->GetGeometry());
                (conditions_buffer[thread_id]).insert((conditions_buffer[thread_id]).begin(), p_new_cond);
            }
        }

        if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS)) {
            #pragma omp for
            for(int i = 0; i < num_constraints; ++i) {
                auto it_const = mrOriginModelPart.MasterSlaveConstraintsBegin() + i;
                MasterSlaveConstraint::Pointer p_new_const = it_const->Clone(total_num_constraints + i + 1);
                (constraints_buffer[thread_id]).insert((constraints_buffer[thread_id]).begin(), p_new_const);
            }
        }

        // We add to the model part
        #pragma omp single
        {
            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                for( auto& node_buffer : nodes_buffer)
                    mrDestinationModelPart.AddNodes(node_buffer.begin(),node_buffer.end());
        }

        #pragma omp single
        {
            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                for( auto& element_buffer : elements_buffer)
                    mrDestinationModelPart.AddElements(element_buffer.begin(),element_buffer.end());
        }

        #pragma omp single
        {
            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                for( auto& condition_buffer : conditions_buffer)
                    mrDestinationModelPart.AddConditions(condition_buffer.begin(),condition_buffer.end());
        }

        #pragma omp single
        {
            if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
                for( auto& constraint_buffer : constraints_buffer)
                    mrDestinationModelPart.AddMasterSlaveConstraints(constraint_buffer.begin(),constraint_buffer.end());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::ReplicateWithFlags()
{
    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<NodesArrayType> nodes_buffer(num_threads);
    std::vector<ElementsArrayType> elements_buffer(num_threads);
    std::vector<ConditionsArrayType> conditions_buffer(num_threads);
    std::vector<MasterSlaveConstraintArrayType> constraints_buffer(num_threads);

    // Reordering ids (necessary for consistency)
    ModelPart& root_model_part = mrOriginModelPart.GetRootModelPart();
    ReorderAllIds(root_model_part);

    // Getting the auxiliar values
    const SizeType total_num_nodes = root_model_part.Nodes().size();
    const int num_nodes = static_cast<int>(mrOriginModelPart.Nodes().size());
    const SizeType total_num_elements = root_model_part.Elements().size();
    const int num_elements = static_cast<int>(mrOriginModelPart.Elements().size());
    const SizeType total_num_conditions = root_model_part.Conditions().size();
    const int num_conditions = static_cast<int>(mrOriginModelPart.Conditions().size());
    const SizeType total_num_constraints = root_model_part.MasterSlaveConstraints().size();
    const int num_constraints = static_cast<int>(mrOriginModelPart.MasterSlaveConstraints().size());

    #pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            #pragma omp for
            for(int i = 0; i < num_nodes; ++i) {
                auto it_node = mrOriginModelPart.NodesBegin() + i;
                if (it_node->Is(mFlag)) {
                    NodeType::Pointer p_new_node = it_node->Clone();
                    p_new_node->SetId(total_num_nodes + i + 1);
                    (nodes_buffer[thread_id]).insert((nodes_buffer[thread_id]).begin(), p_new_node);
                }
            }
        }

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
            #pragma omp for
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = mrOriginModelPart.ElementsBegin() + i;
                if (it_elem->Is(mFlag)) {
                    Element::Pointer p_new_elem = it_elem->Clone(total_num_elements + i + 1, it_elem->GetGeometry());
                    (elements_buffer[thread_id]).insert((elements_buffer[thread_id]).begin(), p_new_elem);
                }
            }
        }

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            #pragma omp for
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
                if (it_cond->Is(mFlag)) {
                    Condition::Pointer p_new_cond = it_cond->Clone(total_num_conditions + i + 1, it_cond->GetGeometry());
                    (conditions_buffer[thread_id]).insert((conditions_buffer[thread_id]).begin(), p_new_cond);
                }
            }
        }

        if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS)) {
            #pragma omp for
            for(int i = 0; i < num_constraints; ++i) {
                auto it_const = mrOriginModelPart.MasterSlaveConstraintsBegin() + i;
                if (it_const->Is(mFlag)) {
                    MasterSlaveConstraint::Pointer p_new_const = it_const->Clone(total_num_constraints + i + 1);
                    (constraints_buffer[thread_id]).insert((constraints_buffer[thread_id]).begin(), p_new_const);
                }
            }
        }

        // We add to the model part
        #pragma omp single
        {
            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                for( auto& node_buffer : nodes_buffer)
                    mrDestinationModelPart.AddNodes(node_buffer.begin(),node_buffer.end());
        }

        #pragma omp single
        {
            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                for( auto& element_buffer : elements_buffer)
                    mrDestinationModelPart.AddElements(element_buffer.begin(),element_buffer.end());
        }

        #pragma omp single
        {
            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                for( auto& condition_buffer : conditions_buffer)
                    mrDestinationModelPart.AddConditions(condition_buffer.begin(),condition_buffer.end());
        }

        #pragma omp single
        {
            if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
                for( auto& constraint_buffer : constraints_buffer)
                    mrDestinationModelPart.AddMasterSlaveConstraints(constraint_buffer.begin(),constraint_buffer.end());
        }
    }
}

}
