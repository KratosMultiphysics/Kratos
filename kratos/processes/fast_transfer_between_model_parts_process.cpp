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
    const Flags Flag
    ) : Process(),
        mrDestinationModelPart(rDestinationModelPart),
        mrOriginModelPart(rOriginModelPart),
        mEntity(Entity),
        mFlag(Flag)
{
    KRATOS_TRY

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

    // In case of not flag defined we transfer all the elements
    if (mFlag == Flags()) {
        const SizeType num_nodes = mrOriginModelPart.Nodes().size();

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS))
            mrDestinationModelPart.AddNodes(mrOriginModelPart.NodesBegin(),mrOriginModelPart.NodesEnd());

        const SizeType num_elements = mrOriginModelPart.Elements().size();

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
            mrDestinationModelPart.AddElements(mrOriginModelPart.ElementsBegin(),mrOriginModelPart.ElementsEnd());

        const SizeType num_conditions = mrOriginModelPart.Conditions().size();

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
            mrDestinationModelPart.AddConditions(mrOriginModelPart.ConditionsBegin(),mrOriginModelPart.ConditionsEnd());
    } else {
        NodesArrayType vector_nodes;
        ElementsArrayType vector_elements;
        ConditionsArrayType vector_conditions;

        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<NodesArrayType> nodes_buffer(num_threads);
        std::vector<ElementsArrayType> elements_buffer(num_threads);
        std::vector<ConditionsArrayType> conditions_buffer(num_threads);

        // Auxiliar sizes
        const int num_nodes = static_cast<int>(mrOriginModelPart.Nodes().size());
        const int num_elements = static_cast<int>(mrOriginModelPart.Elements().size());
        const int num_conditions = static_cast<int>(mrOriginModelPart.Conditions().size());

        #pragma omp parallel
        {
            const int id = OpenMPUtils::ThisThread();

            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
                #pragma omp for
                for(int i = 0; i < num_nodes; ++i) {
                    auto it_node = mrOriginModelPart.NodesBegin() + i;
                    if (it_node->Is(mFlag)) {
                        (nodes_buffer[id]).insert((nodes_buffer[id]).begin(), *(it_node.base()));
                    }
                }
            }

            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
                #pragma omp for
                for(int i = 0; i < num_elements; ++i) {
                    auto it_elem = mrOriginModelPart.ElementsBegin() + i;
                    if (it_elem->Is(mFlag)) {
                        (elements_buffer[id]).insert((elements_buffer[id]).begin(), *(it_elem.base()));
                    }
                }
            }

            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
                #pragma omp for
                for(int i = 0; i < num_conditions; ++i) {
                    auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
                    if (it_cond->Is(mFlag)) {
                        (conditions_buffer[id]).insert((conditions_buffer[id]).begin(), *(it_cond.base()));
                    }
                }
            }

            // We transfer
            #pragma omp single
            {
                if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                    for( auto& node_buffer : nodes_buffer)
                        mrDestinationModelPart.AddNodes(node_buffer.begin(),node_buffer.end());

                if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                    for( auto& element_buffer : elements_buffer)
                        mrDestinationModelPart.AddElements(element_buffer.begin(),element_buffer.end());

                if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                    for( auto& condition_buffer : conditions_buffer)
                        mrDestinationModelPart.AddConditions(condition_buffer.begin(),condition_buffer.end());
            }
        }
    }

    KRATOS_CATCH("");
}

}
