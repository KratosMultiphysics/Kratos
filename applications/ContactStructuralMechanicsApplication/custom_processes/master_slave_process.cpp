// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/master_slave_process.h"

namespace Kratos
{
    void MasterSlaveProcess::Execute()
    {
        KRATOS_TRY;
        
        ModelPart& contact_model_part = mrThisModelPart.GetSubModelPart("Contact");
        
        // Now we iterate over the nodes
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        // Now we iterate over the conditions
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<std::size_t> index_node;
        std::vector<std::vector<std::size_t>> index_node_buffer(num_threads);
        std::vector<std::size_t> index_cond;
        std::vector<std::vector<std::size_t>> index_cond_buffer(num_threads);
        
        #pragma omp parallel
        {
            const int thread_id = OpenMPUtils::ThisThread();
            
            #pragma omp for
            for(int i = 0; i < num_nodes; ++i) 
            {
                auto it_node = nodes_array.begin() + i;
                if (it_node->IsDefined(INTERFACE) == true)
                {
                    if (it_node->Is(INTERFACE) == true)
                    {
                        (index_node_buffer[thread_id]).push_back(it_node->Id());
                    }
                }
            }
            
            #pragma omp for
            for(int i = 0; i < num_conditions; ++i) 
            {
                auto it_cond = conditions_array.begin() + i;
                const auto& this_geometry = it_cond->GetGeometry();
                
                // We set the condition as master or slave (master by default)
                bool is_interface = true;
                bool is_slave = true;
                for (unsigned int i_node = 0; i_node < this_geometry.size(); ++i_node)
                {
                    auto& this_node = this_geometry[i_node];
                    if (this_node.IsDefined(INTERFACE) == true)
                    {
                        if (!this_node.Is(INTERFACE)) is_interface = false;
                        else if (!this_node.Is(SLAVE)) is_slave = false;
                    }
                }
                if (is_interface == true)
                {
                    (index_cond_buffer[thread_id]).push_back(it_cond->Id());
                    if (is_slave == true)  it_cond->Set(SLAVE, true);
                    else it_cond->Set(MASTER, true);
                }
            }
            
            // Combine buffers together
            #pragma omp single
            {
                for( auto& node_buffer : index_node_buffer)
                {
                    std::move(node_buffer.begin(),node_buffer.end(),back_inserter(index_node));
                }
                for( auto& cond_buffer : index_cond_buffer)
                {
                    std::move(cond_buffer.begin(),cond_buffer.end(),back_inserter(index_cond));
                }
            }
        }
        
        // Adding nodes and conditions
        contact_model_part.AddNodes(index_node);
        contact_model_part.AddConditions(index_cond);

        KRATOS_CATCH("");
    }
}
