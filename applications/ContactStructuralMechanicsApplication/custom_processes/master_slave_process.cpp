// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
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
    
    ModelPart& r_contact_model_part = mrThisModelPart.GetSubModelPart("Contact");
    
    // Now we iterate over the nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    const int num_nodes = static_cast<int>(r_nodes_array.size());
    
    // Now we iterate over the conditions
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    const int num_conditions = static_cast<int>(r_conditions_array.size());
    
    std::vector<IndexType> index_node, index_cond;
    
    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        std::vector<IndexType> index_node_buffer, index_cond_buffer;
        
        #pragma omp for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;
            if (it_node->IsDefined(INTERFACE)) {
                if (it_node->Is(INTERFACE))
                    (index_node_buffer).push_back(it_node->Id());
            }
        }
        
        #pragma omp for
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = it_cond_begin + i;
            const auto& r_geometry = it_cond->GetGeometry();
            
            // We set the condition as master or slave (master by default)
            bool is_interface = true;
            bool is_slave = true;
            for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
                const auto& r_node = r_geometry[i_node];
                if (r_node.IsDefined(INTERFACE)) {
                    if (r_node.IsNot(INTERFACE)) is_interface = false;
                    else if (r_node.IsNot(SLAVE)) is_slave = false;
                }
            }
            if (is_interface) {
                index_cond_buffer.push_back(it_cond->Id());
                if (is_slave) it_cond->Set(SLAVE, true);
                else it_cond->Set(MASTER, true);
            }
        }
        
        // Combine buffers together
        #pragma omp critical
        {
            std::move(index_node_buffer.begin(),index_node_buffer.end(),back_inserter(index_node));
            std::move(index_cond_buffer.begin(),index_cond_buffer.end(),back_inserter(index_cond));
        }
    }
    
    // Adding nodes and conditions
    r_contact_model_part.AddNodes(index_node);
    r_contact_model_part.AddConditions(index_cond);

    KRATOS_CATCH("");
} // class MasterSlaveProcess
} // namespace Kratos
