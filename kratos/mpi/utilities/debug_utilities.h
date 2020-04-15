//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Carlos A. Roig
//                   Riccardo Rossi
//

#pragma once

#include <vector>
#include "includes/model_part.h"
#include "includes/parallel_environment.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/retrieve_global_pointers_by_index_functor.h"
#include "utilities/get_value_functor.h"

namespace Kratos {

class MpiDebugUtilities {
public:

    KRATOS_CLASS_POINTER_DEFINITION(MpiDebugUtilities);

    MpiDebugUtilities() {}

    static void CheckNodalHistoricalDatabaseConsistency(ModelPart::NodesContainerType & rNodes) { 
        // Build the list of indices and the pointer communicator
        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
        std::vector<int> indices;

        for(auto& node : rNodes) {
            indices.push_back(node.Id());
        }

        auto gp_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(rNodes, indices, r_default_comm );
        auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(rNodes, indices, r_default_comm );

        GlobalPointerCommunicator< Node<3>> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

        // If there are nodes, asume they have the same variable list and get it.
        if(rNodes.size()) {
            auto & variable_list = rNodes.begin()->pGetVariablesList();

            // Note: I think this can be change for an access function with std::tuple of the variables instead of a loop, but not in C++11/14
            for(auto & variable: variable_list) {
                CheckNodalHistoricalVariableConsistency(rNodes, variable, pointer_comm, gp_map);
            }
        }
    }

    template<class TVarType>
    static void CheckNodalHistoricalVariableConsistency(
        ModelPart::NodesContainerType & rNodes,
        const TVarType & rVariable,
        GlobalPointerCommunicator< Node<3>> & rPointerCommunicator, 
        std::unordered_map<int, GlobalPointer<Node<3>>>) 
    {
        // Create the data functior
        auto data_proxy = pointer_comm.Apply(
            [rVariable](GlobalPointer< Node<3> >& gp)->std::pair<typename TVarType::Type, bool> {
                return {gp->FastGetSolutionStepValue(rVariable),gp->IsFixed(rVariable)};
            }
        );

        // Check all nodes.
        for(auto& node : rNodes) {
            auto& gp = gp_map[node.Id()];

            // Check Variable
            if(data_proxy.Get(gp).first != node.FastGetSolutionStepValue(rVariable)) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable Val for Id: " << node.Id() <<  " Expected: " << node.FastGetSolutionStepValue(rVariable) << " Obtained " << data_proxy.Get(gp).first << std::endl;
            }

            // Check Fixity
            if(data_proxy.Get(gp).second != node.IsFixed(rVariable)) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable Fix for Id: " << node.Id() <<  " Expected: " << node.IsFixed(rVariable) << " Obtained " << data_proxy.Get(gp).second << std::endl;
            }
        }
    }

    template<class TVarType>
    static void CheckNodalDatabaseConsistency(ModelPart::NodesContainerType & rNodes, const TVarType & rVariable) {
        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
        std::vector<int> indices;

        for(auto& node : rNodes) {
            indices.push_back(node.Id());
        }

        auto gp_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(rNodes, indices, r_default_comm );
        auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(rNodes, indices, r_default_comm );

        GlobalPointerCommunicator< Node<3>> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

        auto data_proxy = pointer_comm.Apply(
            [rVariable](GlobalPointer< Node<3> >& gp)->std::pair<typename TVarType::Type, bool> {
                return {gp->FastGetSolutionStepValue(rVariable),gp->IsFixed(rVariable)};
            }
        );

        for(auto& node : rNodes) {
            auto& gp = gp_map[node.Id()];

            // Check Variable
            if(data_proxy.Get(gp).first != node.FastGetSolutionStepValue(rVariable)) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable Val for Id: " << node.Id() <<  " Expected: " << node.FastGetSolutionStepValue(rVariable) << " Obtained " << data_proxy.Get(gp).first << std::endl;
            }

            // Check Fixity
            if(data_proxy.Get(gp).second != node.IsFixed(rVariable)) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable Fix for Id: " << node.Id() <<  " Expected: " << node.IsFixed(rVariable) << " Obtained " << data_proxy.Get(gp).second << std::endl;
            }
        }
    }
};
    
} // namespace Kratos