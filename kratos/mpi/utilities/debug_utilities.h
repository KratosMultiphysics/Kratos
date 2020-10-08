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
#include <sstream>

#include "includes/model_part.h"
#include "includes/parallel_environment.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/retrieve_global_pointers_by_index_functor.h"
#include "utilities/get_value_functor.h"

#define KRATOS_SINGLE_VARIABLE_TYPES bool, int, double, unsigned int
#define KRATOS_BOUNDED_VECTOR_VARIABLE_TYPES array_1d<double,3>, array_1d<double,4>, array_1d<double,6>, array_1d<double,9>
#define KRATOS_UNBOUNDED_VECTOR_VARIABLE_TYPES Vector, Matrix
namespace Kratos {

class MpiDebugUtilities {
public:

    KRATOS_CLASS_POINTER_DEFINITION(MpiDebugUtilities);

    MpiDebugUtilities() {}

    /* ==== Enable this once we move to C++17 ==== */
    // template <class T, class... Ts>
    // struct is_any : std::integral_constant<bool, (std::is_same<T, Ts>::value || ...)> {};

    // template<class TVarType> struct is_kratos_single_variable : is_any<TVarType, KRATOS_SINGLE_VARIABLE_TYPES> {};
    // template<class TVarType> struct is_kratos_bounded_vector_variable : is_any<TVarType, KRATOS_BOUNDED_VECTOR_VARIABLE_TYPES> {};
    // template<class TVarType> struct is_kratos_unbounded_vector_variable : is_any<TVarType, KRATOS_UNBOUNDED_VECTOR_VARIABLE_TYPES> {};

    // template<class TVarType, class TReturnType> using EnalbeIfSingle = typename std::enable_if<is_kratos_single_variable<TVarType>::value, TReturnType>::type;
    // template<class TVarType, class TReturnType> using EnalbeIfBoundedVector = typename std::enable_if<is_kratos_bounded_vector_variable<TVarType>::value, TReturnType>::type;
    // template<class TVarType, class TReturnType> using EnalbeIfUnbounedVector = typename std::enable_if<is_kratos_unbounded_vector_variable<TVarType>::value, TReturnType>::type;

    // static void CheckNodalHistoricalDatabase(ModelPart & rModelPart) { 
    //     // Build the list of indices and the pointer communicator
    //     DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    //     std::vector<int> indices;

    //     auto node_list = rModelPart.Nodes();

    //     for(auto& node : node_list) {
    //         indices.push_back(node.Id());
    //     }

    //     auto gp_map  = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(node_list, indices, r_default_comm );
    //     auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(node_list, indices, r_default_comm );

    //     GlobalPointerCommunicator<Node<3>> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

    //     // Note: I think this can be change for an access function with std::tuple of the variables instead of a loop, but not in C++11/14
    //     for(auto & variableData: rModelPart.GetNodalSolutionStepVariablesList()) {
    //         auto & variable = KratosComponents<VariableComponentType>::Get(variableData.Name());
    //         CheckNodalHistoricalVariable(rModelPart, variable, pointer_comm, gp_map);
    //     }
    // }
    /* ==== Enable this once we move to C++17 ==== */

    // Non historical Variables
    static bool InteralCmpEq(const bool& rVar1, const bool& rVar2) 
    {
        return rVar1 == rVar2;
    }

    static bool InteralCmpEq(const int& rVar1, const int& rVar2) 
    {
        return rVar1 == rVar2;
    }

    static bool InteralCmpEq(const unsigned int& rVar1, const unsigned int& rVar2) 
    {
        return rVar1 == rVar2;
    }

    static bool InteralCmpEq(const double& rVar1, const double& rVar2) 
    {
        return rVar1 == rVar2;
    }

    static bool InteralCmpEq(const array_1d<double,3>& rVar1, const array_1d<double,3>& rVar2) 
    {
        for(int i = 0; i < 3; i++) {
            if(rVar1[i] != rVar2[i]) return false;
        }
        return true;
    }

    static bool InteralCmpEq(const array_1d<double,4>& rVar1, const array_1d<double,4>& rVar2) 
    {
        for(int i = 0; i < 4; i++) {
            if(rVar1[i] != rVar2[i]) return false;
        }
        return true;
    }

    static bool InteralCmpEq(const array_1d<double,6>& rVar1, const array_1d<double,6>& rVar2) 
    {
        for(int i = 0; i < 6; i++) {
            if(rVar1[i] != rVar2[i]) return false;
        }
        return true;
    }

    static bool InteralCmpEq(const array_1d<double,9>& rVar1, const array_1d<double,9>& rVar2) 
    {
        for(int i = 0; i < 9; i++) {
            if(rVar1[i] != rVar2[i]) return false;
        }
        return true;
    }

    /* ==== Enable this once we move to C++17 ==== */
    // This will prevent having to specialize the functions inline and also will make them generic to all types having traits.
    // template<class TVarType>
    // static EnalbeIfSingle<TVarType, bool> InteralCmpEq(const TVarType& rVar1, const TVarType& rVar2) 
    // {
    //     return rVar1 == rVar2;
    // }

    // template<class TVarType>
    // static EnalbeIfBoundedVector<TVarType, bool> InteralCmpEq(const TVarType& rVar1, const TVarType& rVar2) 
    // {
    //     for(typename TVarType::array_type::size_type i = 0; i < std::tuple_size<typename TVarType::array_type>::value; i++) {
    //         if(rVar1[i] != rVar2[i]) return false;
    //     }
    //     return true;
    // }
    /* ==== Enable this once we move to C++17 ==== */

    template<class TVarType, class TContainerType>
    static void CheckNonHistoricalVariable(
        ModelPart & rModelPart,
        const TContainerType & rContainer,
        const TVarType & rVariable) 
    {
        // Build the list of indices and the pointer communicator
        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
        std::vector<int> indices;

        for(auto& entity : rContainer) {
            indices.push_back(entity.Id());
        }

        auto gp_map  = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(rContainer, indices, r_default_comm );
        auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(rContainer, indices, r_default_comm );

        GlobalPointerCommunicator<typename TContainerType::data_type> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

        CheckNonHistoricalVariable(rModelPart, rContainer, rVariable, pointer_comm, gp_map);
    }

    template<class TContainerType, class TVarType>
    static void CheckNonHistoricalVariable(
        ModelPart & rModelPart,
        const TContainerType & rContainer,
        const TVarType & rVariable,
        GlobalPointerCommunicator<typename TContainerType::data_type> & rPointerCommunicator, 
        std::unordered_map<int, GlobalPointer<typename TContainerType::data_type>> & gp_map) 
    {
        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
        
        bool val_error_detected = false;

        std::stringstream error_stream;

        // Create the data functior
        auto data_proxy = rPointerCommunicator.Apply(
            [rVariable](GlobalPointer<typename TContainerType::data_type>& gp)-> typename TVarType::Type {
                return gp->GetValue(rVariable);
            }
        );

        // Check variable for all entities.
        for(auto& entity : rContainer) {
            auto& gp = gp_map[entity.Id()];

            if(!InteralCmpEq(entity.GetValue(rVariable),data_proxy.Get(gp))) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable value for Id: " << entity.Id() << " Expected: " << entity.GetValue(rVariable) << " Obtained " << data_proxy.Get(gp) << std::endl;
                val_error_detected = true;
            }
        }

        if(val_error_detected) {
            error_stream << "Value error(s) found" << std::endl;
        }

        if(error_stream.rdbuf()->in_avail())
        {
            KRATOS_ERROR << error_stream.str() << std::endl;
        }
    }

    // Historical Variables

    template<class TVarType>
    static void CheckHistoricalVariable(
        ModelPart & rModelPart,
        const TVarType & rVariable) 
    {
        // Build the list of indices and the pointer communicator
        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
        std::vector<int> indices;

        auto node_list = rModelPart.Nodes();

        for(auto& node : node_list) {
            indices.push_back(node.Id());
        }

        auto gp_map  = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(node_list, indices, r_default_comm );
        auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(node_list, indices, r_default_comm );

        GlobalPointerCommunicator<Node<3>> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

        CheckHistoricalVariable(rModelPart, rVariable, pointer_comm, gp_map);
    }

    template<class TVarType>
    static void CheckHistoricalVariable(
        ModelPart & rModelPart,
        const TVarType & rVariable,
        GlobalPointerCommunicator<Node<3>> & rPointerCommunicator, 
        std::unordered_map<int, GlobalPointer<Node<3>>> & gp_map) 
    {
        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
        
        bool val_error_detected = false;
        bool fix_error_detected = false;

        std::stringstream error_stream;

        // Create the data functior
        auto data_proxy = rPointerCommunicator.Apply(
            [rVariable](GlobalPointer< Node<3> >& gp)-> std::pair<typename TVarType::Type, bool> {
                return {gp->FastGetSolutionStepValue(rVariable),gp->IsFixed(rVariable)};
            }
        );

        // Check variable for all nodes.
        for(auto& node : rModelPart.Nodes()) {
            auto& gp = gp_map[node.Id()];

            // Check Variable
            if(!InteralCmpEq(node.FastGetSolutionStepValue(rVariable),data_proxy.Get(gp).first)) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable value for Id: " << node.Id() << " Expected: " << node.FastGetSolutionStepValue(rVariable) << " Obtained " << data_proxy.Get(gp).first << std::endl;
                val_error_detected = true;
            }
        }

        // Check fixity for all nodes.
        for(auto& node : rModelPart.Nodes()) {
            auto& gp = gp_map[node.Id()];

            // Check Fixity
            if(node.IsFixed(rVariable) != data_proxy.Get(gp).second) {
                std::cout << r_default_comm.Rank() << " Inconsistent variable Fix for Id: " << node.Id() << " Expected: " << node.IsFixed(rVariable) << " Obtained " << data_proxy.Get(gp).second << std::endl;
                fix_error_detected = true;
            }
        }

        if(val_error_detected) {
            error_stream << "Value error(s) found" << std::endl;
        }

        if(fix_error_detected) {
            error_stream << "Fixity error(s) found" << std::endl;
        }

        if(error_stream.rdbuf()->in_avail())
        {
            KRATOS_ERROR << error_stream.str() << std::endl;
        }
    }
};
    
} // namespace Kratos