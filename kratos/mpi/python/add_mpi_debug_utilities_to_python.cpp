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
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "mpi/utilities/debug_utilities.h"

namespace Kratos {
namespace Python {

// template<class TVariableTpe>
// void CheckScalarVariables(pybind11::module& m) {

// }

/* ==== Enable this once we move to C++17 ==== */
// template<class TContainerType, class TVariableType>
// void AddNonHistoricalVariableCheck(pybind11::module& m) {
//     m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variable<TVariableType>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
// }

// // CheckNonHistoricalNodeVariableConsistency()
// template<class TContainerType, typename... args>
// void AddNonHistoricalVariableCheckFolder(pybind11::module& m) {
//     (AddNonHistoricalVariableCheck<TContainerType, args>(m), ...);
// }

// template<class TContainerType>
// void AddNonHistoricalVariableCheck(pybind11::module& m) {
//     AddNonHistoricalVariableCheckFolder<TContainerType, 
//         KRATOS_SINGLE_VARIABLE_TYPES,
//         KRATOS_BOUNDED_VECTOR_VARIABLE_TYPES
//     >(m);
// }

// template<class TContainerType, class TVariableType>
// void AddHistoricalVariableCheck(pybind11::module& m) {
//     m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variable<TVariableType>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
// }

// // CheckNonHistoricalNodeVariableConsistency()
// template<class TContainerType, typename... args>
// void AddHistoricalVariableCheck(pybind11::module& m) {
//     (AddHistoricalVariableCheck<TContainerType, args>(m), ...);
// }
/* ==== Enable this once we move to C++17 ==== */

void AddMPIDebugUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MpiDebugUtilities>(m,"MPIDebugUtilities");

    /* ==== Enable this once we move to C++17 ==== */
    // AddNonHistoricalVariableCheck<ModelPart::NodesContainerType>(m);
    // AddNonHistoricalVariableCheck<ModelPart::ElementsContainerType>(m);
    // AddNonHistoricalVariableCheck<ModelPart::ConditionsContainerType>(m);

    // AddHistoricalVariableCheck<ModelPart::NodesContainerType, KRATOS_SINGLE_VARIABLE_TYPES, KRATOS_BOUNDED_VECTOR_VARIABLE_TYPES>(m);

    /* ==== Enable this once we move to C++17 ==== */

    // .def("SetMPICommunicator",&MpiDebugUtilities::SetMPICommunicator)
    // ;

    /* ==== Disable this once we move to C++17 ==== */
    // Non Historical
    // Nodes
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<int>& rVariable){});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<unsigned int>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<double>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<bool>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<array_1d<double, 3>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<array_1d<double, 4>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<array_1d<double, 6>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::NodesContainerType& rContainer, const Kratos::Variable<array_1d<double, 9>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});

    // Elements
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<int>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<unsigned int>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<double>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<bool>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<array_1d<double, 3>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<array_1d<double, 4>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<array_1d<double, 6>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ElementsContainerType& rContainer, const Kratos::Variable<array_1d<double, 9>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});

    // Conditions
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<int>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<unsigned int>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<double>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<bool>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<array_1d<double, 3>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<array_1d<double, 4>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<array_1d<double, 6>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const ModelPart::ConditionsContainerType& rContainer, const Kratos::Variable<array_1d<double, 9>>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});

    // Historical
    // Nodes
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<int>& rVariable){});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<unsigned int>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<double>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<bool>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<array_1d<double, 3>>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<array_1d<double, 4>>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<array_1d<double, 6>>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<array_1d<double, 9>>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
    /* ==== Disable this once we move to C++17 ==== */
}

} // namespace Python
} // namespace Kratos

