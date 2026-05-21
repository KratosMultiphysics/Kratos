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

namespace Kratos::Python {

// CheckNonHistoricalNodeVariableConsistency()
template<class TContainerType, class TVariableType>
void AddNonHistoricalVariableCheck(pybind11::module& m) {
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variable<TVariableType>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
}

template<class TContainerType, typename... args>
void AddNonHistoricalVariableCheckFolder(pybind11::module& m) {
    (AddNonHistoricalVariableCheck<TContainerType, args>(m), ...);
}

template<class TContainerType>
void AddNonHistoricalVariableCheck(pybind11::module& m) {
    AddNonHistoricalVariableCheckFolder<TContainerType, 
        KRATOS_SINGLE_VARIABLE_TYPES,
        KRATOS_BOUNDED_VECTOR_VARIABLE_TYPES
    >(m);
}

// CheckNonHistoricalNodeVariableConsistency()
template<class TContainerType, class TVariableType>
void AddHistoricalVariableCheck(pybind11::module& m) {
    m.def("CheckHistoricalNodeVariable",[](ModelPart& rModelPart, const Kratos::Variable<TVariableType>& rVariable){MpiDebugUtilities::CheckHistoricalVariable(rModelPart, rVariable);});
}

template<class TContainerType, typename... args>
void AddHistoricalVariableCheckFolder(pybind11::module& m) {
    (AddHistoricalVariableCheck<TContainerType, args>(m), ...);
}

template<class TContainerType>
void AddHistoricalVariableCheck(pybind11::module& m) {
    AddHistoricalVariableCheckFolder<TContainerType, 
        KRATOS_SINGLE_VARIABLE_TYPES,
        KRATOS_BOUNDED_VECTOR_VARIABLE_TYPES
    >(m);
}

void AddMPIDebugUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MpiDebugUtilities>(m,"MPIDebugUtilities");

    AddNonHistoricalVariableCheck<ModelPart::NodesContainerType>(m);
    AddNonHistoricalVariableCheck<ModelPart::ElementsContainerType>(m);
    AddNonHistoricalVariableCheck<ModelPart::ConditionsContainerType>(m);

    AddHistoricalVariableCheck<ModelPart::NodesContainerType>(m);
}

} // namespace Kratos::Python

