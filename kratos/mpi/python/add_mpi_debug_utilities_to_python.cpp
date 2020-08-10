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

#define KRATOS_BASIC_VARIABLE_TYPES bool, int, double, unsigned int
#define KRATOS_ARRAY_VARIABLE_TYPES Vector, Matrix
#define KRATOS_COMP_VARIABLE_TYPES array_1d<double,3>, array_1d<double,4>, array_1d<double,6>, array_1d<double,9>
#define KRATOS_COMP_ADAPTOR_VARIABLE_TYPES VectorComponentAdaptor<array_1d<double,3>>, VectorComponentAdaptor<array_1d<double,4>>, VectorComponentAdaptor<array_1d<double,6>>, VectorComponentAdaptor<array_1d<double,9>>

template<class TContainerType, class TVariableType>
void AddNonHistoricalVariableCheck(pybind11::module& m) {
    m.def("CheckNonHistoricalNodeVariable",[](ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variable<TVariableType>& rVariable){MpiDebugUtilities::CheckNonHistoricalVariable(rModelPart, rContainer, rVariable);});
}

// CheckNonHistoricalNodeVariableConsistency()
template<class TContainerType, typename... args>
void AddNonHistoricalVariableCheckFolder(pybind11::module& m) {
    (AddNonHistoricalVariableCheck<TContainerType, args>(m), ...);
}

template<class TContainerType>
void AddNonHistoricalVariableCheck(pybind11::module& m) {
    AddNonHistoricalVariableCheckFolder<TContainerType, 
        KRATOS_BASIC_VARIABLE_TYPES,
        KRATOS_COMP_VARIABLE_TYPES
    >(m);
}

//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<bool> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<int> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<double> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<unsigned int> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<Vector> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<Matrix> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<array_1d<double,3>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<array_1d<double,4>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<array_1d<double,6>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<array_1d<double,9>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<VectorComponentAdaptor<array_1d<double,3>>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<VectorComponentAdaptor<array_1d<double,4>>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<VectorComponentAdaptor<array_1d<double,6>>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
//     m.def("CheckNonHistoricalNodeVariable",[](const ModelPart& rModelPart, const TContainerType& rContainer, const Kratos::Variabel<VectorComponentAdaptor<array_1d<double,9>>> rVariable){CheckNonHistoricalNodeVariable(rModelPart, rContainer, rVariable);});
// }

void AddMPIDebugUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MpiDebugUtilities>(m,"MPIDebugUtilities");

    AddNonHistoricalVariableCheck<ModelPart::NodesContainerType>(m);
    AddNonHistoricalVariableCheck<ModelPart::ElementsContainerType>(m);
    AddNonHistoricalVariableCheck<ModelPart::ConditionsContainerType>(m);

    // CheckScalarVariables<int>(m);
    // CheckScalarVariables<double>(m);
    // CheckScalarVariables<bool>(m);

    // .def("SetMPICommunicator",&MpiDebugUtilities::SetMPICommunicator)
    // ;
}

} // namespace Python
} // namespace Kratos

