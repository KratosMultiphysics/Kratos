//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_communicator_to_python.h"
#include "includes/communicator.h"

namespace Kratos::Python
{

template<class TDataType>
bool CommunicatorSynchronizeVariable(Communicator& rCommunicator, Variable<TDataType> const& ThisVariable)
{
    return rCommunicator.SynchronizeVariable(ThisVariable);
}

template<class TDataType>
bool CommunicatorSynchronizeNonHistoricalVariable(Communicator& rCommunicator, Variable<TDataType> const& ThisVariable)
{
    return rCommunicator.SynchronizeNonHistoricalVariable(ThisVariable);
}

template<class TDataType>
bool CommunicatorAssembleCurrentData(Communicator& rCommunicator, Variable<TDataType> const& ThisVariable)
{
    return rCommunicator.AssembleCurrentData(ThisVariable);
}

template<class TDataType>
bool CommunicatorAssembleNonHistoricalData(Communicator& rCommunicator, Variable<TDataType> const& ThisVariable)
{
    return rCommunicator.AssembleNonHistoricalData(ThisVariable);
}

void AddCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<Communicator > (m,"Communicator")
        .def(py::init<>())
        .def("MyPID", &Communicator::MyPID)
        .def("TotalProcesses", &Communicator::TotalProcesses)
        .def("GlobalNumberOfNodes", &Communicator::GlobalNumberOfNodes)
        .def("GlobalNumberOfElements", &Communicator::GlobalNumberOfElements)
        .def("GlobalNumberOfConditions", &Communicator::GlobalNumberOfConditions)
        .def("GlobalNumberOfMasterSlaveConstraints", &Communicator::GlobalNumberOfMasterSlaveConstraints)
        .def("GetNumberOfColors", &Communicator::GetNumberOfColors)
        .def("NeighbourIndices", [](Communicator& rCommunicator){return rCommunicator.NeighbourIndices();}, py::return_value_policy::reference_internal )
        .def("SynchronizeNodalSolutionStepsData", &Communicator::SynchronizeNodalSolutionStepsData)
        .def("SynchronizeNodalFlags", &Communicator::SynchronizeNodalFlags)
        .def("SynchronizeOrNodalFlags", &Communicator::SynchronizeOrNodalFlags)
        .def("SynchronizeAndNodalFlags", &Communicator::SynchronizeAndNodalFlags)
        .def("SynchronizeDofs", &Communicator::SynchronizeDofs)
        .def("LocalMesh", [](Communicator& rCommunicator){return rCommunicator.LocalMesh();}, py::return_value_policy::reference_internal )
        .def("LocalMesh", [](Communicator& rCommunicator, Communicator::IndexType Index){return rCommunicator.LocalMesh(Index);}, py::return_value_policy::reference_internal )
        .def("GhostMesh", [](Communicator& rCommunicator){return rCommunicator.GhostMesh();}, py::return_value_policy::reference_internal )
        .def("GhostMesh", [](Communicator& rCommunicator, Communicator::IndexType Index){return rCommunicator.GhostMesh(Index);}, py::return_value_policy::reference_internal )
        .def("InterfaceMesh", [](Communicator& rCommunicator){return rCommunicator.InterfaceMesh();}, py::return_value_policy::reference_internal )
        .def("InterfaceMesh", [](Communicator& rCommunicator, Communicator::IndexType Index){return rCommunicator.InterfaceMesh(Index);}, py::return_value_policy::reference_internal )
        .def("GetDataCommunicator", &Communicator::GetDataCommunicator, py::return_value_policy::reference_internal )
        .def("SynchronizeVariable", CommunicatorSynchronizeVariable<int> )
        .def("SynchronizeVariable", CommunicatorSynchronizeVariable<double> )
        .def("SynchronizeVariable", CommunicatorSynchronizeVariable<array_1d<double,3> > )
        .def("SynchronizeVariable", CommunicatorSynchronizeVariable<Vector> )
        .def("SynchronizeVariable", CommunicatorSynchronizeVariable<Matrix> )
        .def("SynchronizeNonHistoricalVariable", CommunicatorSynchronizeNonHistoricalVariable<int> )
        .def("SynchronizeNonHistoricalVariable", CommunicatorSynchronizeNonHistoricalVariable<double> )
        .def("SynchronizeNonHistoricalVariable", CommunicatorSynchronizeNonHistoricalVariable<array_1d<double,3> > )
        .def("SynchronizeNonHistoricalVariable", CommunicatorSynchronizeNonHistoricalVariable<Vector> )
        .def("SynchronizeNonHistoricalVariable", CommunicatorSynchronizeNonHistoricalVariable<Matrix> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<int> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<double> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<array_1d<double,3> > )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<Vector> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<Matrix> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<int> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<double> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<array_1d<double,3> > )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<Vector> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<Matrix> )
        .def("__str__", PrintObject<Communicator>);
        ;
}

} // namespace Kratos::Python.
