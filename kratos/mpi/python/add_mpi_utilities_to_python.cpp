//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "add_mpi_utilities_to_python.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/synchronize_constraints_utility.h"
namespace Kratos {
namespace Python {

// Master slave constraints

void AddConstraint1(MpiConstraintsUtility& rConstraintsUtility,
                    std::string ConstraintName,
                    const ModelPart::IndexType Id,
                    const ModelPart::IndexType MasterNodeId,
                    ModelPart::DoubleVariableType& rMasterVariable,
                    const ModelPart::IndexType SlaveNodeId,
                    ModelPart::DoubleVariableType& rSlaveVariable,
                    double Weight,
                    double Constant)
{
    rConstraintsUtility.AddConstraint(
        ConstraintName, Id, MasterNodeId, rMasterVariable, SlaveNodeId,
        rSlaveVariable, Weight, Constant);
}

void AddConstraint2(MpiConstraintsUtility& rConstraintsUtility,
                    std::string ConstraintName,
                    const ModelPart::IndexType Id,
                    const ModelPart::IndexType MasterNodeId,
                    ModelPart::VariableComponentType& rMasterVariable,
                    const ModelPart::IndexType SlaveNodeId,
                    ModelPart::VariableComponentType& rSlaveVariable,
                    double Weight,
                    double Constant)
{
    rConstraintsUtility.AddConstraint(
        ConstraintName, Id, MasterNodeId, rMasterVariable, SlaveNodeId,
        rSlaveVariable, Weight, Constant);
}

void AddMPIUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ModelPartCommunicatorUtilities>(m,
                                               "ModelPartCommunicatorUtilities")
        .def_static("SetMPICommunicator", &ModelPartCommunicatorUtilities::SetMPICommunicator);

    py::class_<ParallelFillCommunicator>(m, "ParallelFillCommunicator")
        .def(py::init<ModelPart&>())
        .def("Execute", &ParallelFillCommunicator::Execute)
        .def("PrintDebugInfo", &ParallelFillCommunicator::PrintDebugInfo);

    py::class_<MpiConstraintsUtility>(m, "MpiConstraintsUtility")
        .def(py::init<ModelPart&>())
        .def("AddConstraint", AddConstraint1)
        .def("AddConstraint", AddConstraint2)
        .def("SynchronizeAndCreateConstraints",
             &MpiConstraintsUtility::SynchronizeAndCreateConstraints);
}

} // namespace Python
} // namespace Kratos
