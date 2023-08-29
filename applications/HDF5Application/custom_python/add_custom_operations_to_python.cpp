//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya
//

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes

// Application includes
#include "custom_operations/hdf5_xdmf_connectivities_writer_operation.h"

// Include base h
#include "add_custom_operations_to_python.h"

namespace Kratos {
namespace Python {

void AddCustomOperationsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<HDF5::XdmfConnectivitiesWriterOperation, HDF5::XdmfConnectivitiesWriterOperation::Pointer, Operation>(m,"HDF5XdmfConnectivitiesWriterOperation")
        .def(py::init<const std::string&>(), py::arg("file_name"))
        ;
}

} // namespace Python.
} // Namespace Kratos
