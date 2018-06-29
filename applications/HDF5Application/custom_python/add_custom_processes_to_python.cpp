// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes

// Application includes
#include "custom_processes/hdf5_xdmf_connectivities_writer_process.h"

namespace Kratos
{
namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<HDF5::XdmfConnectivitiesWriterProcess, HDF5::XdmfConnectivitiesWriterProcess::Pointer, Process>(m,"HDF5XdmfConnectivitiesWriterProcess")
    .def(init<const std::string&, const std::string&>())
    ;
}
} // namespace Python.
} // Namespace Kratos
