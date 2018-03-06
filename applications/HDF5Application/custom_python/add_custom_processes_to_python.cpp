// System includes

// External includes
#include <boost/python.hpp>

// Project includes

// Application includes
#include "custom_processes/hdf5_xdmf_connectivities_writer_process.h"

namespace Kratos
{
namespace Python
{

void AddCustomProcessesToPython()
{
    using namespace boost::python;

    class_<HDF5::XdmfConnectivitiesWriterProcess, bases<Process>, boost::noncopyable>(
        "HDF5XdmfConnectivitiesWriterProcess", init<const std::string&, const std::string&>())
        ;
}
} // namespace Python.
} // Namespace Kratos
