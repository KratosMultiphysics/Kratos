// System includes

// External includes
#include <boost/python.hpp>

// Project includes

// Application includes
#include "custom_utilities/hdf5_sorted_coordinates_process.h"

namespace Kratos
{
namespace Python
{

void AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    class_<HDF5::SortedCoordinatesProcess, HDF5::SortedCoordinatesProcess::Pointer, bases<Process>, boost::noncopyable>(
        "HDF5SortedCoordinatesProcess", init<std::string, std::string>());
}
} // namespace Python.
} // Namespace Kratos
