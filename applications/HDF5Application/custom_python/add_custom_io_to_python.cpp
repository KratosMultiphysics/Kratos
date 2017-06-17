//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/io.h"

// Application includes
#include "custom_io/hdf5_io.h"

namespace Kratos
{

namespace Python
{

void AddCustomIOToPython()
{
    using namespace boost::python;

    class_<HDF5IO, HDF5IO::Pointer, bases<IO>, boost::noncopyable >("HDF5IO", init<std::string, Flags>())
        .def("WriteModelPart",&HDF5IO::WriteModelPart)
    ;

}

} // namespace Python.

} // Namespace Kratos
