//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "custom_python/add_custom_io_to_python.h"
#include "custom_io/gid_io_eigen.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomIOToPython()
{
        class_<GidIOEigen, GidIOEigen::Pointer, bases<GidIO<>>, boost::noncopyable>(
            "GidIOEigen",init<std::string const&, 
                              GiD_PostMode,
                              MultiFileFlag,
                              WriteDeformedMeshFlag,
                              WriteConditionsFlag>())
        .def("WriteEigenResults",&GidIOEigen::WriteEigenResults)
        ;
    
}

}  // namespace Python.

} // Namespace Kratos