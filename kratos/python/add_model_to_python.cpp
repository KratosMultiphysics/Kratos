//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pooyan Dadvand
//

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "python/add_model_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  AddModelToPython()
{
    class_<Model, Model::Pointer, boost::noncopyable >("Model", init<>())
    .def("AddModelPart", &Model::AddModelPart)
    .def("GetModelPart", &Model::GetModelPart, return_internal_reference<>())
//     .def("__setitem__", &Model::AddModelPart)
    .def("__getitem__", &Model::GetModelPart, return_internal_reference<>())
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

