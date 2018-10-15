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

// Project includes
#include "includes/define_python.h"
#include "containers/model.h"
#include "python/add_model_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  AddModelToPython(pybind11::module& m)
{
    class_<Model, Model::Pointer >(m,"Model")
    .def(init<>())
    .def("AddModelPart", &Model::AddModelPart)
    .def("GetModelPart", &Model::GetModelPart, return_value_policy::reference_internal)
    .def("HasModelPart", &Model::HasModelPart)
    .def("__getitem__", &Model::GetModelPart, return_value_policy::reference_internal)
    .def("__str__", KRATOS_DEF_PYTHON_STR(Model))
    ;
}

}  // namespace Python.

} // Namespace Kratos

