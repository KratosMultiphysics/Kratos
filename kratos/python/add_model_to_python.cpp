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

namespace Kratos {
namespace Python {

void  AddModelToPython(pybind11::module& m)
{
    namespace py = pybind11;
    py::class_<Model >(m,"Model")
    .def(py::init<>())
    .def("Reset", &Model::Reset)
    .def("CreateModelPart", [&](Model& self, const std::string& Name){return &self.CreateModelPart(Name);}, py::return_value_policy::reference_internal )
    .def("CreateModelPart", [&](Model& self, const std::string& Name, unsigned int BufferSize){return &self.CreateModelPart(Name, BufferSize);}, py::return_value_policy::reference_internal )
    .def("DeleteModelPart", &Model::DeleteModelPart)
    .def("GetModelPart", &Model::GetModelPart, py::return_value_policy::reference_internal)
    .def("HasModelPart", &Model::HasModelPart)
    .def("__getitem__", &Model::GetModelPart, py::return_value_policy::reference_internal)
    .def("__str__", PrintObject<Model>)
    ;
}

}  // namespace Python.
} // Namespace Kratos

