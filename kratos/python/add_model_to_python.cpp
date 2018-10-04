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
#include "includes/kernel.h"
#include "containers/model.h"
#include "python/add_model_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  AddModelToPython(pybind11::module& m)
{
    class_<Model >(m,"Model")
    .def(init<>())
    .def("Reset", &Model::Reset)
    .def("CreateModelPart", [&](Model& self, const std::string& Name){return &self.CreateModelPart(Name);}, return_value_policy::reference_internal )
    .def("CreateModelPart", [&](Model& self, const std::string& Name, unsigned int BufferSize){return &self.CreateModelPart(Name, BufferSize);}, return_value_policy::reference_internal )
    .def("DeleteModelPart", &Model::DeleteModelPart)
    .def("GetModelPart", &Model::GetModelPart, return_value_policy::reference_internal)
    .def("HasModelPart", &Model::HasModelPart)
    .def("__getitem__", &Model::GetModelPart, return_value_policy::reference_internal)
    .def("__repr__", [](const Model& self) -> const std::string { std::stringstream ss;  ss << self; return ss.str(); })
    ;
}

}  // namespace Python.

} // Namespace Kratos

