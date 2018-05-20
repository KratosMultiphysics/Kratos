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
    m.def("Model", &Kernel::GetModel, return_value_policy::reference);
    
    //NOTE: we call this class "ModelInterface" instead of "Model" since the cosntructor is emulated as a standalone function which gets it from the kernel
    class_<Model >(m,"ModelInterface")
    .def("Reset", &Model::Reset)
    .def("AddModelPart", &Model::AddModelPart)
    .def("GetModelPart", &Model::GetModelPart, return_value_policy::reference_internal)
    .def("HasModelPart", &Model::HasModelPart)
    .def("__getitem__", &Model::GetModelPart, return_value_policy::reference_internal)
    .def("__repr__", [](const Model& self) -> const std::string { std::stringstream ss;  ss << self; return ss.str(); })
    ;
}

}  // namespace Python.

} // Namespace Kratos

