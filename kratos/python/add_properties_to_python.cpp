//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"
#include "python/containers_interface.h"
#include "utilities/interface_folder.h"

namespace Kratos::Python {

namespace py = pybind11;

using MeshType = Mesh<Node, Properties, Element, Condition>;
using IndexType = std::size_t;
using PropertiesContainerType = PointerVectorSet<Properties, IndexedObject>;

void  AddPropertiesToPython(pybind11::module& m)
{
    auto binding = py::class_<Properties, Properties::Pointer, Properties::BaseType >(m,"Properties")
    .def(py::init<Kratos::Properties::IndexType>())
    .def(py::init<const Properties&>())

    .def("IsEmpty", &Properties::IsEmpty)
    .def("HasTables", &Properties::HasTables)
    .def("HasVariables", &Properties::HasVariables)
    .def("AddSubProperties", &Properties::AddSubProperties)
    .def("SetSubProperties", &Properties::SetSubProperties)
    .def("NumberOfSubproperties", &Properties::NumberOfSubproperties)

    .def("HasSubProperties", [](Properties& rProp, IndexType Index){return rProp.HasSubProperties(Index);})
    .def("GetSubProperties", [](Properties& rProp, IndexType Index){return rProp.pGetSubProperties(Index);})
    .def("GetSubProperties", [](Properties& rProp){return rProp.GetSubProperties();})
    .def("GetSubProperties", [](const Properties& rProp){return rProp.GetSubProperties();})

    .def("GetTable", [](Properties& rProp, Variable<double>& rXVar, Variable<double>& rYVar){return rProp.GetTable(rXVar, rYVar);}, py::return_value_policy::reference_internal)
    .def("SetTable", [](Properties& rProp, Variable<double>& rXVar, Variable<double>& rYVar, typename Properties::TableType& rData){rProp.SetTable(rXVar, rYVar, rData);})
    .def("HasTable", [](Properties& rProp, Variable<double>& rXVar, Variable<double>& rYVar){return rProp.HasTable(rXVar, rYVar);})

    .def("__str__", PrintObject<Properties>)
    ;

    AddToInterface<int, double, bool, std::string, Matrix, Vector, array_1d<double,3>, array_1d<double,6>, ConstitutiveLaw::Pointer>(
        binding, []<typename TDataType>(auto binding) {

        using TVariableType = Variable<TDataType>;

        binding
            .def("Has", [](Properties &rProp, TVariableType &rVar){return rProp.Has(rVar);})
            .def("Erase", [](Properties &rProp, TVariableType &rVar){rProp.Erase(rVar);})
            .def("SetValue", [](Properties &rProp, TVariableType &rVar, TDataType &rData){rProp.SetValue(rVar,rData);})
            .def("GetValue", [](Properties &rProp, TVariableType &rVar){return rProp.GetValue(rVar);})
            .def("__setitem__", [](Properties &rProp, TVariableType &rVar, TDataType &rData){rProp.SetValue(rVar,rData);})
            .def("__getitem__", [](Properties &rProp, TVariableType &rVar){return rProp.GetValue(rVar);})
            ;
    });

    PointerVectorSetPythonInterface<MeshType::PropertiesContainerType>().CreateInterface(m,"PropertiesArray");
}

}  // namespace Kratos::Python.
