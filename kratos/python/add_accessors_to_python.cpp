//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "add_accessors_to_python.h"
#include "includes/define_python.h"
#include "includes/accessor.h"
#include "includes/properties.h"

namespace Kratos::Python
{

template<class TDataType, class TAccessorType, class... TAccessorArgs>
auto& InternalCreateAndSetAccessorToProperties(
    Properties& rProperties,
    const Variable<TDataType>& rVariable,
    TAccessorArgs... rAccessorArgs)
{
    auto p_new = Kratos::make_unique<TAccessorType>(std::forward<TAccessorArgs>(rAccessorArgs)...);
    rProperties.SetAccessor(rVariable, std::move(p_new));
    return rProperties.GetAccessor(rVariable);
}

void AddAccessorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Accessor, Accessor::UniquePointer>(m, "Accessor")
        .def("Create", [](Properties& rProp, Variable<double>& rVariable) -> auto& 
        { return InternalCreateAndSetAccessorToProperties<double, Accessor>(rProp,rVariable); },
        py::return_value_policy::reference_internal)
        ;
}

}  // namespace Kratos::Python.

