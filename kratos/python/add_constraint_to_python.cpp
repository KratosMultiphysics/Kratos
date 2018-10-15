//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_constraint_to_python.h"
#include "includes/master_slave_constraint.h"

namespace Kratos
{
namespace Python
{

template <class TContainerType, class TVariableType>
bool HasHelperFunctionConstraint(TContainerType &el, const TVariableType &rVar)
{
    return el.Has(rVar);
}

template <class TContainerType, class TVariableType>
void SetValueHelperFunctionConstraint(TContainerType &el, const TVariableType &rVar, const typename TVariableType::Type &Data)
{
    el.SetValue(rVar, Data);
}

template <class TContainerType, class TVariableType>
typename TVariableType::Type GetValueHelperFunctionConstraint(TContainerType &el, const TVariableType &rVar)
{
    return el.GetValue(rVar);
}

void AddConstraintToPython(pybind11::module &m)
{
    pybind11::class_<MasterSlaveConstraint, MasterSlaveConstraint::Pointer, MasterSlaveConstraint::BaseType>(m, "MasterSlaveConstraint")
        .def(pybind11::init<>())
        .def(pybind11::init<int>())
        .def("__setitem__", SetValueHelperFunctionConstraint< MasterSlaveConstraint, Variable< bool > >)
        .def("__getitem__", GetValueHelperFunctionConstraint< MasterSlaveConstraint, Variable< bool > >)
        .def("Has", HasHelperFunctionConstraint< MasterSlaveConstraint, Variable< bool > >)
        .def("SetValue", SetValueHelperFunctionConstraint< MasterSlaveConstraint, Variable< bool > >)
        .def("GetValue", GetValueHelperFunctionConstraint< MasterSlaveConstraint, Variable< bool > >)
        .def("__str__", KRATOS_DEF_PYTHON_STR(MasterSlaveConstraint))
        ;
}

} // namespace Python.

} // Namespace Kratos
