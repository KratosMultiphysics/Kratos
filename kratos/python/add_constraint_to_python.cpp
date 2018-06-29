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
using namespace pybind11;

template <class TContainerType, class TVariableType>
bool HasHelperFunction(TContainerType &el, const TVariableType &rVar)
{
    return el.Has(rVar);
}

template <class TContainerType, class TVariableType>
void SetValueHelperFunction(TContainerType &el, const TVariableType &rVar, const typename TVariableType::Type &Data)
{
    el.SetValue(rVar, Data);
}

template <class TContainerType, class TVariableType>
typename TVariableType::Type GetValueHelperFunction(TContainerType &el, const TVariableType &rVar)
{
    return el.GetValue(rVar);
}

void AddConstraintToPython(pybind11::module &m)
{
    class_<MasterSlaveConstraint, MasterSlaveConstraint::Pointer>(m, "MasterSlaveConstraint")
        .def(init<>())
        .def(init<int>())
        .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
        .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
        .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< bool > >)
        .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
        .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
        ;
}

} // namespace Python.

} // Namespace Kratos
