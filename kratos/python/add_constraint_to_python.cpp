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
    pybind11::class_<MasterSlaveConstraint, MasterSlaveConstraint::Pointer, MasterSlaveConstraint::BaseType>(m, "MasterSlaveConstraint")
    .def(pybind11::init<>())
    .def(pybind11::init<int>())

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 3>  > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 4>  > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 4>  > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 4>  > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 4>  > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 4>  > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 6>  > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 6>  > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 6>  > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 6>  > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 6>  > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 9>  > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 9>  > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 9>  > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 9>  > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< array_1d<double, 9>  > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< Vector > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< DenseVector<int> > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< DenseVector<int> > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< DenseVector<int> > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< DenseVector<int> > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< DenseVector<int> > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< Matrix > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< int > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< int > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< int > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< double > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< double > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< double > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< bool > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< bool > >)

    .def("__setitem__", SetValueHelperFunction< MasterSlaveConstraint, Variable< std::string > >)
    .def("__getitem__", GetValueHelperFunction< MasterSlaveConstraint, Variable< std::string > >)
    .def("Has", HasHelperFunction< MasterSlaveConstraint, Variable< std::string > >)
    .def("SetValue", SetValueHelperFunction< MasterSlaveConstraint, Variable< std::string > >)
    .def("GetValue", GetValueHelperFunction< MasterSlaveConstraint, Variable< std::string > >)

    .def("__str__", PrintObject<MasterSlaveConstraint>)
    ;
}

} // namespace Python.

} // Namespace Kratos
