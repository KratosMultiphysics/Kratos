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
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_constraint_to_python.h"
#include "includes/master_slave_constraint.h"
#include "includes/linear_master_slave_constraint.h"
#include "constraints/slip_constraint.h"

namespace Kratos
{
namespace Python
{

template <class TVariableType>
bool HasHelperFunctionMasterSlaveConstraint(MasterSlaveConstraint &el, const TVariableType &rVar)
{
    return el.Has(rVar);
}

template <class TVariableType>
void SetValueHelperFunctionMasterSlaveConstraint(MasterSlaveConstraint &el, const TVariableType &rVar, const typename TVariableType::Type &Data)
{
    el.SetValue(rVar, Data);
}

template <class TVariableType>
typename TVariableType::Type GetValueHelperFunctionMasterSlaveConstraint(MasterSlaveConstraint &el, const TVariableType &rVar)
{
    return el.GetValue(rVar);
}

void AddConstraintToPython(pybind11::module &m)
{
    pybind11::class_<MasterSlaveConstraint, MasterSlaveConstraint::Pointer, MasterSlaveConstraint::BaseType, Flags>(m, "MasterSlaveConstraint")
    .def(pybind11::init<>())
    .def(pybind11::init<int>())

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 3>  > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 4>  > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 4>  > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 4>  > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 4>  > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 4>  > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 6>  > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 6>  > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 6>  > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 6>  > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 6>  > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 9>  > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 9>  > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 9>  > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 9>  > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< array_1d<double, 9>  > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< Vector > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< Vector > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< Vector > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< DenseVector<int> > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< DenseVector<int> > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< DenseVector<int> > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< DenseVector<int> > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< DenseVector<int> > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< Matrix > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< int > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< int > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< int > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< int > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< int > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< double > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< double > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< double > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< double > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< double > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< bool > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< bool > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< bool > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< bool > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< bool > >)

    .def("__setitem__", SetValueHelperFunctionMasterSlaveConstraint< Variable< std::string > >)
    .def("__getitem__", GetValueHelperFunctionMasterSlaveConstraint< Variable< std::string > >)
    .def("Has", HasHelperFunctionMasterSlaveConstraint< Variable< std::string > >)
    .def("SetValue", SetValueHelperFunctionMasterSlaveConstraint< Variable< std::string > >)
    .def("GetValue", GetValueHelperFunctionMasterSlaveConstraint< Variable< std::string > >)

    .def("GetMasterDofsVector", &MasterSlaveConstraint::GetMasterDofsVector, pybind11::return_value_policy::reference_internal)
    .def("GetSlaveDofsVector", &MasterSlaveConstraint::GetSlaveDofsVector, pybind11::return_value_policy::reference_internal)
    .def("CalculateLocalSystem", &MasterSlaveConstraint::CalculateLocalSystem)

    .def("__str__", PrintObject<MasterSlaveConstraint>)
    ;

    pybind11::class_<LinearMasterSlaveConstraint, LinearMasterSlaveConstraint::Pointer, MasterSlaveConstraint, Flags>
            (m, "LinearMasterSlaveConstraint")
    .def(pybind11::init<
        MasterSlaveConstraint::IndexType,
        MasterSlaveConstraint::DofPointerVectorType&,
        MasterSlaveConstraint::DofPointerVectorType&,
        const MasterSlaveConstraint::MatrixType&,
        const MasterSlaveConstraint::VectorType&
        >())
    ;


    pybind11::class_<SlipConstraint, SlipConstraint::Pointer, LinearMasterSlaveConstraint, Flags>
            (m, "SlipConstraint")
    .def(pybind11::init<
        SlipConstraint::IndexType,
        Dof<double>* ,
        Dof<double>* ,
        array_1d<double,3>
        >())
    .def(pybind11::init<
        SlipConstraint::IndexType,
        Dof<double>* ,
        Dof<double>* ,
        Dof<double>* ,
        array_1d<double,3>
        >())
    ;


}

} // namespace Python.

} // Namespace Kratos
