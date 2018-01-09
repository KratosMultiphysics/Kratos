//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/table.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef Table<double,double> DoubleTableType;

double TableGetNearestValue(DoubleTableType& ThisTable, double X)
{
    return ThisTable.GetNearestValue(X);
}


void  AddTableToPython()
{
    class_<DoubleTableType, DoubleTableType::Pointer>("PiecewiseLinearTable")
    .def(init<Matrix const&>())
//    .def(init<Variable<double> const&, Variable<double> const&>())
    .def("GetValue", &DoubleTableType::GetValue)
    .def("GetDerivative",&DoubleTableType::GetDerivative)
    .def("GetNearestValue", TableGetNearestValue)
    .def("AddRow", &DoubleTableType::PushBack)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

