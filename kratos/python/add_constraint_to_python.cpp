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
#include "utilities/constraint.h" 

namespace Kratos
{
namespace Python
{
    using namespace pybind11;
    
void  AddConstraintToPython(pybind11::module& m)
{
    class_<MasterSlaveConstraint>(m, "MasterSlaveConstraint" )
	.def(init<>())
    .def("Clear", &MasterSlaveConstraint::Clear)
    .def("AddMaster", &MasterSlaveConstraint::AddMaster)
    .def("AddSlave", &MasterSlaveConstraint::AddSlave)
    .def("GetNumberOfMasters", &MasterSlaveConstraint::GetNumberOfMasters)
    ;

}

}  // namespace Python.

} // Namespace Kratos

