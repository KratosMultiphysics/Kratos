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
#include "utilities/constrain.h"

namespace Kratos
{
namespace Python
{

using namespace pybind11;

void  AddConstraintToPython(pybind11::module& m)
{

    class_<MasterSlaveConstraint>(m, "MasterSlaveConstraint" )
    .def(init<MasterSlaveConstraint::IndexType const&>())
	.def(init<>())
    .def("Clear", &MasterSlaveConstraint::Clear)
    .def("AddMaster", &MasterSlaveConstraint::AddMaster)
    .def("AddSlave", &MasterSlaveConstraint::AddSlave)
    .def("GetNumberOfMasters", &MasterSlaveConstraint::GetNumberOfMasters)
    ;

/*
	class_<ModelPart, ModelPart::Pointer, DataValueContainer, Flags >(m,"ModelPart")
		.def(init<std::string const&>())
		.def(init<>())
		.def_property("Name", GetModelPartName, SetModelPartName)
		//  .def_property("ProcessInfo", GetProcessInfo, SetProcessInfo)
		.def_property("ProcessInfo", pointer_to_get_process_info, pointer_to_set_process_info)
		.def("CreateSolutionStep", &ModelPart::CreateSolutionStep)
		.def("CloneSolutionStep", &ModelPart::CloneSolutionStep)
		.def("CreateTimeStep", &ModelPart::CreateTimeStep)
		.def("ReduceTimeStep", &ModelPart::ReduceTimeStep)
		.def("CloneTimeStep", pointer_to_clone_time_step_1)
		.def("CloneTimeStep", pointer_to_clone_time_step_2)    */

}

}  // namespace Python.

} // Namespace Kratos

