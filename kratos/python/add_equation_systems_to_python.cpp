//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
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
#include "python/add_equation_systems_to_python.h"
#include "equation_systems/equation_system.h"
#include "includes/dof.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

namespace Python
{
void  AddEquationSystemsToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef EquationSystem<SpaceType,  LocalSpaceType, Dof<double> > EquationSystemType;

    void (EquationSystemType::*pointer_to_initalize)() = &EquationSystemType::Initialize;
    const EquationSystemType::SystemMatrixType& (EquationSystemType::*pointer_to_get_system_matrix)() const = &EquationSystemType::GetSystemMatrix;
    const EquationSystemType::SystemVectorType& (EquationSystemType::*pointer_to_get_results)() const = &EquationSystemType::GetResults;
    const EquationSystemType::SystemVectorType& (EquationSystemType::*pointer_to_get_rhs)() const = &EquationSystemType::GetRightHandSide;
    const EquationSystemType::SystemMatrixType& (EquationSystemType::*pointer_to_get_dirichlet_matrix)() const = &EquationSystemType::GetDirichletMatrix;


    using namespace boost::python;

    class_<EquationSystemType, EquationSystemType::Pointer>("EquationSystem")
    //.def(init<std::string const&>())
    .def(init<EquationSystemType::SizeType>())
    .def("Initialize",pointer_to_initalize)
    .def("ApplyDirichletConditions",&EquationSystemType::ApplyDirichletConditions)
    .def("Size",&EquationSystemType::Size)
    .def("DirichletSize",&EquationSystemType::DirichletSize)
    .add_property("SystemMatrix", make_function(pointer_to_get_system_matrix, return_internal_reference<>()),&EquationSystemType::SetSystemMatrix)
    .add_property("Results", make_function(pointer_to_get_results, return_internal_reference<>()),&EquationSystemType::SetResults)
    .add_property("RightHandSide", make_function(pointer_to_get_rhs, return_internal_reference<>()),&EquationSystemType::SetRightHandSide)
    .add_property("DirichletMatrix", make_function(pointer_to_get_dirichlet_matrix, return_internal_reference<>()))
    //.def("",&EquationSystemType::)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

