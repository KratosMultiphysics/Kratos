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
//


// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "python/add_factories_to_python.h"
#include "includes/linear_solver_factory.h"
#include "includes/preconditioner_factory.h"

namespace Kratos
{

namespace Python
{

void  AddFactoriesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
    typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;

    //////////////////////////////////////////////////////////////7
    //HERE THE TOOLS TO REGISTER LINEAR SOLVERS

    typedef LinearSolverFactory< SpaceType, LocalSpaceType > LinearSolverFactoryType;
    typedef LinearSolverFactory< ComplexSpaceType, ComplexLocalSpaceType > ComplexLinearSolverFactoryType;
    typedef PreconditionerFactory< SpaceType, LocalSpaceType > PreconditionerFactoryType;

    py::class_<LinearSolverFactoryType, LinearSolverFactoryType::Pointer>(m, "LinearSolverFactory")
     .def( py::init< >() )
     .def("Create",&LinearSolverFactoryType::Create)
     .def("Has",&LinearSolverFactoryType::Has)
    ;

    py::class_<ComplexLinearSolverFactoryType, ComplexLinearSolverFactoryType::Pointer>(m, "ComplexLinearSolverFactory")
     .def( py::init< >() )
     .def("Create",&ComplexLinearSolverFactoryType::Create)
     .def("Has",&ComplexLinearSolverFactoryType::Has)
    ;

    py::class_<PreconditionerFactoryType, PreconditionerFactoryType::Pointer >(m, "PreconditionerFactory")
     .def( py::init< >() )
     .def("Create",&PreconditionerFactoryType::Create)
    ;

}

}  // namespace Python.

} // Namespace Kratos

