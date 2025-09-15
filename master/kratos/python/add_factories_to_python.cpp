//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_factories_to_python.h"
#include "factories/linear_solver_factory.h"
#include "factories/preconditioner_factory.h"
#include "factories/register_factories.h"

namespace Kratos::Python
{

using SpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SpaceType, LocalSpaceType>;
using ComplexSpaceType = TUblasSparseSpace<std::complex<double>>;
using ComplexLocalSpaceType = TUblasDenseSpace<std::complex<double>>;

using LinearSolverFactoryType = LinearSolverFactory<SpaceType, LocalSpaceType>;
using ComplexLinearSolverFactoryType = LinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType>;
using PreconditionerFactoryType = PreconditionerFactory<SpaceType, LocalSpaceType>;
using ExplicitBuilderType = ExplicitBuilder<SpaceType, LocalSpaceType>;
using ExplicitBuilderFactoryType = Factory<ExplicitBuilderType>;

void  AddFactoriesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER LINEAR SOLVERS
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
     .def("Has",&PreconditionerFactoryType::Has)
    ;

    //////////////////////////////////////////////////////////////
    //HERE WE REGISTER SOME COMMON METHODS
    py::class_<FactoryBase, FactoryBase::Pointer >(m, "FactoryBase")
     .def("Has",&FactoryBase::Has)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER EXPLICIT BUILDER
    py::class_<ExplicitBuilderFactoryType, ExplicitBuilderFactoryType::Pointer, FactoryBase>(m, "ExplicitBuilderFactory")
     .def( py::init< >() )
     .def("Create",[](ExplicitBuilderFactoryType& rExplicitBuilderFactory, Kratos::Parameters Settings) {return rExplicitBuilderFactory.Create(Settings);})
     .def("__str__", PrintObject<ExplicitBuilderFactoryType>)
    ;

}

}  // namespace Kratos::Python.
