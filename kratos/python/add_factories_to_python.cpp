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

namespace Kratos
{

namespace Python
{

typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;

typedef LinearSolverFactory< SpaceType, LocalSpaceType > LinearSolverFactoryType;
typedef LinearSolverFactory< ComplexSpaceType, ComplexLocalSpaceType > ComplexLinearSolverFactoryType;
typedef PreconditionerFactory< SpaceType, LocalSpaceType > PreconditionerFactoryType;
typedef ExplicitBuilder< SpaceType, LocalSpaceType > ExplicitBuilderType;
typedef Factory< ExplicitBuilderType > ExplicitBuilderFactoryType;
typedef ExplicitSolvingStrategy< SpaceType, LocalSpaceType > ExplicitSolvingStrategyType;
typedef Factory< ExplicitSolvingStrategyType > ExplicitStrategyFactoryType;

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

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER EXPLICIT STRATEGIES
    py::class_<ExplicitStrategyFactoryType, ExplicitStrategyFactoryType::Pointer, FactoryBase>(m, "ExplicitStrategyFactory")
     .def( py::init< >() )
     .def("Create",[](ExplicitStrategyFactoryType& rExplicitStrategyFactory, ModelPart& rModelPart, Kratos::Parameters Settings) {return rExplicitStrategyFactory.Create(rModelPart, Settings);})
     .def("__str__", PrintObject<ExplicitStrategyFactoryType>)
    ;

}

}  // namespace Python.

} // Namespace Kratos
