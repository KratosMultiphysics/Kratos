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
#include "factories/base_factory.h"

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
typedef ConvergenceCriteria< SpaceType, LocalSpaceType > ConvergenceCriteriaType;
typedef BaseFactory< ConvergenceCriteriaType > ConvergenceCriteriaFactoryType;
typedef Scheme< SpaceType, LocalSpaceType > SchemeType;
typedef BaseFactory< SchemeType > SchemeFactoryType;
typedef BuilderAndSolver< SpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
typedef BaseFactory< BuilderAndSolverType, LinearSolverType > BuilderAndSolverFactoryType;
typedef ExplicitBuilder< SpaceType, LocalSpaceType > ExplicitBuilderType;
typedef BaseFactory< ExplicitBuilderType > ExplicitBuilderFactoryType;
typedef SolvingStrategy< SpaceType, LocalSpaceType, LinearSolverType > SolvingStrategyType;
typedef BaseFactory< SolvingStrategyType, ModelPart > StrategyFactoryType;
typedef ExplicitSolvingStrategy< SpaceType, LocalSpaceType > ExplicitSolvingStrategyType;
typedef BaseFactory< ExplicitSolvingStrategyType, ModelPart > ExplicitStrategyFactoryType;

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
    py::class_<BaseFactoryMethods, BaseFactoryMethods::Pointer >(m, "BaseFactoryMethods")
     .def( py::init< >() )
     .def("Has",&BaseFactoryMethods::Has)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER CONVERGENCE CRITERIA
    py::class_<ConvergenceCriteriaFactoryType, ConvergenceCriteriaFactoryType::Pointer, BaseFactoryMethods>(m, "ConvergenceCriteriaFactory")
     .def( py::init< >() )
     .def("Create",[](ConvergenceCriteriaFactoryType& rConvergenceCriteriaFactory, Kratos::Parameters Settings) {return rConvergenceCriteriaFactory.Create(Settings);})
     .def("__str__", PrintObject<ConvergenceCriteriaFactoryType>)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER SCHEMES
    py::class_<SchemeFactoryType, SchemeFactoryType::Pointer, BaseFactoryMethods>(m, "SchemeFactory")
     .def( py::init< >() )
     .def("Create",[](SchemeFactoryType& rSchemeFactory, Kratos::Parameters Settings) {return rSchemeFactory.Create(Settings);})
     .def("__str__", PrintObject<SchemeFactoryType>)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER BUILDER AND SOLVERS
    py::class_<BuilderAndSolverFactoryType, BuilderAndSolverFactoryType::Pointer, BaseFactoryMethods>(m, "BuilderAndSolverFactory")
     .def( py::init< >() )
     .def("Create",[](BuilderAndSolverFactoryType& rBuilderAndSolverFactory, typename LinearSolverType::Pointer pLinearSolver, Kratos::Parameters Settings) {return rBuilderAndSolverFactory.Create(pLinearSolver, Settings);})
     .def("__str__", PrintObject<BuilderAndSolverFactoryType>)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER EXPLICIT BUILDER
    py::class_<ExplicitBuilderFactoryType, ExplicitBuilderFactoryType::Pointer, BaseFactoryMethods>(m, "ExplicitBuilderFactory")
     .def( py::init< >() )
     .def("Create",[](ExplicitBuilderFactoryType& rBuilderAndSolverFactory, Kratos::Parameters Settings) {return rBuilderAndSolverFactory.Create(Settings);})
     .def("__str__", PrintObject<ExplicitBuilderFactoryType>)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER STRATEGIES
    py::class_<StrategyFactoryType, StrategyFactoryType::Pointer, BaseFactoryMethods>(m, "StrategyFactory")
     .def( py::init< >() )
     .def("Create",[](StrategyFactoryType& rStrategyFactory, ModelPart& rModelPart, Kratos::Parameters Settings) {return rStrategyFactory.Create(rModelPart, Settings);})
     .def("__str__", PrintObject<StrategyFactoryType>)
    ;

    //////////////////////////////////////////////////////////////
    //HERE THE TOOLS TO REGISTER EXPLICIT STRATEGIES
    py::class_<ExplicitStrategyFactoryType, ExplicitStrategyFactoryType::Pointer, BaseFactoryMethods>(m, "ExplicitStrategyFactory")
     .def( py::init< >() )
     .def("Create",[](ExplicitStrategyFactoryType& rStrategyFactory, ModelPart& rModelPart, Kratos::Parameters Settings) {return rStrategyFactory.Create(rModelPart, Settings);})
     .def("__str__", PrintObject<ExplicitStrategyFactoryType>)
    ;

}

}  // namespace Python.

} // Namespace Kratos
