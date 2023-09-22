//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//
//


// System includes


// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"


#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_strategies/rom_builder_and_solver.h"
#include "custom_strategies/lspg_rom_builder_and_solver.h"
#include "custom_strategies/petrov_galerkin_rom_builder_and_solver.h"
#include "custom_strategies/global_rom_builder_and_solver.h"

//linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

    //********************************************************************
    //********************************************************************
    typedef ROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ROMBuilderAndSolverType;

     py::class_<ROMBuilderAndSolverType, typename ROMBuilderAndSolverType::Pointer, BuilderAndSolverType>(m, "ROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        .def("getqn", &ROMBuilderAndSolverType::getqn)    // for testing purposes, remove before release
        .def("setqn", &ROMBuilderAndSolverType::setqn)    // for testing purposes, remove before release
        .def("GetQs", &ROMBuilderAndSolverType::GetQs)    // for testing purposes, remove before release
        .def("GetUs", &ROMBuilderAndSolverType::GetUs)    // for testing purposes, remove before release
        .def("TestProjectToFineBasis", &ROMBuilderAndSolverType::TestProjectToFineBasis)   //// for testing purposes, remove before release
        .def("setQuadratic", &ROMBuilderAndSolverType::setQuadratic)   //// for testing purposes, remove before release
        .def("getQuadratic", &ROMBuilderAndSolverType::getQuadratic)   //// for testing purposes, remove before release
        ;

    typedef LeastSquaresPetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> LeastSquaresPetrovGalerkinROMBuilderAndSolverType;

     py::class_<LeastSquaresPetrovGalerkinROMBuilderAndSolverType, typename LeastSquaresPetrovGalerkinROMBuilderAndSolverType::Pointer, ROMBuilderAndSolverType, BuilderAndSolverType>(m, "LeastSquaresPetrovGalerkinROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;
    
    typedef PetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> PetrovGalerkinROMBuilderAndSolverType;

     py::class_<PetrovGalerkinROMBuilderAndSolverType, typename PetrovGalerkinROMBuilderAndSolverType::Pointer, ROMBuilderAndSolverType, BuilderAndSolverType>(m, "PetrovGalerkinROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;

    typedef GlobalROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> GlobalROMBuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;
    
    py::class_<GlobalROMBuilderAndSolverType, typename GlobalROMBuilderAndSolverType::Pointer, ResidualBasedBlockBuilderAndSolverType>(m, "GlobalROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        .def("getqn", &GlobalROMBuilderAndSolverType::getqn)    // for testing purposes, remove before release
        .def("setqn", &GlobalROMBuilderAndSolverType::setqn)    // for testing purposes, remove before release
        .def("setQuadratic", &GlobalROMBuilderAndSolverType::setQuadratic)   //// for testing purposes, remove before release
        .def("getQuadratic", &GlobalROMBuilderAndSolverType::getQuadratic)   //// for testing purposes, remove before release
        .def("GetQs", &GlobalROMBuilderAndSolverType::GetQs)    // for testing purposes
        .def("GetUs", &GlobalROMBuilderAndSolverType::GetUs)    // for testing purposes
        ;


}

} // namespace Python.
} // Namespace Kratos

