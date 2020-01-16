//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes
// #include <complex>

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include "boost/numeric/ublas/vector.hpp"


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/custom_strategies/linear_mor_matrix_output_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_online_second_order_strategy.hpp"
#include "custom_strategies/custom_strategies/frequency_response_analysis_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_second_order_krylov_strategy.hpp"

// Builders and solvers
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using complex = std::complex<double>;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef TUblasSparseSpace<complex> ComplexSpaceType;
    typedef TUblasDenseSpace<complex> ComplexLocalSpaceType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef LinearSolverType::Pointer LinearSolverPointer;
    typedef LinearSolver<ComplexSpaceType, ComplexLocalSpaceType> ComplexLinearSolverType;
    typedef ComplexLinearSolverType::Pointer ComplexLinearSolverPointer;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef SolvingStrategy< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType > BaseComplexSolvingStrategyType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType > BaseMixedSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

    // Custom strategy types
    typedef LinearMorMatrixOutputStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > LinearMorMatrixOutputStrategyType;
    typedef FrequencyResponseAnalysisStrategy < SparseSpaceType, LocalSpaceType, LinearSolverType > FrequencyResponseAnalysisStrategyType;
    typedef MorOfflineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, SparseSpaceType, LocalSpaceType > MorSecondOrderRealInRealOutOfflineStrategyType;
    typedef MorOfflineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderRealInComplexOutOfflineStrategyType;
    // typedef MorOfflineSecondOrderStrategy< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderComplexInComplexOutOfflineStrategyType;
    // typedef MorOfflineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, ComplexLocalSpaceType > MorOfflineSecondOrderComplexStrategyType;
    typedef MorOnlineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, MorSecondOrderRealInRealOutOfflineStrategyType > MorSecondOrderRealOnlineStrategyType;
    typedef MorOnlineSecondOrderStrategy< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType, MorSecondOrderRealInComplexOutOfflineStrategyType > MorSecondOrderComplexOnlineStrategyType;
    // typedef MorOnlineStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, MorOfflineSecondOrderStrategyType > MorOnlineStrategyType;

    typedef MorSecondOrderKrylovStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, SparseSpaceType, LocalSpaceType > MorSecondOrderKrylovStrategyType;
    typedef MorSecondOrderKrylovStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderComplexKrylovStrategyType;

    // Custom builder and solver types
    typedef SystemMatrixBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > SystemMatrixBuilderAndSolverType;
    // typedef SystemMatrixBuilderAndSolver< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType > ComplexSystemMatrixBuilderAndSolverType;


    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
    py::class_< BaseComplexSolvingStrategyType, typename BaseComplexSolvingStrategyType::Pointer >(m,"BaseComplexSolvingStrategy")
        .def("SetEchoLevel", &BaseComplexSolvingStrategyType::SetEchoLevel)
        ;
    py::class_< BaseMixedSolvingStrategyType, typename BaseMixedSolvingStrategyType::Pointer >(m,"BaseMixedSolvingStrategy")
        .def("Check", &BaseMixedSolvingStrategyType::Check)
        .def("SetEchoLevel", &BaseMixedSolvingStrategyType::SetEchoLevel)
        .def("Solve", &BaseMixedSolvingStrategyType::Solve)
        ;
    py::class_< MorSecondOrderRealInRealOutOfflineStrategyType, typename MorSecondOrderRealInRealOutOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorSecondOrderRealInRealOutOfflineStrategy")
        .def("EchoInfo", &MorSecondOrderRealInRealOutOfflineStrategyType::EchoInfo);
    py::class_< MorSecondOrderRealInComplexOutOfflineStrategyType, typename MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorSecondOrderRealInComplexOutOfflineStrategy")
        .def("EchoInfo", &MorSecondOrderRealInComplexOutOfflineStrategyType::EchoInfo)
        .def("GetK", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetSystemMatrix)
        .def("GetD", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetDampingMatrix)
        .def("GetM", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetMassMatrix)
        .def("GetRHS", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetSystemVector)
        .def("GetOutputVector", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetOutputVector)
        .def("GetBasis", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetBasis)
        .def("GetKr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetKr)
        .def("GetDr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetDr)
        .def("GetMr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetMr)
        .def("GetOutputVectorR", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetOVr)
        ;
    // py::class_< MorSecondOrderComplexInComplexOutOfflineStrategyType, typename MorSecondOrderComplexInComplexOutOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorSecondOrderComplexInComplexOutOfflineStrategy");


    py::class_< LinearMorMatrixOutputStrategyType, typename LinearMorMatrixOutputStrategyType::Pointer, BaseSolvingStrategyType >(m,"LinearMorMatrixOutputStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, bool >())
        ;

    py::class_< MorSecondOrderRealOnlineStrategyType, typename MorSecondOrderRealOnlineStrategyType::Pointer, BaseMixedSolvingStrategyType >(m,"MorOnlineStrategy")
        .def(py::init < ModelPart&, ComplexLinearSolverPointer, MorSecondOrderRealInRealOutOfflineStrategyType::Pointer >())
        .def("Check", &MorSecondOrderRealOnlineStrategyType::Check)
        // .def(py::init < ModelPart&, LinearSolverPointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer >())
        ;

    py::class_< MorSecondOrderComplexOnlineStrategyType, typename MorSecondOrderComplexOnlineStrategyType::Pointer, BaseComplexSolvingStrategyType >(m,"MorComplexOnlineStrategy")
        .def(py::init < ModelPart&, ComplexLinearSolverPointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer > ())
        .def(py::init < ModelPart&, ComplexLinearSolverPointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer, bool > ())
        .def("Check", &MorSecondOrderComplexOnlineStrategyType::Check)
        .def("Solve", &MorSecondOrderComplexOnlineStrategyType::Solve)
        .def("GetScalarResult", &MorSecondOrderComplexOnlineStrategyType::GetScalarResult)
        ;

    // py::class_< MorSecondOrderComplexOnlineStrategyType, typename MorSecondOrderComplexOnlineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorComplexOnlineStrategy");
        // .def(py::init < ModelPart&, ComplexLinearSolverPointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer >())
        // .def(py::init < ModelPart&, LinearSolverPointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer >())
        // ;

    // py::class_< MorOfflineSecondOrderStrategyType, typename MorOfflineSecondOrderStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorOfflineSecondOrderStrategy")
    //     .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, bool >())
    //     .def("SetBuilderAndSolver", &MorOfflineSecondOrderStrategyType::SetBuilderAndSolver)
    //     ;

    // py::class_< MorSecondOrderRealOnlineStrategyType, typename MorSecondOrderRealOnlineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorOnlineSecondOrderStrategy")
    //     .def(py::init < ModelPart&, ComplexLinearSolverPointer, MorOfflineSecondOrderStrategyType::Pointer, bool >())
    //     ;


    //TODO: needs a distinct file
    // py::class_< MorSecondOrderKrylovStrategyType, typename MorSecondOrderKrylovStrategyType::Pointer, MorSecondOrderRealInRealOutOfflineStrategyType >(m,"MorSecondOrderRealKrylovStrategy")
    //     .def(py::init < ModelPart&, BaseSchemeType::Pointer, SystemMatrixBuilderAndSolverType::Pointer, LinearSolverPointer, vector<double>, bool >())
    //     ;

    py::class_< MorSecondOrderComplexKrylovStrategyType, typename MorSecondOrderComplexKrylovStrategyType::Pointer, MorSecondOrderRealInComplexOutOfflineStrategyType >(m,"MorSecondOrderComplexKrylovStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, SystemMatrixBuilderAndSolverType::Pointer, ComplexLinearSolverPointer, vector<double>, bool >())
        ;

    py::class_< FrequencyResponseAnalysisStrategyType, typename FrequencyResponseAnalysisStrategyType::Pointer, BaseSolvingStrategyType >(m,"FrequencyResponseAnalysisStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ComplexLinearSolverPointer, bool, bool >())
        .def("GetBuilderAndSolver", &FrequencyResponseAnalysisStrategyType::GetBuilderAndSolver)
        ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    py::class_< SystemMatrixBuilderAndSolverType, typename SystemMatrixBuilderAndSolverType::Pointer, BuilderAndSolverType >(m, "SystemMatrixBuilderAndSolver")
        .def(py::init < LinearSolverPointer >())
        .def(py::init < LinearSolverPointer, Parameters >())
        .def("BuildEquationIdVector", &SystemMatrixBuilderAndSolverType::BuildEquationIdVector)
        .def("BuildOutputStructure", &SystemMatrixBuilderAndSolverType::BuildOutputStructure)
        ;

    // py::class_< ComplexSystemMatrixBuilderAndSolverType, typename ComplexSystemMatrixBuilderAndSolverType::Pointer, BuilderAndSolverType >(m, "SystemMatrixBuilderAndSolver")
    //     .def(py::init < ComplexLinearSolverPointer >())
    //     .def(py::init < ComplexLinearSolverPointer, Parameters >())
    //     .def("BuildEquationIdVector", &SystemMatrixBuilderAndSolverType::BuildEquationIdVector)
    //     .def("BuildOutputStructure", &SystemMatrixBuilderAndSolverType::BuildOutputStructure)
    //     ;

}

} // namespace Python.
} // Namespace Kratos
