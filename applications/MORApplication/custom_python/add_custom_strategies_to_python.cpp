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


// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/custom_strategies/linear_mor_matrix_output_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_online_second_order_strategy.hpp"
#include "custom_strategies/custom_strategies/frequency_response_analysis_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_second_order_krylov_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_second_order_irka_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_second_order_toar_strategy.hpp"
#include "custom_strategies/custom_strategies/modal_analysis_strategy.hpp"

// Builders and solvers
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// Schemes
#include "custom_strategies/custom_schemes/matrix_builder_scheme.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using complex = std::complex<double>;

    // typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    // typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef TUblasSparseSpace<double> SparseSpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;
    typedef TUblasSparseSpace<complex> ComplexSpaceType;
    typedef TUblasDenseSpace<complex> ComplexLocalSpaceType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef LinearSolverType::Pointer LinearSolverPointer;
    typedef LinearSolver<LocalSpaceType, LocalSpaceType> DenseLinearSolverType;
    typedef LinearSolver<ComplexSpaceType, ComplexLocalSpaceType> ComplexLinearSolverType;
    typedef LinearSolver<ComplexLocalSpaceType, ComplexLocalSpaceType> ComplexDenseLinearSolverType;
    typedef LinearSolver<SparseSpaceType, ComplexLocalSpaceType> MixedLinearSolverType;
    typedef ComplexLinearSolverType::Pointer ComplexLinearSolverPointer;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef SolvingStrategy< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType > BaseComplexSolvingStrategyType;
    typedef SolvingStrategy< ComplexLocalSpaceType, ComplexLocalSpaceType, ComplexDenseLinearSolverType > BaseComplexDenseSolvingStrategyType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType > BaseMixedSolvingStrategyType;
    typedef SolvingStrategy< LocalSpaceType, LocalSpaceType, ComplexLinearSolverType > BaseMixedDenseSolvingStrategyType;
    typedef SolvingStrategy< LocalSpaceType, LocalSpaceType, DenseLinearSolverType > BaseDenseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

    // Custom strategy types
    typedef LinearMorMatrixOutputStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > LinearMorMatrixOutputStrategyType;
    typedef FrequencyResponseAnalysisStrategy < SparseSpaceType, LocalSpaceType, LinearSolverType > FrequencyResponseAnalysisStrategyType;
                                                                            //this looks weird
    typedef MorOfflineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, SparseSpaceType, LocalSpaceType > MorSecondOrderRealInRealOutOfflineStrategyType;
    typedef MorOfflineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderRealInComplexOutOfflineStrategyType;
    // typedef MorOfflineSecondOrderStrategy< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderComplexInComplexOutOfflineStrategyType;
    // typedef MorOfflineSecondOrderStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, ComplexLocalSpaceType > MorOfflineSecondOrderComplexStrategyType;
    typedef MorOnlineSecondOrderStrategy< LocalSpaceType, LocalSpaceType, ComplexDenseLinearSolverType, MorSecondOrderRealInRealOutOfflineStrategyType > MorSecondOrderRealOnlineStrategyType;
    typedef MorOnlineSecondOrderStrategy< ComplexLocalSpaceType, ComplexLocalSpaceType, ComplexDenseLinearSolverType, MorSecondOrderRealInComplexOutOfflineStrategyType > MorSecondOrderComplexOnlineStrategyType;
    // typedef MorOnlineStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, MorOfflineSecondOrderStrategyType > MorOnlineStrategyType;

    // typedef MorSecondOrderKrylovStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType, SparseSpaceType, LocalSpaceType > MorSecondOrderKrylovStrategyType;
    typedef MorSecondOrderKrylovStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderComplexKrylovStrategyType;

    typedef MorSecondOrderIRKAStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, SparseSpaceType, LocalSpaceType, false > MorSecondOrderIrkaRealStrategyType;
    typedef MorSecondOrderIRKAStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType, true > MorSecondOrderIrkaComplexStrategyType;

    typedef MorSecondOrderTOARStrategy< SparseSpaceType, LocalSpaceType, ComplexLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType > MorSecondOrderComplexTOARStrategyType;

    typedef ModalAnalysisStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType, SparseSpaceType, LocalSpaceType> ModalAnalysisRealInRealOutStrategyType;
    typedef ModalAnalysisStrategy<SparseSpaceType, LocalSpaceType, MixedLinearSolverType, ComplexSpaceType, ComplexLocalSpaceType> ModalAnalysisRealInComplexOutStrategyType;

    // Custom builder and solver types
    typedef SystemMatrixBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > SystemMatrixBuilderAndSolverType;
    // typedef SystemMatrixBuilderAndSolver< ComplexSpaceType, ComplexLocalSpaceType, ComplexLinearSolverType > ComplexSystemMatrixBuilderAndSolverType;

    // Custom scheme types
    typedef MatrixBuilderScheme< SparseSpaceType, LocalSpaceType > MatrixBuilderSchemeType;

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
    py::class_< BaseComplexDenseSolvingStrategyType, typename BaseComplexDenseSolvingStrategyType::Pointer >(m,"BaseComplexDenseSolvingStrategy")
        .def("SetEchoLevel", &BaseComplexDenseSolvingStrategyType::SetEchoLevel)
        ;
    py::class_< BaseComplexSolvingStrategyType, typename BaseComplexSolvingStrategyType::Pointer >(m,"BaseComplexSolvingStrategy")
        .def("SetEchoLevel", &BaseComplexSolvingStrategyType::SetEchoLevel)
        ;
    py::class_< BaseMixedSolvingStrategyType, typename BaseMixedSolvingStrategyType::Pointer >(m,"BaseMixedSolvingStrategy")
        .def("Check", &BaseMixedSolvingStrategyType::Check)
        .def("SetEchoLevel", &BaseMixedSolvingStrategyType::SetEchoLevel)
        .def("Solve", &BaseMixedSolvingStrategyType::Solve)
        ;
    py::class_< BaseMixedDenseSolvingStrategyType, typename BaseMixedDenseSolvingStrategyType::Pointer >(m,"BaseMixedDenseSolvingStrategy")
        .def("Check", &BaseMixedDenseSolvingStrategyType::Check)
        .def("SetEchoLevel", &BaseMixedDenseSolvingStrategyType::SetEchoLevel)
        .def("Solve", &BaseMixedDenseSolvingStrategyType::Solve)
        ;
    py::class_< BaseDenseSolvingStrategyType, typename BaseDenseSolvingStrategyType::Pointer >(m,"BaseDenseSolvingStrategy")
        .def("Initialize", &BaseDenseSolvingStrategyType::Initialize)
        .def("Check", &BaseDenseSolvingStrategyType::Check)
        .def("SetEchoLevel", &BaseDenseSolvingStrategyType::SetEchoLevel)
        .def("Solve", &BaseDenseSolvingStrategyType::Solve)
        ;

    py::class_< MorSecondOrderRealInRealOutOfflineStrategyType, typename MorSecondOrderRealInRealOutOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorSecondOrderRealInRealOutOfflineStrategy")
        .def("EchoInfo", &MorSecondOrderRealInRealOutOfflineStrategyType::EchoInfo)
        .def("GetK", &MorSecondOrderRealInRealOutOfflineStrategyType::GetK)
        .def("GetD", &MorSecondOrderRealInRealOutOfflineStrategyType::GetC)
        .def("GetM", &MorSecondOrderRealInRealOutOfflineStrategyType::GetM)
        .def("GetBasis", &MorSecondOrderRealInRealOutOfflineStrategyType::GetBasis)
        .def("GetKr", &MorSecondOrderRealInRealOutOfflineStrategyType::GetKr)
        .def("GetDr", &MorSecondOrderRealInRealOutOfflineStrategyType::GetDr)
        .def("GetMr", &MorSecondOrderRealInRealOutOfflineStrategyType::GetMr)
        // .def("ImportSystem", &MorSecondOrderRealInRealOutOfflineStrategyType::ImportSystem)
        ;
    py::class_< MorSecondOrderRealInComplexOutOfflineStrategyType, typename MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorSecondOrderRealInComplexOutOfflineStrategy")
        .def("EchoInfo", &MorSecondOrderRealInComplexOutOfflineStrategyType::EchoInfo)
        .def("GetK", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetK)
        .def("GetKi", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetKi)
        .def("GetD", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetC)
        .def("GetM", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetM)
        .def("GetRHS", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetSystemVector)
        .def("GetOutputVector", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetOutputVector)
        .def("GetBasis", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetBasis)
        .def("GetKr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetKr)
        .def("GetDr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetDr)
        .def("GetMr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetMr)
        .def("GetRHSr", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetRHSr)
        .def("GetOutputVectorR", &MorSecondOrderRealInComplexOutOfflineStrategyType::GetOVr)
        .def("ImportSystem", (void (MorSecondOrderRealInComplexOutOfflineStrategyType::*)
            (typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemVectorType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemVectorType&))
            &MorSecondOrderRealInComplexOutOfflineStrategyType::ImportSystem)
        .def("ImportSystem", (void (MorSecondOrderRealInComplexOutOfflineStrategyType::*)
            (typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemMatrixType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemVectorType&,
             typename MorSecondOrderRealInComplexOutOfflineStrategyType::TSystemVectorType&))
            &MorSecondOrderRealInComplexOutOfflineStrategyType::ImportSystem)
        ;
    // py::class_< MorSecondOrderComplexInComplexOutOfflineStrategyType, typename MorSecondOrderComplexInComplexOutOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorSecondOrderComplexInComplexOutOfflineStrategy");


    py::class_< LinearMorMatrixOutputStrategyType, typename LinearMorMatrixOutputStrategyType::Pointer, BaseSolvingStrategyType >(m,"LinearMorMatrixOutputStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, bool >())
        ;

    py::class_< MorSecondOrderRealOnlineStrategyType, typename MorSecondOrderRealOnlineStrategyType::Pointer, BaseDenseSolvingStrategyType >(m,"MorRealOnlineStrategy")
        .def(py::init < ModelPart&, ComplexDenseLinearSolverType::Pointer, MorSecondOrderRealInRealOutOfflineStrategyType::Pointer >())
        .def(py::init < ModelPart&, ComplexDenseLinearSolverType::Pointer, MorSecondOrderRealInRealOutOfflineStrategyType::Pointer, bool >())
        .def("Check", &MorSecondOrderRealOnlineStrategyType::Check)
        .def("Solve", &MorSecondOrderRealOnlineStrategyType::Solve)
        .def("GetScalarResult", &MorSecondOrderRealOnlineStrategyType::GetScalarResult)
        // .def(py::init < ModelPart&, LinearSolverPointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer >())
        ;

    //THIS HAS TO WORK!!!
    py::class_< MorSecondOrderComplexOnlineStrategyType, typename MorSecondOrderComplexOnlineStrategyType::Pointer, BaseComplexDenseSolvingStrategyType >(m,"MorComplexOnlineStrategy")
        .def(py::init < ModelPart&, ComplexDenseLinearSolverType::Pointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer > ())
        .def(py::init < ModelPart&, ComplexDenseLinearSolverType::Pointer, MorSecondOrderRealInComplexOutOfflineStrategyType::Pointer, bool > ())
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
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, vector<double> >())
        ;

    py::class_< MorSecondOrderIrkaRealStrategyType, typename MorSecondOrderIrkaRealStrategyType::Pointer, MorSecondOrderRealInRealOutOfflineStrategyType >(m,"MorSecondOrderRealIrkaStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexVector, size_t, double >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexLinearSolverPointer, ComplexVector, size_t, double >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexVector, complex, complex, size_t, double >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexLinearSolverPointer, ComplexVector, complex, complex, size_t, double >())
        .def("GetSamplingPoints", &MorSecondOrderIrkaRealStrategyType::GetSamplingPoints)
        ;

    py::class_< MorSecondOrderIrkaComplexStrategyType, typename MorSecondOrderIrkaComplexStrategyType::Pointer, MorSecondOrderRealInComplexOutOfflineStrategyType >(m,"MorSecondOrderComplexIrkaStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexVector, size_t, double >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexLinearSolverPointer, ComplexVector, size_t, double >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexVector, complex, complex, size_t, double >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, ComplexLinearSolverPointer, ComplexVector, complex, complex, size_t, double >())
        .def("GetSamplingPoints", &MorSecondOrderIrkaComplexStrategyType::GetSamplingPoints)
        ;

    py::class_< MorSecondOrderComplexTOARStrategyType, typename MorSecondOrderComplexTOARStrategyType::Pointer, MorSecondOrderRealInComplexOutOfflineStrategyType >(m,"MorSecondOrderComplexTOARStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, std::vector<complex>, std::vector<int> >())
        ;

    //********************************************************************
    //*************************FREQUENCY RESPONSE*************************
    //********************************************************************

    py::class_< FrequencyResponseAnalysisStrategyType, typename FrequencyResponseAnalysisStrategyType::Pointer, BaseSolvingStrategyType >(m,"FrequencyResponseAnalysisStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, ComplexLinearSolverPointer, bool, bool >())
        .def("GetBuilderAndSolver", &FrequencyResponseAnalysisStrategyType::GetBuilderAndSolver)
        .def("EchoInfo", &FrequencyResponseAnalysisStrategyType::EchoInfo)
        ;

    //********************************************************************
    //***************************MODAL ANALYSIS***************************
    //********************************************************************

    py::class_< ModalAnalysisRealInRealOutStrategyType, typename ModalAnalysisRealInRealOutStrategyType::Pointer, BaseSolvingStrategyType >(m,"ModalAnalysisRealInRealOutStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, LinearSolverPointer>())
        .def("GetEigenvalueVector", &ModalAnalysisRealInRealOutStrategyType::GetEigenvalueVector)
        .def("GetEigenvectorMatrix", &ModalAnalysisRealInRealOutStrategyType::GetEigenvectorMatrix)
        ;

    py::class_< ModalAnalysisRealInComplexOutStrategyType, typename ModalAnalysisRealInComplexOutStrategyType::Pointer, BaseSolvingStrategyType >(m,"ModalAnalysisRealInComplexOutStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, MixedLinearSolverType::Pointer>())
        .def("GetEigenvalueVector", &ModalAnalysisRealInComplexOutStrategyType::GetEigenvalueVector)
        .def("GetEigenvectorMatrix", &ModalAnalysisRealInComplexOutStrategyType::GetEigenvectorMatrix)
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

    //********************************************************************
    //******************************SCHEMES*******************************
    //********************************************************************

    py::class_< MatrixBuilderSchemeType, typename MatrixBuilderSchemeType::Pointer, BaseSchemeType >(m, "MatrixBuilderScheme")
        .def(py::init<>())
        ;
}

} // namespace Python.
} // Namespace Kratos
