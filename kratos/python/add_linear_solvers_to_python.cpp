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
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_linear_solvers_to_python.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "includes/dof.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_complex_interface.h"

#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "linear_solvers/fallback_linear_solver.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"
#include "linear_solvers/power_iteration_highest_eigenvalue_solver.h"
#include "linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h"

namespace Kratos::Python
{
    template <class TDataType>
    using TSpaceType = UblasSpace<TDataType, boost::numeric::ublas::compressed_matrix<TDataType>, boost::numeric::ublas::vector<TDataType>>;
    template <class TDataType>
    using TLocalSpaceType = UblasSpace<TDataType, DenseMatrix<TDataType>, DenseVector<TDataType>>;
    template <class TDataType, class TOtherDataType>
    using TLinearSolverType = LinearSolver<TSpaceType<TDataType>, TLocalSpaceType<TOtherDataType>>;
    template <class TDataType>
    using TDirectSolverType = DirectSolver<TUblasSparseSpace<TDataType>, TUblasDenseSpace<TDataType>>;

void  AddLinearSolversToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using ComplexSpaceType = TUblasSparseSpace<std::complex<double>>;
    using ComplexLocalSpaceType = TUblasDenseSpace<std::complex<double>>;

    using LinearSolverType = LinearSolver<SpaceType, LocalSpaceType>;
    using IterativeSolverType = IterativeSolver<SpaceType, LocalSpaceType>;
    using CGSolverType = CGSolver<SpaceType, LocalSpaceType>;
    using DeflatedCGSolverType = DeflatedCGSolver<SpaceType, LocalSpaceType>;
    using BICGSTABSolverType = BICGSTABSolver<SpaceType, LocalSpaceType>;
    using TFQMRSolverType = TFQMRSolver<SpaceType, LocalSpaceType>;
    using ScalingSolverType = ScalingSolver<SpaceType, LocalSpaceType>;
    using PowerIterationEigenvalueSolverType = PowerIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType>;
    using PowerIterationHighestEigenvalueSolverType = PowerIterationHighestEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType>;
    using RayleighQuotientIterationEigenvalueSolverType = RayleighQuotientIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType>;

    using ComplexLinearSolverType = TLinearSolverType<std::complex<double>, std::complex<double>>;
    using MixedLinearSolverType = TLinearSolverType<double, std::complex<double>>;
    using ComplexDirectSolverType = TDirectSolverType<std::complex<double>>;
    using ComplexSkylineLUSolverType = SkylineLUCustomScalarSolver<ComplexSpaceType, ComplexLocalSpaceType>;

    bool (LinearSolverType::*pointer_to_solve)(LinearSolverType::SparseMatrixType& rA, LinearSolverType::VectorType& rX, LinearSolverType::VectorType& rB) = &LinearSolverType::Solve;
    void (LinearSolverType::*pointer_to_solve_eigen)(LinearSolverType::SparseMatrixType& rK, LinearSolverType::SparseMatrixType& rM,LinearSolverType::DenseVectorType& Eigenvalues, LinearSolverType::DenseMatrixType& Eigenvectors) = &LinearSolverType::Solve;
    bool (ComplexLinearSolverType::*pointer_to_complex_solve)(ComplexLinearSolverType::SparseMatrixType& rA, ComplexLinearSolverType::VectorType& rX, ComplexLinearSolverType::VectorType& rB) = &ComplexLinearSolverType::Solve;
    void (ComplexLinearSolverType::*pointer_to_complex_solve_eigen)(ComplexLinearSolverType::SparseMatrixType& rK, ComplexLinearSolverType::SparseMatrixType& rM, ComplexLinearSolverType::DenseVectorType& Eigenvalues, ComplexLinearSolverType::DenseMatrixType& Eigenvectors) = &ComplexLinearSolverType::Solve;
    void (MixedLinearSolverType::*pointer_to_mixed_solve_eigen)(MixedLinearSolverType::SparseMatrixType& rK, MixedLinearSolverType::SparseMatrixType& rM,MixedLinearSolverType::DenseVectorType& Eigenvalues, MixedLinearSolverType::DenseMatrixType& Eigenvectors) = &MixedLinearSolverType::Solve;

    //****************************************************************************************************
    //preconditioners
    //****************************************************************************************************
    using PreconditionerType = Preconditioner<SpaceType,  LocalSpaceType>;

    py::class_<PreconditionerType, PreconditionerType::Pointer>(m,"Preconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<PreconditionerType>)
    ;

    using DiagonalPreconditionerType = DiagonalPreconditioner<SpaceType,  LocalSpaceType>;
    py::class_<DiagonalPreconditionerType, DiagonalPreconditionerType::Pointer, PreconditionerType>(m,"DiagonalPreconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<DiagonalPreconditionerType>)
    ;

    using ILUPreconditionerType = ILUPreconditioner<SpaceType,  LocalSpaceType>;
    py::class_<ILUPreconditionerType, ILUPreconditionerType::Pointer, PreconditionerType>(m,"ILUPreconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<ILUPreconditionerType>)
    ;

    using ILU0PreconditionerType = ILU0Preconditioner<SpaceType,  LocalSpaceType>;
    py::class_<ILU0PreconditionerType, ILU0PreconditionerType::Pointer, PreconditionerType>(m,"ILU0Preconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<ILU0PreconditionerType>)
    ;

    //****************************************************************************************************
    //linear solvers
    //****************************************************************************************************
    py::class_<LinearSolverType, LinearSolverType::Pointer>(m,"LinearSolver")
    .def(py::init< >() )
    .def("Initialize",&LinearSolverType::Initialize)
    .def("Solve",pointer_to_solve)
    .def("Solve",pointer_to_solve_eigen)
    .def("Clear",&LinearSolverType::Clear)
    .def("__str__", PrintObject<LinearSolverType>)
    .def( "GetIterationsNumber",&LinearSolverType::GetIterationsNumber)
    ;

    py::class_<ComplexLinearSolverType, ComplexLinearSolverType::Pointer>(m,"ComplexLinearSolver")
    .def(py::init< >() )
    .def("Initialize",&ComplexLinearSolverType::Initialize)
    .def("Solve",pointer_to_complex_solve)
    .def("Solve",pointer_to_complex_solve_eigen)
    .def("Clear",&ComplexLinearSolverType::Clear)
    .def("__str__", PrintObject<ComplexLinearSolverType>)
    ;

    py::class_<MixedLinearSolverType, MixedLinearSolverType::Pointer>(m,"MixedLinearSolver")
    .def(py::init< >() )
    .def("Initialize",&MixedLinearSolverType::Initialize)
    .def("Solve",pointer_to_mixed_solve_eigen)
    .def("Clear",&MixedLinearSolverType::Clear)
    .def("__str__", PrintObject<MixedLinearSolverType>)
    ;

    py::class_<IterativeSolverType, IterativeSolverType::Pointer, LinearSolverType>(m,"IterativeSolver")
    .def(py::init< >() )
    .def("__str__", PrintObject<IterativeSolverType>)
    ;

    py::class_<CGSolverType, CGSolverType::Pointer,IterativeSolverType>(m,"CGSolver")
    .def(py::init<double>())
    .def(py::init<double, unsigned int>())
    .def(py::init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(py::init<Parameters,  PreconditionerType::Pointer>())
    .def("__str__", PrintObject<CGSolverType>)
    ;

    py::class_<BICGSTABSolverType, BICGSTABSolverType::Pointer,IterativeSolverType>(m,"BICGSTABSolver")
    .def(py::init<double>())
    .def(py::init<double, unsigned int>())
    .def("__str__", PrintObject<BICGSTABSolverType>)
    .def(py::init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(py::init<Parameters,  PreconditionerType::Pointer>())
    .def("SetTolerance",&BICGSTABSolverType::SetTolerance)
    ;

    py::class_<TFQMRSolverType, TFQMRSolverType::Pointer,IterativeSolverType>(m,"TFQMRSolver")
    .def(py::init<double>())
    .def(py::init<double, unsigned int>())
    .def("__str__", PrintObject<TFQMRSolverType>)
    .def(py::init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(py::init<Parameters,  PreconditionerType::Pointer>())
    ;

    py::class_<ScalingSolverType, ScalingSolverType::Pointer, LinearSolverType>(m,"ScalingSolver")
    .def(py::init<LinearSolverType::Pointer>())
    .def(py::init<LinearSolverType::Pointer, bool >())
    .def(py::init<Parameters >())
    ;

    py::class_<PowerIterationEigenvalueSolverType, PowerIterationEigenvalueSolverType::Pointer, LinearSolverType>(m,"PowerIterationEigenvalueSolver")
    .def(py::init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
    .def(py::init<Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<PowerIterationHighestEigenvalueSolverType, PowerIterationHighestEigenvalueSolverType::Pointer, LinearSolverType>(m,"PowerIterationHighestEigenvalueSolver")
    .def(py::init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
    .def(py::init<Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<RayleighQuotientIterationEigenvalueSolverType, RayleighQuotientIterationEigenvalueSolverType::Pointer, LinearSolverType>(m,"RayleighQuotientIterationEigenvalueSolver")
    .def(py::init<double, unsigned int, unsigned int, LinearSolverType::Pointer, double>())
    .def(py::init<Parameters, LinearSolverType::Pointer>())
    ;

    using ReordererType = Reorderer<SpaceType, LocalSpaceType>;
    using DirectSolverType = DirectSolver<SpaceType, LocalSpaceType, ReordererType>;
    using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SpaceType, LocalSpaceType, ReordererType>;

    py::class_<ReordererType, ReordererType::Pointer>(m,"Reorderer")
    .def(py::init< >() )
    .def("__str__", PrintObject<ReordererType>)
    .def( "Initialize",&ReordererType::Initialize)
    .def( "Reorder",&ReordererType::Reorder)
    .def( "InverseReorder",&ReordererType::InverseReorder)
    ;

    py::class_<DirectSolverType, DirectSolverType::Pointer, LinearSolverType>(m,"DirectSolver")
    .def(py::init< >() )
    .def(py::init<Parameters>())
    .def("__str__", PrintObject<DirectSolverType>)
    ;

    py::class_<ComplexDirectSolverType, ComplexDirectSolverType::Pointer, ComplexLinearSolverType>(m,"ComplexDirectSolver")
    .def(py::init< >() )
    .def(py::init<Parameters>())
    .def("__str__", PrintObject<ComplexDirectSolverType>)
    ;

    py::class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, DirectSolverType>(m,"SkylineLUFactorizationSolver")
    .def(py::init< >())
    .def(py::init<Parameters>())
    .def("__str__", PrintObject<SkylineLUFactorizationSolverType>)
    ;

    py::class_<ComplexSkylineLUSolverType, typename ComplexSkylineLUSolverType::Pointer, ComplexDirectSolverType>(m,"ComplexSkylineLUSolver")
    .def(py::init< >())
    .def(py::init<Parameters&>())
    .def("__str__", PrintObject<ComplexSkylineLUSolverType>)
    ;

    py::class_<DeflatedCGSolverType, DeflatedCGSolverType::Pointer,IterativeSolverType>(m,"DeflatedCGSolver")
    .def(py::init<double,bool,int>())
    .def(py::init<double, unsigned int,bool,int>())
    .def(py::init<double, unsigned int,  PreconditionerType::Pointer,bool,int>())
    .def(py::init<Parameters>())
// 		  .def(py::init<double, unsigned int,  PreconditionerType::Pointer, ModelPart*>())
    //.def("",&LinearSolverType::)
    .def("__str__", PrintObject<DeflatedCGSolverType>)
    ;

    using FallbackLinearSolverType = FallbackLinearSolver<SpaceType, LocalSpaceType>;
    py::class_<FallbackLinearSolverType, FallbackLinearSolverType::Pointer, LinearSolverType>(m, "FallbackLinearSolver")
    .def(py::init<Parameters>())
    .def(py::init<const std::vector<LinearSolverType::Pointer>&, Parameters>())
    // .def("AddSolver", [](FallbackLinearSolverType& rSelf, LinearSolverType::Pointer pSolver) {
    //     rSelf.AddSolver(pSolver);
    // })
    // .def("AddSolver", [](FallbackLinearSolverType& rSelf, const Parameters ThisParameters) {
    //     rSelf.AddSolver(ThisParameters);
    // })
    .def("GetSolvers", &FallbackLinearSolverType::GetSolvers)
    // .def("SetSolvers", &FallbackLinearSolverType::SetSolvers)
    .def("GetResetSolverEachTry", &FallbackLinearSolverType::GetResetSolverEachTry)
    // .def("SetResetSolverIndexEachTry", &FallbackLinearSolverType::SetResetSolverIndexEachTry)
    .def("GetParameters", &FallbackLinearSolverType::GetParameters)
    .def("GetCurrentSolverIndex", &FallbackLinearSolverType::GetCurrentSolverIndex)
    .def("ClearCurrentSolverIndex", &FallbackLinearSolverType::ClearCurrentSolverIndex)
    ;

}

}  // namespace Kratos::Python.
