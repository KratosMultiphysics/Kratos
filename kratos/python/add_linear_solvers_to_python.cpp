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

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"
#include "linear_solvers/power_iteration_highest_eigenvalue_solver.h"
#include "linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h"

namespace Kratos
{

namespace Python
{
    template <class TDataType>
    using TSpaceType = UblasSpace<TDataType, boost::numeric::ublas::compressed_matrix<TDataType>, boost::numeric::ublas::vector<TDataType>>;
    template <class TDataType>
    using TLocalSpaceType = UblasSpace<TDataType, DenseMatrix<TDataType>, DenseVector<TDataType>>;
    template <class TDataType>
    using TLinearSolverType = LinearSolver<TSpaceType<TDataType>, TLocalSpaceType<TDataType>>;
    template <class TDataType>
    using TDirectSolverType = DirectSolver<TUblasSparseSpace<TDataType>, TUblasDenseSpace<TDataType>>;

void  AddLinearSolversToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
    typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;

    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;
    typedef CGSolver<SpaceType,  LocalSpaceType> CGSolverType;
    typedef DeflatedCGSolver<SpaceType,  LocalSpaceType> DeflatedCGSolverType;
    typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
    typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;
    typedef ScalingSolver<SpaceType,  LocalSpaceType> ScalingSolverType;
    typedef PowerIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> PowerIterationEigenvalueSolverType;
    typedef PowerIterationHighestEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> PowerIterationHighestEigenvalueSolverType;
    typedef RayleighQuotientIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> RayleighQuotientIterationEigenvalueSolverType;

    typedef TLinearSolverType<std::complex<double>> ComplexLinearSolverType;
    typedef TDirectSolverType<std::complex<double>> ComplexDirectSolverType;
    typedef SkylineLUCustomScalarSolver<ComplexSpaceType, ComplexLocalSpaceType> ComplexSkylineLUSolverType;

    bool (LinearSolverType::*pointer_to_solve)(LinearSolverType::SparseMatrixType& rA, LinearSolverType::VectorType& rX, LinearSolverType::VectorType& rB) = &LinearSolverType::Solve;
    void (LinearSolverType::*pointer_to_solve_eigen)(LinearSolverType::SparseMatrixType& rK, LinearSolverType::SparseMatrixType& rM,LinearSolverType::DenseVectorType& Eigenvalues, LinearSolverType::DenseMatrixType& Eigenvectors) = &LinearSolverType::Solve;
    bool (ComplexLinearSolverType::*pointer_to_complex_solve)(ComplexLinearSolverType::SparseMatrixType& rA, ComplexLinearSolverType::VectorType& rX, ComplexLinearSolverType::VectorType& rB) = &ComplexLinearSolverType::Solve;

    //****************************************************************************************************
    //preconditioners
    //****************************************************************************************************
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

    py::class_<PreconditionerType, PreconditionerType::Pointer>(m,"Preconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<PreconditionerType>)
    ;

    typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
    py::class_<DiagonalPreconditionerType, DiagonalPreconditionerType::Pointer, PreconditionerType>(m,"DiagonalPreconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<DiagonalPreconditionerType>)
    ;

    typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;
    py::class_<ILUPreconditionerType, ILUPreconditionerType::Pointer, PreconditionerType>(m,"ILUPreconditioner")
    .def(py::init< >() )
    .def("__str__", PrintObject<ILUPreconditionerType>)
    ;

    typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
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
    ;

    py::class_<ComplexLinearSolverType, ComplexLinearSolverType::Pointer>(m,"ComplexLinearSolver")
    .def(py::init< >() )
    .def("Initialize",&ComplexLinearSolverType::Initialize)
    .def("Solve",pointer_to_complex_solve)
    .def("Clear",&ComplexLinearSolverType::Clear)
    .def("__str__", PrintObject<ComplexLinearSolverType>)
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

    typedef Reorderer<SpaceType,  LocalSpaceType > ReordererType;
    typedef DirectSolver<SpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
    typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

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

}

}  // namespace Python.

} // Namespace Kratos
