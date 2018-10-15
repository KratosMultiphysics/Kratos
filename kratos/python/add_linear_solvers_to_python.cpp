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
#include "python/add_equation_systems_to_python.h"
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
#include "linear_solvers/mixedup_linear_solver.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"
#include "linear_solvers/power_iteration_highest_eigenvalue_solver.h"
#include "linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h"
#include "linear_solvers/deflated_gmres_solver.h"

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
    using TDirectSolverType = DirectSolver<TSpaceType<TDataType>, TLocalSpaceType<TDataType>>;

void  AddLinearSolversToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;
    typedef CGSolver<SpaceType,  LocalSpaceType> CGSolverType;
    typedef DeflatedCGSolver<SpaceType,  LocalSpaceType> DeflatedCGSolverType;
    typedef MixedUPLinearSolver<SpaceType,  LocalSpaceType> MixedUPLinearSolverType;
    typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
    typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;
    typedef ScalingSolver<SpaceType,  LocalSpaceType> ScalingSolverType;
    typedef PowerIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> PowerIterationEigenvalueSolverType;
    typedef PowerIterationHighestEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> PowerIterationHighestEigenvalueSolverType;
    typedef RayleighQuotientIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> RayleighQuotientIterationEigenvalueSolverType;
    typedef DeflatedGMRESSolver<SpaceType,  LocalSpaceType> DeflatedGMRESSolverType;

    typedef TSpaceType<std::complex<double>> ComplexSparseSpaceType;
    typedef TLocalSpaceType<std::complex<double>> ComplexDenseSpaceType;
    typedef TLinearSolverType<std::complex<double>> ComplexLinearSolverType;
    typedef TDirectSolverType<std::complex<double>> ComplexDirectSolverType;
    typedef SkylineLUCustomScalarSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexSkylineLUSolverType;

    bool (LinearSolverType::*pointer_to_solve)(LinearSolverType::SparseMatrixType& rA, LinearSolverType::VectorType& rX, LinearSolverType::VectorType& rB) = &LinearSolverType::Solve;
    void (LinearSolverType::*pointer_to_solve_eigen)(LinearSolverType::SparseMatrixType& rK, LinearSolverType::SparseMatrixType& rM,LinearSolverType::DenseVectorType& Eigenvalues, LinearSolverType::DenseMatrixType& Eigenvectors) = &LinearSolverType::Solve;
    bool (ComplexLinearSolverType::*pointer_to_complex_solve)(ComplexLinearSolverType::SparseMatrixType& rA, ComplexLinearSolverType::VectorType& rX, ComplexLinearSolverType::VectorType& rB) = &ComplexLinearSolverType::Solve;

    using namespace pybind11;

    //****************************************************************************************************
    //preconditioners
    //****************************************************************************************************
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

    class_<PreconditionerType, PreconditionerType::Pointer>(m,"Preconditioner")
    .def( init< >() )
    .def("__str__", KRATOS_DEF_PYTHON_STR(PreconditionerType))
    ;

    typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
    class_<DiagonalPreconditionerType, DiagonalPreconditionerType::Pointer, PreconditionerType>(m,"DiagonalPreconditioner")
    .def( init< >() )
    .def("__str__", KRATOS_DEF_PYTHON_STR(DiagonalPreconditionerType))
    ;

    typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;
    class_<ILUPreconditionerType, ILUPreconditionerType::Pointer, PreconditionerType>(m,"ILUPreconditioner")
    .def( init< >() )
    .def("__str__", KRATOS_DEF_PYTHON_STR(ILUPreconditionerType))
    ;

    typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
    class_<ILU0PreconditionerType, ILU0PreconditionerType::Pointer, PreconditionerType>(m,"ILU0Preconditioner")
    .def( init< >() )
    .def("__str__", KRATOS_DEF_PYTHON_STR(ILU0PreconditionerType))
    ;

    //****************************************************************************************************
    //linear solvers
    //****************************************************************************************************
    class_<LinearSolverType, LinearSolverType::Pointer>(m,"LinearSolver")
    .def( init< >() )
    .def("Initialize",&LinearSolverType::Initialize)
    .def("Solve",pointer_to_solve)
    .def("Solve",pointer_to_solve_eigen)
    .def("Clear",&LinearSolverType::Clear)
    .def("__str__", KRATOS_DEF_PYTHON_STR(LinearSolverType))
    ;

    class_<ComplexLinearSolverType, ComplexLinearSolverType::Pointer>(m,"ComplexLinearSolver")
    .def( init< >() )
    .def("Initialize",&ComplexLinearSolverType::Initialize)
    .def("Solve",pointer_to_complex_solve)
    .def("Clear",&ComplexLinearSolverType::Clear)
    .def("__str__", KRATOS_DEF_PYTHON_STR(ComplexLinearSolverType))
    ;

    class_<IterativeSolverType, IterativeSolverType::Pointer, LinearSolverType>(m,"IterativeSolver")
    .def( init< >() )
    .def("__str__", KRATOS_DEF_PYTHON_STR(IterativeSolverType))
    ;

    class_<CGSolverType, CGSolverType::Pointer,IterativeSolverType>(m,"CGSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(init<Parameters,  PreconditionerType::Pointer>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(CGSolverType))
    ;

    class_<BICGSTABSolverType, BICGSTABSolverType::Pointer,IterativeSolverType>(m,"BICGSTABSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(BICGSTABSolverType))
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(init<Parameters,  PreconditionerType::Pointer>())
    .def("SetTolerance",&BICGSTABSolverType::SetTolerance)
    ;

    class_<TFQMRSolverType, TFQMRSolverType::Pointer,IterativeSolverType>(m,"TFQMRSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(TFQMRSolverType))
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(init<Parameters,  PreconditionerType::Pointer>())
    ;

    class_<ScalingSolverType, ScalingSolverType::Pointer, LinearSolverType>(m,"ScalingSolver")
    .def(init<LinearSolverType::Pointer, bool >())
    ;

    class_<PowerIterationEigenvalueSolverType, PowerIterationEigenvalueSolverType::Pointer, LinearSolverType>(m,"PowerIterationEigenvalueSolver")
    .def(init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
    .def(init<Parameters, LinearSolverType::Pointer>())
    ;

    class_<PowerIterationHighestEigenvalueSolverType, PowerIterationHighestEigenvalueSolverType::Pointer, LinearSolverType>(m,"PowerIterationHighestEigenvalueSolver")
    .def(init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
    .def(init<Parameters, LinearSolverType::Pointer>())
    ;

    class_<RayleighQuotientIterationEigenvalueSolverType, RayleighQuotientIterationEigenvalueSolverType::Pointer, LinearSolverType>(m,"RayleighQuotientIterationEigenvalueSolver")
    .def(init<double, unsigned int, unsigned int, LinearSolverType::Pointer, double>())
    .def(init<Parameters, LinearSolverType::Pointer>())
    ;

    typedef Reorderer<SpaceType,  LocalSpaceType > ReordererType;
    typedef DirectSolver<SpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
    typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

    class_<ReordererType, ReordererType::Pointer>(m,"Reorderer")
    .def( init< >() )
    .def("__str__", KRATOS_DEF_PYTHON_STR(ReordererType))
    .def( "Initialize",&ReordererType::Initialize)
    .def( "Reorder",&ReordererType::Reorder)
    .def( "InverseReorder",&ReordererType::InverseReorder)
    ;

    class_<DirectSolverType, DirectSolverType::Pointer, LinearSolverType>(m,"DirectSolver")
    .def( init< >() )
    .def(init<Parameters>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(DirectSolverType))
    ;

    class_<ComplexDirectSolverType, ComplexDirectSolverType::Pointer, ComplexLinearSolverType>(m,"ComplexDirectSolver")
    .def( init< >() )
    .def(init<Parameters>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(ComplexDirectSolverType))
    ;

    class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, DirectSolverType>(m,"SkylineLUFactorizationSolver")
    .def(init< >())
    .def(init<Parameters>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(SkylineLUFactorizationSolverType))
    ;

    class_<ComplexSkylineLUSolverType, typename ComplexSkylineLUSolverType::Pointer, ComplexDirectSolverType>(m,"ComplexSkylineLUSolver")
    .def(init< >())
    .def(init<Parameters&>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(ComplexSkylineLUSolverType))
    ;

    class_<DeflatedCGSolverType, DeflatedCGSolverType::Pointer,IterativeSolverType>(m,"DeflatedCGSolver")
    .def(init<double,bool,int>())
    .def(init<double, unsigned int,bool,int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer,bool,int>())
    .def(init<Parameters>())
// 		  .def(init<double, unsigned int,  PreconditionerType::Pointer, ModelPart*>())
    //.def("",&LinearSolverType::)
    .def("__str__", KRATOS_DEF_PYTHON_STR(DeflatedCGSolverType))
    ;

    class_<MixedUPLinearSolverType, MixedUPLinearSolverType::Pointer,IterativeSolverType>(m,"MixedUPLinearSolver")
    .def(init<LinearSolverType::Pointer, LinearSolverType::Pointer ,double, unsigned int, unsigned int >())
    .def(init<Parameters,LinearSolverType::Pointer, LinearSolverType::Pointer >())
    .def("__str__", KRATOS_DEF_PYTHON_STR(MixedUPLinearSolverType))
    ;

    class_<DeflatedGMRESSolverType, DeflatedGMRESSolverType::Pointer,IterativeSolverType>(m,"DeflatedGMRESSolver")
    .def(init<LinearSolverType::Pointer ,double, unsigned int, unsigned int, unsigned int >())
    .def("__str__", KRATOS_DEF_PYTHON_STR(DeflatedGMRESSolverType))
    ;

}

}  // namespace Python.

} // Namespace Kratos

