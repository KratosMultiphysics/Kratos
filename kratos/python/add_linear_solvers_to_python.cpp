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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "includes/dof.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "linear_solvers/mixedup_linear_solver.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"
#include "linear_solvers/deflated_gmres_solver.h"



namespace Kratos
{

namespace Python
{
void  AddLinearSolversToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
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
    typedef DeflatedGMRESSolver<SpaceType,  LocalSpaceType> DeflatedGMRESSolverType;

    bool (LinearSolverType::*pointer_to_solve)(LinearSolverType::SparseMatrixType& rA, LinearSolverType::VectorType& rX, LinearSolverType::VectorType& rB) = &LinearSolverType::Solve;

    using namespace boost::python;

    //****************************************************************************************************
    //preconditioners
    //****************************************************************************************************
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

    class_<PreconditionerType, PreconditionerType::Pointer, boost::noncopyable>("Preconditioner")
    .def(self_ns::str(self))
    ;

    typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
    class_<DiagonalPreconditionerType, DiagonalPreconditionerType::Pointer, bases<PreconditionerType>, boost::noncopyable >("DiagonalPreconditioner")
    .def(self_ns::str(self))
    ;

    typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;
    class_<ILUPreconditionerType, ILUPreconditionerType::Pointer, bases<PreconditionerType>, boost::noncopyable >("ILUPreconditioner")
    .def(self_ns::str(self))
    ;

    typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
    class_<ILU0PreconditionerType, ILU0PreconditionerType::Pointer, bases<PreconditionerType>, boost::noncopyable >("ILU0Preconditioner")
    .def(self_ns::str(self))
    ;

    //****************************************************************************************************
    //linear solvers
    //****************************************************************************************************
    class_<LinearSolverType, LinearSolverType::Pointer, boost::noncopyable>("LinearSolver")
    .def("Initialize",&LinearSolverType::Initialize)
    .def("Solve",pointer_to_solve)
    .def("Clear",&LinearSolverType::Clear)
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<IterativeSolverType, IterativeSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable >("IterativeSolver")
    .def(self_ns::str(self))
    ;

    class_<CGSolverType, CGSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("CGSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(init<Parameters,  PreconditionerType::Pointer>())
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<BICGSTABSolverType, BICGSTABSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("BICGSTABSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(self_ns::str(self))
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(init<Parameters,  PreconditionerType::Pointer>())
    .def("SetTolerance",&BICGSTABSolverType::SetTolerance)
    ;

    class_<TFQMRSolverType, TFQMRSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("TFQMRSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(self_ns::str(self))
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(init<Parameters,  PreconditionerType::Pointer>())
    ;

    class_<ScalingSolverType, ScalingSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable >("ScalingSolver")
    .def(init<LinearSolverType::Pointer, bool >())
    ;

    class_<PowerIterationEigenvalueSolverType, PowerIterationEigenvalueSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable >("PowerIterationEigenvalueSolver")
    .def(init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
    ;

    typedef Reorderer<SpaceType,  LocalSpaceType > ReordererType;
    typedef DirectSolver<SpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
    typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

    class_<ReordererType, ReordererType::Pointer, boost::noncopyable >("Reorderer")
    .def( init< >() )
    .def(self_ns::str(self))
    .def( "Initialize",&ReordererType::Initialize)
    .def( "Reorder",&ReordererType::Reorder)
    .def( "InverseReorder",&ReordererType::InverseReorder)
    ;

    class_<DirectSolverType, DirectSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable >("DirectSolver")
    .def( init< >() )
    .def(init<Parameters>())
    .def(self_ns::str(self))
    ;

    class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >("SkylineLUFactorizationSolver")
    .def(init< >())
    .def(init<Parameters>())
    .def(self_ns::str(self))
    ;

    class_<DeflatedCGSolverType, DeflatedCGSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("DeflatedCGSolver")
    .def(init<double,bool,int>())
    .def(init<double, unsigned int,bool,int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer,bool,int>())
    .def(init<Parameters>())
// 		  .def(init<double, unsigned int,  PreconditionerType::Pointer, ModelPart::Pointer>())
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<MixedUPLinearSolverType, MixedUPLinearSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("MixedUPLinearSolver",init<LinearSolverType::Pointer, LinearSolverType::Pointer ,double, unsigned int, unsigned int >())
    .def(init<Parameters,LinearSolverType::Pointer, LinearSolverType::Pointer >())
    .def(self_ns::str(self))
    ;

    class_<DeflatedGMRESSolverType, DeflatedGMRESSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("DeflatedGMRESSolver",init<LinearSolverType::Pointer ,double, unsigned int, unsigned int, unsigned int >())
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

