// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
//#include "linear_solvers/superlu_solver.h"
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

    class_<PreconditionerType, PreconditionerType::Pointer>("Preconditioner")
    .def(self_ns::str(self))
    ;

    typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
    class_<DiagonalPreconditionerType, DiagonalPreconditionerType::Pointer, bases<PreconditionerType> >("DiagonalPreconditioner")
    .def(self_ns::str(self))
    ;

    typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;
    class_<ILUPreconditionerType, ILUPreconditionerType::Pointer, bases<PreconditionerType> >("ILUPreconditioner")
    .def(self_ns::str(self))
    ;

    typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
    class_<ILU0PreconditionerType, ILU0PreconditionerType::Pointer, bases<PreconditionerType> >("ILU0Preconditioner")
    .def(self_ns::str(self))
    ;

    //****************************************************************************************************
    //linear solvers
    //****************************************************************************************************
    class_<LinearSolverType, LinearSolverType::Pointer>("LinearSolver")
    .def("Initialize",&LinearSolverType::Initialize)
    .def("Solve",pointer_to_solve)
    .def("Clear",&LinearSolverType::Clear)
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<IterativeSolverType, IterativeSolverType::Pointer, bases<LinearSolverType> >("IterativeSolver")
    .def(self_ns::str(self))
    ;

    class_<CGSolverType, CGSolverType::Pointer, bases<IterativeSolverType> >("CGSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<BICGSTABSolverType, BICGSTABSolverType::Pointer, bases<IterativeSolverType> >("BICGSTABSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(self_ns::str(self))
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def("SetTolerance",&BICGSTABSolverType::SetTolerance)
    ;

    class_<TFQMRSolverType, TFQMRSolverType::Pointer, bases<IterativeSolverType> >("TFQMRSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(self_ns::str(self))
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    ;

    class_<ScalingSolverType, ScalingSolverType::Pointer, bases<LinearSolverType> >("ScalingSolver")
    .def(init<LinearSolverType::Pointer, bool >())
    ;

    class_<PowerIterationEigenvalueSolverType, PowerIterationEigenvalueSolverType::Pointer, bases<LinearSolverType> >("PowerIterationEigenvalueSolver")
    .def(init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
    ;

    typedef Reorderer<SpaceType,  LocalSpaceType > ReordererType;
    typedef DirectSolver<SpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
    typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

    class_<ReordererType, ReordererType::Pointer >("Reorderer")
    .def( init< >() )
    .def(self_ns::str(self))
    .def( "Initialize",&ReordererType::Initialize)
    .def( "Reorder",&ReordererType::Reorder)
    .def( "InverseReorder",&ReordererType::InverseReorder)
    ;

    class_<DirectSolverType, DirectSolverType::Pointer, bases<LinearSolverType> >("DirectSolver")
    .def( init< >() )
    .def(self_ns::str(self))
    ;

    class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, bases<DirectSolverType> >("SkylineLUFactorizationSolver")
    .def(init< >())
    .def(self_ns::str(self))
    ;

    class_<DeflatedCGSolverType, DeflatedCGSolverType::Pointer, bases<IterativeSolverType> >("DeflatedCGSolver")
    .def(init<double,bool,int>())
    .def(init<double, unsigned int,bool,int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer,bool,int>())
// 		  .def(init<double, unsigned int,  PreconditionerType::Pointer, ModelPart::Pointer>())
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<MixedUPLinearSolverType, MixedUPLinearSolverType::Pointer, bases<IterativeSolverType> >("MixedUPLinearSolver",init<LinearSolverType::Pointer, LinearSolverType::Pointer ,double, unsigned int, unsigned int >())
    .def(self_ns::str(self))
    ;

    class_<DeflatedGMRESSolverType, DeflatedGMRESSolverType::Pointer, bases<IterativeSolverType> >("DeflatedGMRESSolver",init<LinearSolverType::Pointer ,double, unsigned int, unsigned int, unsigned int >())
    .def(self_ns::str(self))
    ;

//	typedef SuperLUSolver<SparseSpaceType, LocalSpaceType> SuperLUSolverType;
//	class_<SuperLUSolverType, bases<DirectSolverType>, boost::noncopyable >
//		( "SuperLUSolver",
//			init<>() )
}

}  // namespace Python.

} // Namespace Kratos

