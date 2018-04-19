/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_solvers_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigensystem_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomSolversToPython(pybind11::module& m)
{
	using namespace pybind11;

	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
	typedef DirectSolver<SparseSpaceType, LocalSpaceType> DirectSolverType;

	// --- direct solvers

	using SparseLUSolver = EigenDirectSolver<SparseLU, SparseSpaceType, LocalSpaceType>;
	class_<SparseLUSolver, SparseLUSolver::Pointer, DirectSolverType>
		(m, "SparseLUSolver")
		.def(init<>())
		.def(init<Parameters>())
	;

	#if defined USE_EIGEN_MKL
	using PardisoLLTSolver = EigenDirectSolver<PardisoLLT, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLLTSolver, PardisoLLTSolver::Pointer, DirectSolverType>
		(m, "PardisoLLTSolver")
		.def(init<>())
		.def(init<Parameters>())
	;

	using PardisoLDLTSolver = EigenDirectSolver<PardisoLDLT, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLDLTSolver, PardisoLDLTSolver::Pointer, DirectSolverType>
		(m, "PardisoLDLTSolver")
		.def(init<>())
		.def(init<Parameters>())
	;

	using PardisoLUSolver = EigenDirectSolver<PardisoLU, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLUSolver, PardisoLUSolver::Pointer, DirectSolverType>
		(m, "PardisoLUSolver")
		.def(init<>())
		.def(init<Parameters>())
	;
	#endif // defined USE_EIGEN_MKL

	using SparseQRSolver = EigenDirectSolver<SparseQR, SparseSpaceType, LocalSpaceType>;
	class_<SparseQRSolver, SparseQRSolver::Pointer, DirectSolverType>
		(m, "SparseQRSolver")
		.def(init<>())
		.def(init<Parameters>())
	;

	// --- eigensystem solver

	#if defined USE_EIGEN_MKL
	using EigensystemSolverType = EigensystemSolver<PardisoLDLT, SparseSpaceType, LocalSpaceType>;
	#else  // defined USE_EIGEN_MKL
	using EigensystemSolverType = EigensystemSolver<SparseLU, SparseSpaceType, LocalSpaceType>;
	#endif // defined USE_EIGEN_MKL
	class_<EigensystemSolverType, EigensystemSolverType::Pointer, LinearSolverType>
    	(m, "EigensystemSolver")
		.def(init<Parameters>())
    	.def("Solve", &EigensystemSolverType::Solve)
    	.def("GetEigenValue", &EigensystemSolverType::GetEigenValue)
	;
;
}

} // namespace Python

} // namespace Kratos
