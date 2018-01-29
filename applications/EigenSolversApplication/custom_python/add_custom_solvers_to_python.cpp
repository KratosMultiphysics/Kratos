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
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_solvers_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/sparse_eigensystem_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomSolversToPython()
{
	using namespace boost::python;

	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
	typedef DirectSolver<SparseSpaceType, LocalSpaceType> DirectSolverType;

	using SparseLUSolver = EigenDirectSolver<SparseLU, SparseSpaceType, LocalSpaceType>;
	class_<SparseLUSolver, bases<DirectSolverType>, boost::noncopyable>
		("SparseLUSolver", init<>())
		.def(init<Parameters>());

	#if defined EIGEN_USE_MKL_ALL
	using PardisoLLTSolver = EigenDirectSolver<PardisoLLT, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLLTSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLLTSolver", init<>())
		.def(init<Parameters>());

	using PardisoLDLTSolver = EigenDirectSolver<PardisoLDLT, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLDLTSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLDLTSolver", init<>())
		.def(init<Parameters>());

	using PardisoLUSolver = EigenDirectSolver<PardisoLU, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLUSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLUSolver", init<>())
		.def(init<Parameters>());
	#endif

	using SparseEigensystemSolverType = SparseEigensystemSolver<SparseSpaceType, LocalSpaceType>;
	class_<SparseEigensystemSolverType, SparseEigensystemSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable>
    	("SparseEigensystemSolver", init<Parameters>())
    	.def("Solve", &SparseEigensystemSolverType::Solve)
    	.def("GetEigenValue", &SparseEigensystemSolverType::GetEigenValue)
	;
;
}

} // namespace Python

} // namespace Kratos
