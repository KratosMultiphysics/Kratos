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
#include "custom_solvers/eigensystem_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomSolversToPython()
{
	using namespace boost::python;

	using complex_t = std::complex<double>;
	
	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	
    typedef UblasSpace<complex_t, ComplexCompressedMatrix, ComplexVector> ComplexSparseSpaceType;
    typedef UblasSpace<complex_t, ComplexMatrix, ComplexVector> ComplexLocalSpaceType;

	typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
	typedef DirectSolver<SparseSpaceType, LocalSpaceType> DirectSolverType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;

	// --- direct solvers

	using SparseLUSolver = EigenDirectSolver<SparseLU<double>, SparseSpaceType, LocalSpaceType>;
	class_<SparseLUSolver, bases<DirectSolverType>, boost::noncopyable>
		("SparseLUSolver", init<>())
		.def(init<Parameters>())
	;
	
	using ComplexSparseLUSolver = EigenDirectSolver<SparseLU<complex_t>, ComplexSparseSpaceType, ComplexLocalSpaceType>;
	class_<ComplexSparseLUSolver, bases<ComplexLinearSolverType>, boost::noncopyable>
		("ComplexSparseLUSolver", init<>())
		.def(init<Parameters>())
	;

	#if defined USE_EIGEN_MKL
	using PardisoLLTSolver = EigenDirectSolver<PardisoLLT<double>, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLLTSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLLTSolver", init<>())
		.def(init<Parameters>())
	;

	using PardisoLDLTSolver = EigenDirectSolver<PardisoLDLT<double>, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLDLTSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLDLTSolver", init<>())
		.def(init<Parameters>())
	;

	using PardisoLUSolver = EigenDirectSolver<PardisoLU<double>, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLUSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLUSolver", init<>())
		.def(init<Parameters>())
	;	
	#endif // defined USE_EIGEN_MKL

	using SparseQRSolver = EigenDirectSolver<SparseQR<double>, SparseSpaceType, LocalSpaceType>;
	class_<SparseQRSolver, bases<DirectSolverType>, boost::noncopyable>
		("SparseQRSolver", init<>())
		.def(init<Parameters>())
	;

	// --- eigensystem solver

	#if defined USE_EIGEN_MKL
	using EigensystemSolverType = EigensystemSolver<PardisoLDLT<double>, SparseSpaceType, LocalSpaceType>;
	#else  // defined USE_EIGEN_MKL
	using EigensystemSolverType = EigensystemSolver<SparseLU<double>, SparseSpaceType, LocalSpaceType>;
	#endif // defined USE_EIGEN_MKL
	class_<EigensystemSolverType, EigensystemSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable>
    	("EigensystemSolver", init<Parameters>())
    	.def("Solve", &EigensystemSolverType::Solve)
    	.def("GetEigenValue", &EigensystemSolverType::GetEigenValue)
	;
;
}

} // namespace Python

} // namespace Kratos
