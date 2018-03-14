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
#include "linear_solvers/linear_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigensystem_solver.h"

namespace Kratos
{

namespace Python
{

template <typename SolverType>
void register_solver(const std::string& name)
{
	using namespace boost::python;

	using Space = typename SpaceType<typename SolverType::TScalar>;

	using Base = DirectSolver<typename Space::Global, typename Space::Local>;
	
	class_<EigenDirectSolver<SolverType>, bases<Base>, boost::noncopyable>
		(name.c_str(), init<>())
		.def(init<Parameters>())
	;
}

template <typename SolverType>
void register_eigensystem_solver(const std::string& name = SolverType::Name)
{
	using namespace boost::python;

	using Space = typename SpaceType<typename SolverType::TScalar>;

	using Base = LinearSolver<typename Space::Global, typename Space::Local>;
	
	using EigenSolver = EigensystemSolver<SolverType>;

	class_<EigenSolver, bases<Base>, boost::noncopyable>
		(name.c_str(), init<Parameters>())
    	.def("Solve", &EigenSolver::Solve)
    	.def("GetEigenValue", &EigenSolver::GetEigenValue)
	;
}

void AddCustomSolversToPython()
{
	using namespace boost::python;

	using complex = std::complex<double>;

	// --- direct solvers

	register_solver<SparseLU<double>>("SparseLUSolver");
	register_solver<SparseLU<complex>>("ComplexSparseLUSolver");
	
	register_solver<SparseQR<double>>("SparseQRSolver");
	register_solver<SparseQR<complex>>("ComplexSparseQRSolver");

	#if defined USE_EIGEN_MKL
	register_solver<PardisoLLT<double>>("PardisoLLTSolver");
	register_solver<PardisoLLT<complex>>("ComplexPardisoLLTSolver");

	register_solver<PardisoLDLT<double>>("PardisoLDLTSolver");
	register_solver<PardisoLDLT<complex>>("ComplexPardisoLDLTSolver");
	
	register_solver<PardisoLU<double>>("PardisoLUSolver");
	register_solver<PardisoLU<complex>>("ComplexPardisoLUSolver");
	#endif // defined USE_EIGEN_MKL

	// --- eigensystem solver

	#if !defined USE_EIGEN_MKL
	register_eigensystem_solver<SparseLU<double>>("EigensystemSolver");
	#else  // !defined USE_EIGEN_MKL
	register_eigensystem_solver<PardisoLDLT<double>>("EigensystemSolver");
	#endif // !defined USE_EIGEN_MKL
}

} // namespace Python

} // namespace Kratos
