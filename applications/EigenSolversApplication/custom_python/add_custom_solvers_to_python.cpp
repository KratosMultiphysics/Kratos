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
#include "linear_solvers/linear_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigensystem_solver.h"

namespace Kratos
{

namespace Python
{

template <typename SolverType>
void register_solver(pybind11::module& m, const std::string& name)
{
	using namespace pybind11;

	using Space = SpaceType<typename SolverType::TScalar>;

	using Base = DirectSolver<typename Space::Global, typename Space::Local>;

	class_<EigenDirectSolver<SolverType>, Base >
		(m, name.c_str())
		.def(init<>())
		.def(init<Parameters>())
	;
}

template <typename SolverType>
void register_eigensystem_solver(pybind11::module& m, const std::string& name = SolverType::Name)
{
	using namespace pybind11;

	using Space = SpaceType<typename SolverType::TScalar>;

	using Base = LinearSolver<typename Space::Global, typename Space::Local>;

	class_<EigensystemSolver<SolverType>, Base >
		(m, name.c_str())
		.def(init<Parameters>())
	;
}

void AddCustomSolversToPython(pybind11::module& m)
{
	using namespace pybind11;

	using complex = std::complex<double>;

	// --- direct solvers

	register_solver<SparseLU<double>>(m, "SparseLUSolver");
	register_solver<SparseLU<complex>>(m, "ComplexSparseLUSolver");

	register_solver<SparseQR<double>>(m, "SparseQRSolver");
	register_solver<SparseQR<complex>>(m, "ComplexSparseQRSolver");

	#if defined USE_EIGEN_MKL
	// The commented complex solvers need to be tested first
	register_solver<PardisoLLT<double>>(m, "PardisoLLTSolver");
	// register_solver<PardisoLLT<complex>>(m, "ComplexPardisoLLTSolver");

	register_solver<PardisoLDLT<double>>(m, "PardisoLDLTSolver");
	// register_solver<PardisoLDLT<complex>>(m, "ComplexPardisoLDLTSolver");

	register_solver<PardisoLU<double>>(m, "PardisoLUSolver");
	register_solver<PardisoLU<complex>>(m, "ComplexPardisoLUSolver");
	#endif // defined USE_EIGEN_MKL

	// --- eigensystem solver

	#if !defined USE_EIGEN_MKL
	register_eigensystem_solver<SparseLU<double>>(m, "EigensystemSolver");
	#else  // !defined USE_EIGEN_MKL
	register_eigensystem_solver<PardisoLDLT<double>>(m, "EigensystemSolver");
	#endif // !defined USE_EIGEN_MKL
}

} // namespace Python

} // namespace Kratos
