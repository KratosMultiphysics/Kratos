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

	using Base = DirectSolver<typename SolverType::TGlobalSpace, typename SolverType::TLocalSpace>;

	using EigenDirectSolverType = EigenDirectSolver<SolverType>;

	class_<EigenDirectSolverType, typename EigenDirectSolverType::Pointer, Base >
		(m, name.c_str())
		.def(init<>())
		.def(init<Parameters>())
	;
}

template <typename SolverType>
void register_eigensystem_solver(pybind11::module& m, const std::string& name = SolverType::Name)
{
	using namespace pybind11;

	using Base = LinearSolver<typename SolverType::TGlobalSpace, typename SolverType::TLocalSpace>;

	using EigenSystemSolverType = EigensystemSolver<SolverType>;

	class_<EigenSystemSolverType, typename EigenSystemSolverType::Pointer, Base >
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

	{
		using Type = boost::numeric::ublas::vector<double>;

		class_<Type>(m, "UblasVector")
			.def(init<>())
			.def(init<int>())
			.def("__getitem__",
				[](const Type& self, std::size_t index) -> double {
					return self[index];
				})
			.def("__setitem__",
				[](Type& self, std::size_t index, double value) -> void {
					self[index] = value;
				})
		;
	}
}

} // namespace Python

} // namespace Kratos
