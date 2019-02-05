// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//         -        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//         -        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//                  in the documentation and/or other materials provided with the distribution.
//         -        All advertising materials mentioning features or use of this software must display the following acknowledgement:
//                         This product includes Kratos Multi-Physics technology.
//         -        Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/hexahedra_newton_raphson_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
	using namespace pybind11;

	void  AddCustomStrategiesToPython(pybind11::module& m)
	{
		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

		// Base types
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
		//typedef LinearSolverType::Pointer LinearSolverPointer;
		typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
		typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
		//typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
		//typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
		//typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;
		//typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;
		typedef HexahedraNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > HexahedraNewtonRaphsonStrategyType;

		// class_<HexahedraNewtonRaphsonStrategyType, typename HexahedraNewtonRaphsonStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType>(m, "HexahedraNewtonRaphsonStrategy")
		// 	.def(init<ModelPart& , BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, BuilderAndSolverType, int, bool, bool, bool, bool>());
		class_< HexahedraNewtonRaphsonStrategyType,
				typename HexahedraNewtonRaphsonStrategyType::Pointer,
				BaseSolvingStrategyType  >  (m, "HexahedraNewtonRaphsonStrategy")
				.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
				;
	}

}  // namespace Python.

} // Namespace Kratos
