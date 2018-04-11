/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// External includes
#include <boost/timer.hpp>

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

// Project includes
#include "includes/define.h"
#include "custom_strategies/bfeccStrategy.h"
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m) {

	typedef UblasSpace<double,CompressedMatrix,Vector> 		SparseSpaceType;
	typedef UblasSpace<double,Matrix,Vector>           		LocalSpaceType;

	typedef LinearSolver<SparseSpaceType,LocalSpaceType>	LinearSolverType;
    typedef Scheme<SparseSpaceType,LocalSpaceType>          BaseSchemeType;

	typedef SolvingStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>
    BaseSolvingStrategyType;

	typedef BfeccSolverStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>
		BfeccSolverStrategyType;

	using namespace pybind11;

    class_< BfeccSolverStrategyType, BfeccSolverStrategyType::Pointer, BaseSolvingStrategyType> (m, "BfeccSolverStrategy")
        .def(init<ModelPart&>())
		.def("WriteResults", &BfeccSolverStrategyType::WriteResults)
		.def_property("Dt",
			&BfeccSolverStrategyType::GetDt,
			&BfeccSolverStrategyType::SetDt)
		.def_property("NumCells",
			&BfeccSolverStrategyType::GetNumCells,
			&BfeccSolverStrategyType::SetNumCells)
		.def_property("BorderWidth",
			&BfeccSolverStrategyType::GetBorderWidth,
			&BfeccSolverStrategyType::SetBorderWidth)
	;
}

}  // namespace Python.
} // Namespace Kratos
