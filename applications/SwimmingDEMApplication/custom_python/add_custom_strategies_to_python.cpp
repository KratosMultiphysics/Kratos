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
//   Last modified by:    $Author:  G. Casas (gcasas@cimne.upc.edu)$
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//

// System includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

// External includes
#include "spaces/ublas_space.h"

//strategies
#include "../DEMApplication/custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/strategies/adams_bashforth_strategy.h"
#include "custom_strategies/strategies/residualbased_derivative_recovery_strategy.h"

// convergence criteria
#include "custom_strategies/convergence_criteria/vel_pr_criteria.h"

//schemes
#include "../DEMApplication/custom_strategies/schemes/dem_integration_scheme.h"
#include "solving_strategies/schemes/scheme.h"
#include "add_custom_utilities_to_python.h"
#include "custom_strategies/schemes/hybrid_bashforth_scheme.h"
#include "custom_strategies/schemes/terminal_velocity_scheme.h"
#include "custom_strategies/schemes/symplectic_euler_old_velocity_scheme.h"

//parallel strategies
//parallel schemes
//builder_and_solvers
//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
	namespace Python
	{
		namespace py = pybind11;

		void  AddCustomStrategiesToPython(pybind11::module& m)
		{

            py::class_< HybridBashforthScheme, typename HybridBashforthScheme::Pointer, SymplecticEulerScheme>
            (m,  "HybridBashforthScheme").def(py::init<>());

            py::class_< TerminalVelocityScheme, typename TerminalVelocityScheme::Pointer, HybridBashforthScheme>
            (m,  "TerminalVelocityScheme").def(py::init<>());

            py::class_< SymplecticEulerOldVelocityScheme, typename SymplecticEulerOldVelocityScheme::Pointer, SymplecticEulerScheme>
            (m,  "SymplecticEulerOldVelocityScheme").def(py::init<>());

            py::class_< AdamsBashforthStrategy, typename AdamsBashforthStrategy::Pointer, ExplicitSolverStrategy>
            (m,  "AdamsBashforthStrategy")
            .def(py::init<ExplicitSolverSettings&,
                      double,
                      double,
                      double,
                      int,
                      ParticleCreatorDestructor::Pointer,
                      DEM_FEM_Search::Pointer,
                      SpatialSearch::Pointer,
                      Parameters,
                      const bool>());

            typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
            typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
            typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
            //typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
            //typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ::Pointer TConvergenceCriteriaPointer;
            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
            //typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

            py::class_< ResidualBasedDerivativeRecoveryStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                    typename ResidualBasedDerivativeRecoveryStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                    ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > >
            (m,  "ResidualBasedDerivativeRecoveryStrategy")
            .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool,  bool>())
            .def("GetResidualNorm", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetResidualNorm)
            .def("SetBuilderAndSolver", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetBuilderAndSolver)
            ;
		}

	}  // namespace Python.

} // Namespace Kratos
