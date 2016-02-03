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


// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies


//linear solvers


//convergence criteria


//schemes
#include "solving_strategies/schemes/scheme.h"
#include "../DEM_application/DEM_application.h"
#include "../DEM_application/custom_strategies/schemes/dem_integration_scheme.h"
#include "custom_strategies/schemes/terminal_velocity_scheme.h"

//parallel strategies
//#include "custom_strategies/strategies/mpi_explicit_solver_strategy.h" //MPI CARLOS descomentar aixo

//parallel schemes
//#include "custom_strategies/schemes/mpi_foward_euler_scheme.h" //MPI CARLOS descomentar aixo

//builder_and_solvers


namespace Kratos
{

	namespace Python
	{		
		using namespace boost::python;

		void  AddCustomStrategiesToPython()
		{
            class_< TerminalVelocityScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>
            ("TerminalVelocityScheme", init<>());
		  
           //MPI CARLOS decomentar de aqui....
		  /*
           typedef MpiExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType > MpiExplicitSolverStrategyType;  
           class_< MpiExplicitSolverStrategyType, bases< BaseSolvingStrategyType >,  boost::noncopyable>
           (
           "MpiExplicitSolverStrategy", init< ModelPart&, ModelPart&, int, double, double, double, double, double, double, bool, bool, bool, IntegrationScheme::Pointer>())
                   .def("Initialize", &MpiExplicitSolverStrategyType::Initialized)
                   .def("InitialCriticalTime", &MpiExplicitSolverStrategyType::InitialCriticalTime)
           ;
		  */
		  // MPI CARLOS.... a aqui
	   
		}

	}  // namespace Python.

} // Namespace Kratos
