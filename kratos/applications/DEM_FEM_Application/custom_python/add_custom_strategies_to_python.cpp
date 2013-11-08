/*
==============================================================================
KratosStructuralApplication 
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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:49:07 $
//   Revision:            $Revision: 1.20 $
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
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/DEM_FEM_Explicit_Solver.h"
#include "custom_strategies/FEM_Explicit_Forward_Solver.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"

#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    namespace Python
    {		
        using namespace boost::python;

        void  AddCustomStrategiesToPython()
        {
            typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
            typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

            typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
            typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >
                    BaseSolvingStrategyType;
                    

            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

         
            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;
         
                    
            typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>
                    BuilderAndSolverType;

	
	typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
	//typedef ResidualCriteria < SparseSpaceType, LocalSpaceType >::Pointer TResidual;
	//typedef DisplacementCriteria < SparseSpaceType, LocalSpaceType>::Pointer TDisplacement;

	
// 					;
           //********************************************************************
           //********************************************************************
		   
  
           class_< DemFemExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
                   (
                 "DemFemExplicitSolverStrategy", init< ModelPart&, ModelPart&, double ,int,double, bool, BaseSolvingStrategyType&, BaseSolvingStrategyType&>())
                 .def("Initialize", &DemFemExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Initialize)
				 .def("Solve", &DemFemExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Solve)
              ;  
			  
			class_< FemExplicitForwardStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
                   (
                 "FemExplicitForwardStrategy", init< ModelPart&, int, int ,double, bool, double, bool, bool>() )
                 .def("Initialize", &FemExplicitForwardStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Initialize)
				 .def("Solve", &FemExplicitForwardStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Solve)
              ; 
        }
    }  // namespace Python.

} // Namespace Kratos

