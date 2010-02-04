/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-04-07 09:50:09 $
//   Revision:            $Revision: 1.3 $
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


//builder_and_solvers
#include "custom_strategies/builder_and_solvers/residualbased_elimination_quasiincompresible_builder_and_solver.h"

//linear solvers
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
			
		
			typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
//			typedef ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverType;

			//********************************************************************
			//********************************************************************
			//typedef ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, 2> //ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType2D;
			
			typedef ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, 2> ResidualBasedIncompressibleBuilderType2D;
			
						
			//class_< ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType2D, boost::noncopyable>
//("ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver2D", init< LinearSolverType::Pointer>() )
//
//

			class_< ResidualBasedIncompressibleBuilderType2D, bases< BuilderAndSolverType >, boost::noncopyable>
("ResidualBasedIncompressibleBuilder2D", init< LinearSolverType::Pointer>() )
			.def("AssembleLHS", &ResidualBasedIncompressibleBuilderType2D::AssembleLHS )
			.def("AssembleRHS", &ResidualBasedIncompressibleBuilderType2D::AssembleRHS )
			.def("BuildAndSolve", &ResidualBasedIncompressibleBuilderType2D::BuildAndSolve)
			.def("SetUpDofSet", &ResidualBasedIncompressibleBuilderType2D::SetUpDofSet)
			.def("SetUpSystem", &ResidualBasedIncompressibleBuilderType2D::SetUpSystem)
			.def("ResizeAndInitializeVectors", &ResidualBasedIncompressibleBuilderType2D::ResizeAndInitializeVectors)
			.def("Build", &ResidualBasedIncompressibleBuilderType2D::Build)
			.def("ConstructMatrixStructure", &ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure)
			.def("ConstructMatrixStructure_Mconsistent", &ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure_Mconsistent)
			.def("ConstructMatrixStructure_DivergenceMatrixD", &ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure_DivergenceMatrixD)
			.def("BuildAuxiliaries", &ResidualBasedIncompressibleBuilderType2D::BuildAuxiliaries)
			.def("AssembleMassMatrices", &ResidualBasedIncompressibleBuilderType2D::AssembleMassMatrices)
			.def("calc_GMinvD_prod", &ResidualBasedIncompressibleBuilderType2D::calc_GMinvD_prod)
			.def("CalculatePreconditionerDiagonalMatrix",  &ResidualBasedIncompressibleBuilderType2D::CalculatePreconditionerDiagonalMatrix)
			.def("calc_prod_precond_vec", &ResidualBasedIncompressibleBuilderType2D::calc_prod_precond_vec)   
			.def("ModifyForDirichlet", &ResidualBasedIncompressibleBuilderType2D::ModifyForDirichlet)   
 			.def("UpdatePressures", &ResidualBasedIncompressibleBuilderType2D::UpdatePressures)   
			.def("ReturnDx", &ResidualBasedIncompressibleBuilderType2D::ReturnDx)   
 			.def("UpdatePressuresNew", &ResidualBasedIncompressibleBuilderType2D::UpdatePressuresNew)   
 			.def("CalculateNodalPressureForce", &ResidualBasedIncompressibleBuilderType2D::CalculateNodalPressureForce )	
 			.def("ConvergenceCheck", &ResidualBasedIncompressibleBuilderType2D::ConvergenceCheck)   
			;

			typedef ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, 3> ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D;
			
						
			class_< ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D, bases< BuilderAndSolverType >, boost::noncopyable> ("ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver3D", init< LinearSolverType::Pointer>() )
			.def("AssembleLHS", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::AssembleLHS )
			.def("AssembleRHS", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::AssembleRHS )
			.def("BuildAndSolve", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::BuildAndSolve)
			.def("SetUpDofSet", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::SetUpDofSet)
			.def("SetUpSystem", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::SetUpSystem)
			.def("ResizeAndInitializeVectors", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ResizeAndInitializeVectors)
			.def("Build", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::Build)
			.def("ConstructMatrixStructure", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConstructMatrixStructure)
			.def("ConstructMatrixStructure_Mconsistent", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConstructMatrixStructure_Mconsistent)
			.def("ConstructMatrixStructure_DivergenceMatrixD", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConstructMatrixStructure_DivergenceMatrixD)
			.def("BuildAuxiliaries", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::BuildAuxiliaries)
			.def("AssembleMassMatrices", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::AssembleMassMatrices)
			.def("calc_GMinvD_prod", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::calc_GMinvD_prod)
			.def("CalculatePreconditionerDiagonalMatrix",  &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::CalculatePreconditionerDiagonalMatrix)
			.def("calc_prod_precond_vec", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::calc_prod_precond_vec)
			.def("ModifyForDirichlet", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ModifyForDirichlet)
			.def("UpdatePressures", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::UpdatePressures)
			.def("ReturnDx", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ReturnDx)	
			.def("UpdatePressuresNew", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::UpdatePressuresNew)	
			.def("CalculateNodalPressureForce", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::CalculateNodalPressureForce )	
			.def("ConvergenceCheck", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConvergenceCheck)	
			;
			//********************************************************************
			//********************************************************************
			
			
		}

	}  // namespace Python.

} // Namespace Kratos

