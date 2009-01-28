/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: pooyan $
*   Date:                $Date: 2008-10-27 09:24:43 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/


#if !defined(KRATOS_PETSC_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_PETSC_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER


/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"


namespace Kratos
{

	/**@name Kratos Globals */
	/*@{ */


	/*@} */
	/**@name Type Definitions */       
	/*@{ */

	/*@} */


	/**@name  Enum's */       
	/*@{ */


	/*@} */
	/**@name  Functions */       
	/*@{ */



	/*@} */
	/**@name Kratos Classes */
	/*@{ */

	/** Short class definition.

	Detail class definition.

	Current class provides an implementation for standard builder and solving operations.

	the RHS is constituted by the unbalanced loads (residual)

	Degrees of freedom are reordered putting the restrained degrees of freedom at 
	the end of the system ordered in reverse order with respect to the DofSet.

	Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
	this information.

	Calculation of the reactions involves a cost very similiar to the calculation of the total residual

	\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

	\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


	\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

	\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

	\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

	\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


	*/
	template<class TSparseSpace,
	class TDenseSpace  //= DenseSpace<double>,
	  //	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
	class PetscResidualBasedEliminationBuilderAndSolver 
		: public BuilderAndSolver< TSparseSpace,TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		//typedef boost::shared_ptr< PetscResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION( PetscResidualBasedEliminationBuilderAndSolver );


		typedef BuilderAndSolver<TSparseSpace,TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> > BaseType;

		typedef TSparseSpace SparseSpaceType;

		typedef typename BaseType::TSchemeType TSchemeType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
 
		typedef typename BaseType::NodesArrayType NodesArrayType;
		typedef typename BaseType::ElementsArrayType ElementsArrayType;
		typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

		typedef typename BaseType::ElementsContainerType ElementsContainerType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		PetscResidualBasedEliminationBuilderAndSolver(
			)
		  : BaseType(typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer(new LinearSolver<TSparseSpace,TDenseSpace>))
		{

			/* 			std::cout << "using the standard builder and solver " << std::endl; */

		}


		/** Destructor.
		*/
		virtual ~PetscResidualBasedEliminationBuilderAndSolver(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

		//**************************************************************************
		//**************************************************************************
		void Build(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& b)
		{
			KRATOS_TRY
				if(!pScheme)
					KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

			//getting the elements from the model
			ElementsArrayType& pElements = r_model_part.Elements();

			//getting the array of the conditions
			ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			//resetting to zero the vector of reactions
// 			TSparseSpace::SetToZero(BaseType::mReactionsVector);

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
			LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
// 			for(typename NodesArrayType::iterator i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); ++i_node)
// 			  std::cout << "rank " << rank << " : [Node #" << i_node->Id() << ": coordinates " <<  i_node->X() << ", " <<   i_node->Y() << ", " <<   i_node->Z() << ", reference coordinates" <<  i_node->X0() << ", " <<   i_node->Y0() << ", " <<   i_node->Z0() << "]" << std::endl;
			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
				//calculate elemental contribution
				pScheme->CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);
//KRATOS_WATCH(LHS_Contribution);
//KRATOS_WATCH(RHS_Contribution);

// std::cout << "rank " << rank << " : [element #" << (*it)->Id() << ": Eq. Id " <<  EquationId[0] << ", " <<  EquationId[1] << ", " <<  EquationId[2] << ", " <<  EquationId[3] << ", " <<  EquationId[4] << ", " <<  EquationId[5] << "]" << std::endl;
// std::cout << "rank " << rank << " : [element #" << (*it)->Id() << " LHS_Contribution : " << LHS_Contribution << "]" << std::endl;				//bbb = GetTickCount();
// std::cout << "rank " << rank << " : [element #" << (*it)->Id() << " RHS_Contribution : " << RHS_Contribution << "]" << std::endl;				//bbb = GetTickCount();
				//assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);
				AssembleRHS(b,RHS_Contribution,EquationId);
	//ccc += GetTickCount() - bbb;

				// clean local elemental memory
				pScheme->CleanMemory(*it);
			}
//double EndTime = GetTickCount();
//
//std::cout << "total time " << EndTime - StartTime << std::endl;
//std::cout << "writing in the system matrix " << ccc << std::endl;
//std::cout << "calculating the elemental contrib " << ddd << std::endl;
			LHS_Contribution.resize(0,0,false);
			RHS_Contribution.resize(0,false);

			// assemble all conditions
			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			{
				//calculate elemental contribution
				pScheme->Condition_CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

// std::cout << "rank " << rank << " : [condition #" << (*it)->Id() << ": Eq. Id " <<  EquationId[0] << ", " <<  EquationId[1] << ", " <<  EquationId[2] << ", " <<  EquationId[3] << "]" << std::endl;
// std::cout << "rank " << rank << " : [condition #" << (*it)->Id() << " LHS_Contribution : " << LHS_Contribution << "]" << std::endl;				//bbb = GetTickCount();
// std::cout << "rank " << rank << " : [condition #" << (*it)->Id() << " RHS_Contribution : " << RHS_Contribution << "]" << std::endl;				//bbb = GetTickCount();
				//assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);
				AssembleRHS(b,RHS_Contribution,EquationId);
			}
//for(int i = 0; i<b.size(); i++)
//std::cout << b[i] << std::endl;

//std::cout << "stop" << std::endl;
// 			MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
// 			MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
// 			VecAssemblyBegin(b);
// 			VecAssemblyEnd(b);

			KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************
		void BuildLHS(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A)
		{
			KRATOS_TRY

				//getting the elements from the model
				ElementsArrayType& pElements = r_model_part.Elements();

			//getting the array of the conditions
			ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			//resetting to zero the vector of reactions
// 			TSparseSpace::SetToZero(BaseType::mReactionsVector);

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
				//calculate elemental contribution
				pScheme->Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);

				// clean local elemental memory
				pScheme->CleanMemory(*it);
			}

			LHS_Contribution.resize(0,0,false);

			// assemble all conditions
			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			{
				//calculate elemental contribution
				pScheme->Condition_Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);
			}

 // 			MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//  			MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
			KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************
		void BuildLHS_CompleteOnFreeRows(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A)
		{
			KRATOS_TRY

				//getting the elements from the model
				ElementsArrayType& pElements = r_model_part.Elements();

			//getting the array of the conditions
			ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			//resetting to zero the vector of reactions
// 			TSparseSpace::SetToZero(BaseType::mReactionsVector);

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
				//calculate elemental contribution
				pScheme->Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				AssembleLHS_CompleteOnFreeRows(A,LHS_Contribution,EquationId);

				// clean local elemental memory
				pScheme->CleanMemory(*it);
			}

			LHS_Contribution.resize(0,0,false);
			// assemble all conditions
			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			{
				//calculate elemental contribution
				pScheme->Condition_Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				AssembleLHS_CompleteOnFreeRows(A,LHS_Contribution,EquationId);
			}


			KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************
		void SystemSolve(
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
		{
			KRATOS_TRY

		    KRATOS_WATCH("MatAssemblyBegin");
 			MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		    KRATOS_WATCH("MatAssemblyEnd");
 			MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
		    KRATOS_WATCH("VecAssemblyBegin");
 			VecAssemblyBegin(b);
		    KRATOS_WATCH("VecAssemblyEnd");
 			VecAssemblyEnd(b);
		    KRATOS_WATCH("VecAssemblyBegin");
 			VecAssemblyBegin(Dx);
		    KRATOS_WATCH("VecAssemblyEnd");
 			VecAssemblyEnd(Dx);
		    KRATOS_WATCH("End");
//  			VecAssemblyBegin(BaseType::mReactionsVector);
//  			VecAssemblyEnd(BaseType::mReactionsVector);
// 			VecView(b,PETSC_VIEWER_STDOUT_SELF);
// 			VecView(Dx,PETSC_VIEWER_STDOUT_SELF);
// 			MatView(A,PETSC_VIEWER_STDOUT_SELF);
			double norm_b;
			if(TSparseSpace::Size(b) != 0)
				norm_b = TSparseSpace::TwoNorm(b);
			else
				norm_b = 0.00;

			if(norm_b != 0.00)
			  {
			    KSP ksp;
			    PetscInt ierr;
			    /* 
			       Create linear solver context
			    */

// 			    KRATOS_WATCH("KSPCreate");
			    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);PetscErrorChecking(ierr);
			    
			    /* 
			       Set operators. Here the matrix that defines the linear system
			       also serves as the preconditioning matrix.
			    */
// 			    KRATOS_WATCH("KSPSetOperators");
			    ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);PetscErrorChecking(ierr);
			    
			    /* 
			       Set linear solver defaults for this problem (optional).
			       - By extracting the KSP and PC contexts from the KSP context,
			       we can then directly call any KSP and PC routines to set
			       various options.
			       - The following two statements are optional; all of these
			       parameters could alternatively be specified at runtime via
			       KSPSetFromOptions().  All of these defaults can be
			       overridden at runtime, as indicated below.
			    */
// 			    KRATOS_WATCH("KSPSetTolerances");
			    ierr = KSPSetTolerances(ksp,1.e-9,1.e-10,PETSC_DEFAULT,
						    20000);PetscErrorChecking(ierr);
			    
			    /* 
			       Set runtime options, e.g.,
			       -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
			       These options will override those specified above as long as
			       KSPSetFromOptions() is called _after_ any other customization
			       routines.
			    */
// 			    KRATOS_WATCH("KSPSetFromOptions");
			    ierr = KSPSetFromOptions(ksp);PetscErrorChecking(ierr);

			    KRATOS_WATCH("KSPSetType"); 
 			    KSPSetType(ksp, "bcgs"); 
			    
// 			    const char* ksp_name = new char[255];
// 			    KSPGetType(ksp,&ksp_name);

// 			    std::cout << "Solver Type is : " << ksp_name << std::endl;

// 			    delete [] ksp_name;

			    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
			       Solve the linear system
			       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
			    
			    KRATOS_WATCH("KSPSolve");
			    ierr = KSPSolve(ksp,b,Dx);PetscErrorChecking(ierr);
			    KRATOS_WATCH("KSPSolveFinished");

 			    //KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);
			    

// 			VecView(b,PETSC_VIEWER_STDOUT_SELF);
// 			VecView(Dx,PETSC_VIEWER_STDOUT_SELF);
// 			MatView(A,PETSC_VIEWER_STDOUT_SELF);
			    
			  }
			else
// 			  VecView(b,PETSC_VIEWER_STDOUT_SELF);
				TSparseSpace::SetToZero(Dx);

			//prints informations about the current time
			if (this->GetEchoLevel()>1)
			{
				std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
			}

			KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************
		void BuildAndSolve(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY

			boost::timer building_time;

			Build(pScheme,r_model_part,A,b);

			if(BaseType::GetEchoLevel()>0)
			{
				std::cout << "Building Time : " << building_time.elapsed() << std::endl;
			}
			

//			ApplyPointLoads(pScheme,r_model_part,b);

			//does nothing...dirichlet conditions are naturally dealt with in defining the residual
			ApplyDirichletConditions(pScheme,r_model_part,A,Dx,b);

			if (BaseType::GetEchoLevel()== 3)
			{
				std::cout << "before the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
			}

			boost::timer solve_time;

			SystemSolve(A,Dx,b);

			if(BaseType::GetEchoLevel()>0)
			{
				std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
			}
			if (BaseType::GetEchoLevel()== 3)
			{
				std::cout << "after the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
			}
			
			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void BuildRHSAndSolve(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY

				BuildRHS(pScheme,r_model_part,b);
			SystemSolve(A,Dx,b);

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void BuildRHS(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemVectorType& b)
		{
			KRATOS_TRY		

			  //Getting the Elements
			  ElementsArrayType& pElements = r_model_part.Elements();

			  //getting the array of the conditions
			  ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			  
			  //resetting to zero the vector of reactions
			  // 			TSparseSpace::SetToZero(BaseType::mReactionsVector);
			  
			  //contributions to the system
			  LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
			  LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
			  
			  //vector containing the localization in the system of the different
			  //terms
			  Element::EquationIdVectorType EquationId;
			  
			  // assemble all elements
			  for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			    {
			      //calculate elemental Right Hand Side Contribution
			      pScheme->Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
			      
			      //assemble the elemental contribution
			      AssembleRHS(b,RHS_Contribution,EquationId);
			    }
			  
			  LHS_Contribution.resize(0,0,false);
			  RHS_Contribution.resize(0,false);
			  
			  // assemble all conditions
			  for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			    {
			      //calculate elemental contribution
			      pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
			      
			      //assemble the elemental contribution
			      AssembleRHS(b,RHS_Contribution,EquationId);
			    }
			  
//  			VecAssemblyBegin(b);
//  			VecAssemblyEnd(b);
			  KRATOS_CATCH("")
			    
		}				
		//**************************************************************************
		//**************************************************************************
		void SetUpDofSet(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part
			) 
		{
			KRATOS_TRY

			  //Gets the array of elements from the modeler
 			  ElementsArrayType& pElements = r_model_part.Elements(ModelPart::Kratos_Local);

			  Element::DofsVectorType ElementalDofList;

			  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			  DofsArrayType Doftemp;
			  BaseType::mDofSet = DofsArrayType();

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);

			  for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			    {
			      // gets list of Dof involved on every element
				pScheme->GetElementalDofList(*it,ElementalDofList,CurrentProcessInfo);
				
				for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
				{
// 				  if(( (*i)->Id() > 340) && ( (*i)->Id() < 360))
// 				    std::cout << rank << " : adding dof " << (*i)->Id() << " " << **i << std::endl;
// 					if(( (*i)->Id() > 349) && ( (*i)->Id() < 356))
// 					  std::cout << rank << " : before dof push_back " << (*i)->Id() << " with partition index : " <<(*i)->GetSolutionStepValue(PARTITION_INDEX) << " in " << *i << std::endl;
					Doftemp.push_back(*i);
					//mDofSet.push_back(*i);

// 					if(( (*i)->Id() > 349) && ( (*i)->Id() < 356))
// 					  std::cout << rank << " : dof push_back " << (*i)->Id() << " with partition index : " <<(*i)->GetSolutionStepValue(PARTITION_INDEX) << " in " << *i << std::endl;
				}
			}
	
			//taking in account conditions
			ConditionsArrayType& pConditions = r_model_part.Conditions();
			for (typename ConditionsArrayType::ptr_iterator it=pConditions.ptr_begin(); it!=pConditions.ptr_end(); ++it)
			{
				// gets list of Dof involved on every element
				pScheme->GetConditionDofList(*it,ElementalDofList,CurrentProcessInfo);

				for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
				{
					//mDofSet.push_back(*i);
					Doftemp.push_back(*i);
				}
			}


//  			for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
//  				  if(( dof_iterator->Id() > 349) && ( dof_iterator->Id() < 356))
//  				    std::cout << rank << " : doftemp before unique " << dof_iterator->Id() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
			Doftemp.Unique();

//  			for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
//  				  if(( dof_iterator->Id() > 349) && ( dof_iterator->Id() < 356))
//  				    std::cout << rank << " : doftemp after unique " << dof_iterator->Id() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
			BaseType::mDofSet = Doftemp;

 // 			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
//  				  if(( dof_iterator->Id() > 349) && ( dof_iterator->Id() < 356))
//  				    std::cout << rank << " : mDofSet : " << dof_iterator->Id() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
			//throws an execption if there are no Degrees of freedom involved in the analysis
			if (BaseType::mDofSet.size()==0)
				KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

			// Fixing the ghost degrees of freedom
// 			NodesArrayType& rNodes = r_model_part.Nodes(ModelPart::Kratos_Ghost);
// 			for(typename NodesArrayType::iterator i_node = rNodes.ptr_begin(); i_node != rNodes.end(); ++i_node)
// 			  for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
// 			    i_dof->FixDof();

			BaseType::mDofSetIsInitialized = true;

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void SetUpSystem(
			ModelPart& r_model_part
			)
		{ 
			// Set equation id for degrees of freedom
			// the free degrees of freedom are positioned at the beginning of the system, 
			// while the fixed one are at the end (in opposite order).
			//
			// that means that if the EquationId is greater than "mEquationSystemSize"
			// the pointed degree of freedom is restrained
			//
			int free_size = 0;
			int fixed_size = 0;

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);

			// Calculating number of fixed and free dofs
			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
			  if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
			    {
				if (dof_iterator->IsFixed())
				  fixed_size++;
				else
				  free_size++;
			    }

			// Calculating the total size and required offset
			int free_offset;
			int fixed_offset;
			int global_size;

			// The correspounding offset by the sum of the sizes in thread with inferior rank
			MPI_Scan(&free_size, &free_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			// The total size by the sum of all size in all threads
			MPI_Allreduce(&free_size, &global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			
			// The correspounding fixed offset by the sum of the fixed sizes in thread with inferior rank
			MPI_Scan(&fixed_size, &fixed_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			// finding the offset for the begining of the partition
			free_offset -= free_size;

			fixed_offset += global_size - fixed_size;
			
			std::cout << rank << " : local size = " << BaseType::mDofSet.size() << std::endl; 
			std::cout << rank << " : free_id = " << free_size << std::endl; 
			std::cout << rank << " : fixed_size = " << fixed_size << std::endl; 
			std::cout << rank << " : free_offset = " << free_offset << std::endl; 
			std::cout << rank << " : fixed offset = " << fixed_offset << std::endl;

			// Now setting the equation id with .
			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
			  if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
			    {
				if (dof_iterator->IsFixed())
					dof_iterator->SetEquationId(fixed_offset++);
				else
					dof_iterator->SetEquationId(free_offset++);
 				std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
			    }


			BaseType::mEquationSystemSize = global_size;
			mLocalSystemSize = free_size;
			std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl; 
			std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
			std::cout << rank << " : free_offset = " << free_offset << std::endl; 
			std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;


// 			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
// 			  {
// 			  if(( dof_iterator->Id() == 347) || (dof_iterator->EquationId() == 687))
//  				    std::cout << rank << " : mDofSet : " << dof_iterator->Id() << " with equation id : " <<dof_iterator->EquationId() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
// 			  }
			
			UpdateGhostDofs(r_model_part);

			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
			  {
			  dof_iterator->GetSolutionStepValue(TEMPERATURE) = dof_iterator->EquationId();
// 			  if(( dof_iterator->Id() == 347) || (dof_iterator->EquationId() == 687))
//  				    std::cout << rank << " : Updated mDofSet : " << dof_iterator->Id() << " with equation id : " <<dof_iterator->EquationId() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
			  }
		}
 
	  void UpdateGhostDofs(ModelPart& rThisModelPart)
		{
		  int rank;
		  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		  std::cout << rank << " : Strarting UpdateGhostDofs...." << std::endl;

		  int source=rank;
		  int destination=0;

		  vector<int>& neighbours_indices = rThisModelPart[NEIGHBOURS_INDICES];

		  std::cout << rank << " starting domain loop " << std::endl;
		  for(int i_domain = 0 ; i_domain <  neighbours_indices.size() ; i_domain++)
		    if((destination = neighbours_indices[i_domain]) >= 0)
		      {
			std::cout << rank << " domian #" << i_domain << std::endl;
			unsigned int send_buffer_size = 0;
			unsigned int receive_buffer_size = 0;

			std::cout << rank;
			KRATOS_WATCH(destination);
			// Calculating send and received buffer size
			// The interface meshes are stored after all, local and ghost meshes
			NodesArrayType& r_interface_nodes = rThisModelPart.Nodes(ModelPart::Kratos_Ownership_Size + i_domain);
			NodesArrayType& r_ghost_nodes = rThisModelPart.Nodes(ModelPart::Kratos_Ghost);

			std::cout << rank << " : 2...." << std::endl;
			for(typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
			  send_buffer_size += i_node->GetDofs().size();

			std::cout << rank << " : 3...." << std::endl;
			for(typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
			  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
			    {
			    receive_buffer_size += i_node->GetDofs().size();
			    for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
			      {
				std::cout << rank << " : receving : " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
			      }
			    }
			unsigned int position = 0;
			int* send_buffer = new int[send_buffer_size];
			int* receive_buffer = new int[receive_buffer_size];


			// Filling the buffer
			std::cout << rank << " :  Filling the buffer...." << std::endl;
			for(ModelPart::NodeIterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
 			  for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
			    {
			      send_buffer[position++] = i_dof->EquationId();
 			    std::cout << rank << " : sending equation id : " << i_dof->EquationId() << " for " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
			    }
// 		    {
// 		      std::cout << rank << " : sending DISPLACEMENT_X eq id : " << i_node->GetDof(DISPLACEMENT_X).EquationId() << " for node " << i_node->Id() << std::endl;
// 		      send_buffer[position++] = i_node->GetDof(DISPLACEMENT_X).EquationId();
// 		      std::cout << rank << " : sending DISPLACEMENT_Y eq id : " << i_node->GetDof(DISPLACEMENT_Y).EquationId() << " for node " << i_node->Id() << std::endl;
// 		      send_buffer[position++] = i_node->GetDof(DISPLACEMENT_Y).EquationId();
// 		      std::cout << rank << " : sending DISPLACEMENT_Z eq id : " << i_node->GetDof(DISPLACEMENT_Z).EquationId() << " for node " << i_node->Id() << std::endl;
// 		      send_buffer[position++] = i_node->GetDof(DISPLACEMENT_Z).EquationId();
// 		    }


			MPI_Status status;

////////////////////////// test
/* 			std::cout << rank ; */
/* 			KRATOS_WATCH(source); */
/* 			std::cout << rank ; */
/* 			KRATOS_WATCH(destination); */
/* 			int* test_send_buffer = new int[1]; */
/* 			int* test_recv_buffer = new int[1]; */
/* 			test_send_buffer[0] = rank; */
/* MPI_Sendrecv (test_send_buffer, 1, MPI_INT, destination, 12, test_recv_buffer, 1, MPI_INT, destination, 12, MPI_COMM_WORLD, &status); */
/* 			std::cout << rank ; */
/* 			KRATOS_WATCH(test_send_buffer[0]); */
/* 			std::cout << rank ; */
/* 			KRATOS_WATCH(test_recv_buffer[0]); */
/* 			delete [] test_send_buffer; */
/* 			delete [] test_recv_buffer; */











			if(position > send_buffer_size)
			  std::cout << rank << " Error in estimating send buffer size...." << std::endl;
		  

			int send_tag = 1;//i_domain;
			int receive_tag = 1;//i_domain;
			int rc;

			std::cout << rank << " : Strarting Send and receive...." << std::endl;
			std::cout << rank ;
			KRATOS_WATCH(source);
			std::cout << rank ;
			KRATOS_WATCH(destination);
			std::cout << rank ;
			KRATOS_WATCH(send_buffer_size);
			std::cout << rank ;
			KRATOS_WATCH(receive_buffer_size);
			std::cout << rank ;
			KRATOS_WATCH(send_tag);
			std::cout << rank ;
			KRATOS_WATCH(receive_tag);


			MPI_Sendrecv (send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
 				 MPI_COMM_WORLD, &status);

			std::cout << rank << " : Send and receive Finished" << std::endl;
			
			// Updating nodes
			position = 0;
			for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin(ModelPart::Kratos_Ghost) ; 
			    i_node != rThisModelPart.NodesEnd(ModelPart::Kratos_Ghost) ; i_node++)
			  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
			    for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
			  {
			    i_dof->SetEquationId(receive_buffer[position++]);
			    // 			    std::cout << rank << " : receiving equation id  : " << i_dof->EquationId() <<  " for " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
			  }
			// 		    {
			// 		      i_node->GetDof(DISPLACEMENT_X).SetEquationId(receive_buffer[position++]);
			// 		      i_node->GetDof(DISPLACEMENT_Y).SetEquationId(receive_buffer[position++]);
			// 		      i_node->GetDof(DISPLACEMENT_Z).SetEquationId(receive_buffer[position++]);
			// 		    }
			
			if(position > receive_buffer_size)
			  std::cout << rank << " Error in estimating receive buffer size...." << std::endl;
			
			delete [] send_buffer;
			delete [] receive_buffer;
		      }
		  
		}
		//**************************************************************************
		//**************************************************************************
		void ResizeAndInitializeVectors(
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b,
			ElementsArrayType& rElements,
			ConditionsArrayType& rConditions,
			ProcessInfo& CurrentProcessInfo
			)
		{
			KRATOS_TRY

			  PetscInt ierr;

			  //if needed resize the vector for the calculation of reactions
// 			  if(BaseType::mCalculateReactionsFlag == true)
// 			    {
// 			      unsigned int ReactionsVectorSize = BaseType::mDofSet.size()-BaseType::mEquationSystemSize;
// 			      if(BaseType::mReactionsVector.size() != ReactionsVectorSize)
// 				{
// 				  ierr = VecCreate(PETSC_COMM_WORLD,&BaseType::mReactionsVector);PetscErrorChecking(ierr);
// 				  ierr = VecSetSizes(BaseType::mReactionsVector,PETSC_DECIDE,ReactionsVectorSize);PetscErrorChecking(ierr);
// 				}
// 				BaseType::mReactionsVector.resize(ReactionsVectorSize,false);
// 			    }
			  
				//resizing the system vectors and matrix
//  				if (TSparseSpace::Size1(A) == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
//   				{
				  ierr = MatCreate(PETSC_COMM_WORLD,&A);PetscErrorChecking(ierr);
				  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,BaseType::mEquationSystemSize,BaseType::mEquationSystemSize);PetscErrorChecking(ierr);
				  ConstructMatrixStructure(A,rElements,rConditions,CurrentProcessInfo);
//   				}
//  				else
//  				{
// 				  if((TSparseSpace::Size1(A) != BaseType::mEquationSystemSize) || (TSparseSpace::Size2(A) != BaseType::mEquationSystemSize))
// 				    {
// 				      KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
// 				      ierr = MatCreate(PETSC_COMM_WORLD,&A);PetscErrorChecking(ierr);
// 				      ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,BaseType::mEquationSystemSize,BaseType::mEquationSystemSize);PetscErrorChecking(ierr);
// 				      ConstructMatrixStructure(A,rElements,rConditions,CurrentProcessInfo);
// 				    }
// 				}
				  

			  
				ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE, BaseType::mEquationSystemSize,&Dx);PetscErrorChecking(ierr);
// 				if(TSparseSpace::Size(Dx) != BaseType::mEquationSystemSize)
// 				{
				  ierr = VecSetType(Dx, VECMPI);PetscErrorChecking(ierr);
				  ierr = VecSetFromOptions(Dx);PetscErrorChecking(ierr);
				  ierr = VecDuplicate(Dx,&b);PetscErrorChecking(ierr);  
// 				}



				
// 					Dx.resize(BaseType::mEquationSystemSize,false);
// 				if(b.size() != BaseType::mEquationSystemSize)
// 					b.resize(BaseType::mEquationSystemSize,false);

// 				//


// 				//if needed resize the vector for the calculation of reactions
// 				if(BaseType::mCalculateReactionsFlag == true)
// 				{
// 					unsigned int ReactionsVectorSize = BaseType::mDofSet.size()-BaseType::mEquationSystemSize;
// 					if(BaseType::mReactionsVector.size() != ReactionsVectorSize)
// 						BaseType::mReactionsVector.resize(ReactionsVectorSize,false);
// 				}

			  
				KRATOS_CATCH("")

		}



		//**************************************************************************
		//**************************************************************************
		void InitializeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void FinalizeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
		}


		//**************************************************************************
		//**************************************************************************
		void CalculateReactions(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
// 			//refresh RHS to have the correct reactions
// 			BuildRHS(pScheme,r_model_part,b);

// 			int i;
// 			int systemsize = BaseType::mDofSet.size() - SparseSpaceType::Size(BaseType::mReactionsVector);

// 			typename DofsArrayType::ptr_iterator it2;
// 			//std::set<Dof::Pointer,ComparePDof>::iterator it2;

// 			//updating variables
// 			//for (it2=mDofSet.begin();it2 != mDofSet.end(); ++it2)
// 			for (it2=BaseType::mDofSet.ptr_begin();it2 != BaseType::mDofSet.ptr_end(); ++it2)
// 			{
// 				if ( (*it2)->IsFixed()  )
// 				{
// 					i=(*it2)->EquationId();
// 					i-=systemsize; 

// 					 VecGetValues(BaseType::mReactionsVector, 1, &i, &(*it2)->GetSolutionStepReactionValue());
// 				}
// 			}
		}

		//**************************************************************************
		//**************************************************************************
		void ApplyDirichletConditions(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{}

		//**************************************************************************
		//**************************************************************************
		void ApplyPointLoads(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemVectorType& b)
		{}

		/** 
		this function is intended to be called at the end of the solution step to clean up memory
		storage not needed 
		*/
		void Clear() 
		{
			this->mDofSet = DofsArrayType();
			this->mReactionsVector = TSystemVectorType();

			if (this->GetEchoLevel()>0)
			{

				KRATOS_WATCH("PetscResidualBasedEliminationBuilderAndSolver Clear Function called");
			}
		}




		/*@} */
		/**@name Operations */
		/*@{ */


		/*@} */  
		/**@name Access */
		/*@{ */


		/*@} */
		/**@name Inquiry */
		/*@{ */


		/*@} */      
		/**@name Friends */
		/*@{ */


		/*@} */

	protected:
		/**@name Protected static Member Variables */
		/*@{ */


		/*@} */
		/**@name Protected member Variables */
		/*@{ */


		/*@} */
		/**@name Protected Operators*/
		/*@{ */
		//**************************************************************************
		virtual void ConstructMatrixStructure(
			TSystemMatrixType& A,
			ElementsContainerType& rElements,
			ConditionsArrayType& rConditions,
			ProcessInfo& CurrentProcessInfo)
		{

		  std::size_t equation_size = SparseSpaceType::Size1(A);
		  std::vector<std::vector<std::size_t> > indices(equation_size);
		  //				std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

		  Element::EquationIdVectorType ids(3,0);
		  for(typename ElementsContainerType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; i_element++)
		    {
		      (i_element)->EquationIdVector(ids, CurrentProcessInfo);
		      
		      for(std::size_t i = 0 ; i < ids.size() ; i++)
			if(ids[i] < equation_size)
			  {
			    std::vector<std::size_t>& row_indices = indices[ids[i]];
			    for(std::size_t j = 0 ; j < ids.size() ; j++)
			      if(ids[j] < equation_size)
				{
				  AddUnique(row_indices,ids[j]);
				  //indices[ids[i]].push_back(ids[j]);
				}
			  }
		      
		    }
		  
		  for(typename ConditionsArrayType::iterator i_condition = rConditions.begin() ; i_condition != rConditions.end() ; i_condition++)
		    {
		      (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
		      for(std::size_t i = 0 ; i < ids.size() ; i++)
			if(ids[i] < equation_size)
			  {
			    std::vector<std::size_t>& row_indices = indices[ids[i]];
			    for(std::size_t j = 0 ; j < ids.size() ; j++)
			      if(ids[j] < equation_size)
				{
				  AddUnique(row_indices,ids[j]);
				  //	indices[ids[i]].push_back(ids[j]);
				}
			  }
		    }
		  
		  //allocating the memory needed 
// 		  int data_size = 0;
// 		  for(std::size_t i = 0 ; i < indices.size() ; i++)
// 		    {
// 		      data_size += indices[i].size();
// 		    }
// 		  A.reserve(data_size,false);
		  
// 		  //filling with zero the matrix (creating the structure)
// 		  for(std::size_t i = 0 ; i < indices.size() ; i++)
// 		    {
// 		      std::vector<std::size_t>& row_indices = indices[i];
// 		      std::sort(row_indices.begin(), row_indices.end());
		      
// 		      for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
// 			{
// 			  A.push_back(i,*it,0.00);
// 			  //					A()(i,*it) = 0.00;
// 			}
// 		      //row_indices = std::vector<std::size_t>(); 
// 		      row_indices.clear(); 
// 		    }
		  
			  PetscInt ierr;
			  PetscInt* row_nonzeros = new PetscInt[BaseType::mEquationSystemSize];
			  for(int i = 0 ; i < BaseType::mEquationSystemSize ; i++)
			    row_nonzeros[i] =  indices[i].size();

std::cout <<  "Resizing the Matrix A ..." << std::endl;
			  ierr = MatSetType(A, MATMPIAIJ);PetscErrorChecking(ierr);
			  //ierr = MatSetType(A, MATBAIJ);PetscErrorChecking(ierr);
			  ierr = MatSetFromOptions(A);PetscErrorChecking(ierr);
			  /* ierr = MatMPIAIJSetPreallocation(A,PETSC_DEFAULT , row_nonzeros,PETSC_DEFAULT , row_nonzeros);PetscErrorChecking(ierr); */
			  ierr = MatSeqAIJSetPreallocation(A,PETSC_DEFAULT , row_nonzeros);PetscErrorChecking(ierr);
			  delete row_nonzeros;
			  

		}

		//**************************************************************************
		void AssembleLHS(
			TSystemMatrixType& A,
			LocalSystemMatrixType& LHS_Contribution,
			Element::EquationIdVectorType& EquationId
			)
		{
		  //		  TSparseSpace::AssembleLHS(A, LHS_Contribution, EquationId);
			unsigned int local_size = LHS_Contribution.size1();

			for (unsigned int i_local=0; i_local<local_size; i_local++)
			{
				unsigned int i_global=EquationId[i_local];
				if ( i_global < BaseType::mEquationSystemSize )
				{
					for (unsigned int j_local=0; j_local<local_size; j_local++)
					{
						unsigned int j_global=EquationId[j_local];
						if ( j_global < BaseType::mEquationSystemSize )
						  MatSetValue(A, i_global, j_global, LHS_Contribution(i_local,j_local), ADD_VALUES);
					}
				}
			}
		}



		//**************************************************************************
		void AssembleRHS(
			TSystemVectorType& b,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId
			)
		{
			unsigned int local_size = RHS_Contribution.size();

// 			if (BaseType::mCalculateReactionsFlag==false) //if we don't need to calculate reactions
// 			{
				for (unsigned int i_local=0; i_local<local_size; i_local++)
				{
					unsigned int i_global=EquationId[i_local];
					if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs
					{	// ASSEMBLING THE SYSTEM VECTOR
					  VecSetValue(b, i_global, RHS_Contribution[i_local], ADD_VALUES);
						
					}
				}
/* 			  TSparseSpace::AssembleRHS(b, RHS_Contribution, EquationId); */
// 			}
// 			else //when the calculation of reactions is needed
// 			{
/* 				for (unsigned int i_local=0; i_local<local_size; i_local++) */
/* 				{ */
/* 					unsigned int i_global=EquationId[i_local]; */
/* 					if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs */
/* 					{	// ASSEMBLING THE SYSTEM VECTOR */
/* 						b[i_global] += RHS_Contribution[i_local]; */
/* 					} */
/* 					else //on "fixed" DOFs */
/* 					{	// Assembling the Vector of REACTIONS */
/* 						BaseType::mReactionsVector[i_global-BaseType::mEquationSystemSize] -= RHS_Contribution[i_local]; */
/* 					} */
/* 				} */
// 			}
		}



		/*@} */
		/**@name Protected Operations*/
		/*@{ */


		/*@} */
		/**@name Protected  Access */
		/*@{ */


		/*@} */     
		/**@name Protected Inquiry */
		/*@{ */


		/*@} */   
		/**@name Protected LifeCycle */  
		/*@{ */



		/*@} */    

	private:
		/**@name Static Member Variables */
		/*@{ */


		/*@} */
		/**@name Member Variables */
		/*@{ */

		unsigned int mLocalSystemSize;

		/*@} */
		/**@name Private Operators*/
		/*@{ */


		/*@} */
		/**@name Private Operations*/
		/*@{ */


		//**************************************************************************
		void AssembleLHS_CompleteOnFreeRows(
			TSystemMatrixType& A,
			LocalSystemMatrixType& LHS_Contribution,
			Element::EquationIdVectorType& EquationId
			)
		{
		  KRATOS_ERROR(std::logic_error, "This method is not implemented for petsc", "");
// 			unsigned int local_size = LHS_Contribution.size1();
// 			for (unsigned int i_local=0; i_local<local_size; i_local++)
// 			{
// 				unsigned int i_global=EquationId[i_local];
// 				if ( i_global < BaseType::mEquationSystemSize )
// 				{
// 					for (unsigned int j_local=0; j_local<local_size; j_local++)
// 					{
// 						int j_global=EquationId[j_local];

// 						A(i_global,j_global) += LHS_Contribution(i_local,j_local);
// 					}
// 				}
// 			}
		}		


		//******************************************************************************************
		//******************************************************************************************
		inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
		{
			std::vector<std::size_t>::iterator i = v.begin();
			std::vector<std::size_t>::iterator endit = v.end();
			while ( i != endit && (*i) != candidate)
			{
				i++;
			}
			if( i == endit )
			{
				v.push_back(candidate);
			}

		}

	        int PetscErrorChecking(int ErrorIndex)
	        {
		  CHKERRQ(ErrorIndex);
		  return ErrorIndex;
		}

		/*@} */
		/**@name Private  Access */
		/*@{ */


		/*@} */     
		/**@name Private Inquiry */
		/*@{ */


		/*@} */   
		/**@name Un accessible methods */
		/*@{ */


		/*@} */   

	}; /* Class PetscResidualBasedEliminationBuilderAndSolver */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_PETSC_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

