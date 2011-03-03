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
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2009-01-13 15:39:56 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/


#if !defined(KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_COMPONENTWISE_SLIP )
#define  KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_COMPONENTWISE_SLIP


/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp" 


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "Epetra_MpiComm.h"

//trilinos includes
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
// #include "epetra_test_err.h"


//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"


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
	class TDenseSpace,
	class TLinearSolver,
	class TVariableType
	>
	class TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit 
		: public BuilderAndSolver< TSparseSpace,TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		//typedef boost::shared_ptr< TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION( TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit );


		typedef BuilderAndSolver<TSparseSpace,TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> > BaseType;

		typedef TSparseSpace SparseSpaceType;

		typedef typename BaseType::TSchemeType TSchemeType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
		typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
 
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
		TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit(
			Epetra_MpiComm& Comm,
			int guess_row_size,
			typename TLinearSolver::Pointer pNewLinearSystemSolver,
			unsigned int dim, 
			TVariableType const& Var_x,
			TVariableType const& Var_y, 
			TVariableType const& Var_z)
			: BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
			, mrComm(Comm),mguess_row_size(guess_row_size), mdim(dim),mrVar_x(Var_x),mrVar_y(Var_y),mrVar_z(Var_z)
		{
			

			/* 			std::cout << "using the standard builder and solver " << std::endl; */

		}
// 		TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit(
// 			)
// 		  : BaseType(typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer(new LinearSolver<TSparseSpace,TDenseSpace>))
// 		{
// 
// 			/* 			std::cout << "using the standard builder and solver " << std::endl; */
// 
// 		}


		/** Destructor.
		*/
		virtual ~TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit(){}


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
 			TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
			LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			
			for(unsigned int ddd=1; ddd<mdim+1; ddd++)
			{
			    //setting the active fractional step number
			    CurrentProcessInfo[FRACTIONAL_STEP] = ddd;
			    
			    // assemble all elements
			    for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			    {
				    //calculate elemental contribution
				    pScheme->CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

				    //assemble the elemental contribution
				    TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
				    TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);

				    // clean local elemental memory
				    pScheme->CleanMemory(*it);
			    }

			    LHS_Contribution.resize(0,0,false);
			    RHS_Contribution.resize(0,false);

			    // assemble all conditions
			    for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			    {
				    //calculate elemental contribution
				    pScheme->Condition_CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

				    //assemble the elemental contribution
				    TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
				    TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
			    }
			
			}

			//finalizing the assembly
			A.GlobalAssemble();
			b.GlobalAssemble();

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

			for(unsigned int ddd=1; ddd<mdim+1; ddd++)
			{
			    //setting the active fractional step number
			    CurrentProcessInfo[FRACTIONAL_STEP] = ddd;
			    
			    // assemble all elements
			    for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			    {
				    //calculate elemental contribution
				    pScheme->Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				    //assemble the elemental contribution
				    TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);

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
				    TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
			    }
			    
			}

			//finalizing the assembly
			A.GlobalAssemble();
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


			double norm_b;
			if(TSparseSpace::Size(b) != 0)
				norm_b = TSparseSpace::TwoNorm(b);
			else
				norm_b = 0.00;

			if(norm_b != 0.00)
			  {
                                if (this->GetEchoLevel()>1) 
                                    if(mrComm.MyPID() == 0) KRATOS_WATCH("entering in the solver");
				
				BaseType::mpLinearSystemSolver->Solve(A,Dx,b);

			  }
			else
			{
				TSparseSpace::SetToZero(Dx);
			}

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

                        int rank = r_model_part.GetCommunicator().MyPID();


			Build(pScheme,r_model_part,A,b);

			if(BaseType::GetEchoLevel()>0)
			{
				if(rank == 0) std::cout << "Building Time : " << building_time.elapsed() << std::endl;
			}
			
			//does nothing...dirichlet conditions are naturally dealt with in defining the residual
			ApplyDirichletConditions(pScheme,r_model_part,A,Dx,b);

			if (BaseType::GetEchoLevel()== 3)
			{
                            if(rank == 0)
                            {
				std::cout << "before the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
                            }
			}

			boost::timer solve_time;

			SystemSolve(A,Dx,b);

			if(BaseType::GetEchoLevel()>0)
			{
				if(rank == 0) std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
			}
			if (BaseType::GetEchoLevel()== 3)
			{
                            if(rank == 0)
                            {
				std::cout << "after the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
                            }
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
			  LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
			  
			  //vector containing the localization in the system of the different
			  //terms
			  Element::EquationIdVectorType EquationId;
			  
			for(unsigned int ddd=1; ddd<mdim+1; ddd++)
			{
			    //setting the active fractional step number
			    CurrentProcessInfo[FRACTIONAL_STEP] = ddd;
			    
			    // assemble all elements
			    for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			      {
				//calculate elemental Right Hand Side Contribution
				pScheme->Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
				
				//assemble the elemental contribution
			      TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
			      }
			    
			    RHS_Contribution.resize(0,false);
			    
			    // assemble all conditions
			    for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			      {
				//calculate elemental contribution
				pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
				
				//assemble the elemental contribution
				TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
				
			      }
			  }
			  
			//finalizing the assembly
			b.GlobalAssemble();

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
// 			  ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
/*  			  ElementsArrayType& pElements = r_model_part.Elements(ModelPart::Kratos_Local); */

			  Element::DofsVectorType ElementalDofList;

//			  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			  DofsArrayType Doftemp;
			  BaseType::mDofSet = DofsArrayType();

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			
			if(mdim == 2)
			{
			  for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin();
								      it != r_model_part.NodesEnd(); it++)
			  {
			      Doftemp.push_back(it->pGetDof(mrVar_x)  );
			      Doftemp.push_back(it->pGetDof(mrVar_y)  );
			  }
			}
			else
			{
			  for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin();
								      it != r_model_part.NodesEnd(); it++)
			  {
			      Doftemp.push_back(it->pGetDof(mrVar_x)  );
			      Doftemp.push_back(it->pGetDof(mrVar_y)  );
			      Doftemp.push_back(it->pGetDof(mrVar_z)  );
			  }
			}
			  
			Doftemp.Unique();

			BaseType::mDofSet = Doftemp;

			//throws an execption if there are no Degrees of freedom involved in the analysis
			if (BaseType::mDofSet.size()==0)
				KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");


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

			if(BaseType::GetEchoLevel()>1)
			{
                            if(rank == 0)
                            {
				std::cout << rank << " : local size = " << BaseType::mDofSet.size() << std::endl; 
				std::cout << rank << " : free_id = " << free_size << std::endl; 
				std::cout << rank << " : fixed_size = " << fixed_size << std::endl; 
				std::cout << rank << " : free_offset = " << free_offset << std::endl; 
				std::cout << rank << " : fixed offset = " << fixed_offset << std::endl;
                            }
			}

			// Now setting the equation id with .
			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
			  if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
			    {
				if (dof_iterator->IsFixed())
					dof_iterator->SetEquationId(fixed_offset++);
				else
					dof_iterator->SetEquationId(free_offset++);
//  				std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
			    }


			BaseType::mEquationSystemSize = global_size;
			mLocalSystemSize = free_size;

			if(BaseType::GetEchoLevel()>1)
			{
                            if(rank == 0)
                            {
				std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl; 
				std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
				std::cout << rank << " : free_offset = " << free_offset << std::endl; 
				std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;
                            }
			}

			//by Riccardo ... it may be wrong!	
			mFirstMyId = free_offset-mLocalSystemSize;
			mLastMyId = mFirstMyId+mLocalSystemSize;
			
			r_model_part.GetCommunicator().SynchronizeDofs();

		}
 
	  void UpdateGhostDofs(ModelPart& rThisModelPart)
		{
		  int rank;
		  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		  int destination=0;

		  vector<int>& neighbours_indices = rThisModelPart.GetCommunicator().NeighbourIndices();

		  for(unsigned int i_domain = 0 ; i_domain <  neighbours_indices.size() ; i_domain++)
		    if((destination = neighbours_indices[i_domain]) >= 0)
		      {
			unsigned int send_buffer_size = 0;
			unsigned int receive_buffer_size = 0;

			NodesArrayType& r_interface_nodes = rThisModelPart.GetCommunicator().LocalMesh(i_domain).Nodes();
			NodesArrayType& r_ghost_nodes = rThisModelPart.GetCommunicator().GhostMesh().Nodes();

			for(typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
			  send_buffer_size += i_node->GetDofs().size();

			for(typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
			  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
			    {
			    receive_buffer_size += i_node->GetDofs().size();
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
			    }


			MPI_Status status;

			if(position > send_buffer_size)
			  std::cout << rank << " Error in estimating send buffer size...." << std::endl;
		  

			int send_tag = 1;//i_domain;
			int receive_tag = 1;//i_domain;

			MPI_Sendrecv (send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
 				 MPI_COMM_WORLD, &status);

			// Updating nodes
			position = 0;
			for(ModelPart::NodeIterator i_node = rThisModelPart.GetCommunicator().GhostMesh().NodesBegin() ; 
			    i_node != rThisModelPart.GetCommunicator().GhostMesh().NodesEnd() ; i_node++)
			  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
			    for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
			  {
			    i_dof->SetEquationId(receive_buffer[position++]);
			  }
			
			if(position > receive_buffer_size)
			  std::cout << rank << " Error in estimating receive buffer size...." << std::endl;
			
			delete [] send_buffer;
			delete [] receive_buffer;
		      }
		  
		}

		//**************************************************************************
		//**************************************************************************
		void ResizeAndInitializeVectors(
			TSystemMatrixPointerType& pA,
			TSystemVectorPointerType& pDx,
			TSystemVectorPointerType& pb,
			ElementsArrayType& rElements,
			ConditionsArrayType& rConditions,
			ProcessInfo& CurrentProcessInfo
			)
		{
			KRATOS_TRY
                            if (this->GetEchoLevel()>1)
				std::cout << "entering ResizeAndInitializeVectors" << std::endl;

				//resizing the system vectors and matrix
				if ( pA == NULL || TSparseSpace::Size1(*pA) == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
				{
					//creating a work array
					unsigned int number_of_local_dofs = mLastMyId - mFirstMyId;

					int temp_size = number_of_local_dofs;
					if(temp_size <1000) temp_size = 1000;
                                        int* temp = new int[temp_size]; //
                                        int* assembling_temp = new int[temp_size];

				
					//generate map - use the "temp" array here 
					for(unsigned int i=0; i!=number_of_local_dofs; i++)
						temp[i] = mFirstMyId+i;
					Epetra_Map my_map(-1, number_of_local_dofs, temp, 0, mrComm);		
					
					
					//create and fill the graph of the matrix --> the temp array is reused here with a different meaning
					Epetra_FECrsGraph Agraph(Copy, my_map, mguess_row_size);
					//int ierr;
					Element::EquationIdVectorType EquationId;
					
					for(unsigned int ddd=1; ddd<mdim+1; ddd++)
					{
					    //setting the active fractional step number
					    CurrentProcessInfo[FRACTIONAL_STEP] = ddd;
					    
					    
					    // assemble all elements
					    for (typename ElementsArrayType::ptr_iterator it=rElements.ptr_begin(); it!=rElements.ptr_end(); it++)
					    {
						    //TODO! this should go through the scheme!!!!!!!!!1
						    (*it)->EquationIdVector(EquationId,CurrentProcessInfo);
						    
						    //filling the list of active global indices (non fixed)
						    unsigned int num_active_indices = 0;
						    for(unsigned int i=0; i<EquationId.size(); i++)
							    if ( EquationId[i] < BaseType::mEquationSystemSize ) 
							    {
								    assembling_temp[num_active_indices] =  EquationId[i];
								    num_active_indices += 1;
							    }

						    if(num_active_indices != 0)
						    {
							    int ierr = Agraph.InsertGlobalIndices(num_active_indices,assembling_temp,num_active_indices, assembling_temp);
                                                            if(ierr < 0) KRATOS_ERROR(std::logic_error,"Epetra failure found in Agraph.InsertGlobalIndices --> ln 837","");
						    }						
					    }
	    
					    // assemble all conditions
					    for (typename ConditionsArrayType::ptr_iterator it=rConditions.ptr_begin(); it!=rConditions.ptr_end(); it++)
					    {
						    //TODO! this should go through the scheme!!!!!!!!!1
						    (*it)->EquationIdVector(EquationId,CurrentProcessInfo);
						    
						    //filling the list of active global indices (non fixed)
						    unsigned int num_active_indices = 0;
						    for(unsigned int i=0; i<EquationId.size(); i++)
							    if ( EquationId[i] < BaseType::mEquationSystemSize ) 
							    {
								    assembling_temp[num_active_indices] =  EquationId[i];
								    num_active_indices += 1;
							    }

						    if(num_active_indices != 0)
						    {
							    int ierr = Agraph.InsertGlobalIndices(num_active_indices,assembling_temp,num_active_indices, assembling_temp);
                                                            if(ierr < 0) KRATOS_ERROR(std::logic_error,"Epetra failure found in Agraph.InsertGlobalIndices --> ln 859","");
						    }
					    }
					}

					//finalizing graph construction
					int graph_assemble_ierr = Agraph.GlobalAssemble();
					if(graph_assemble_ierr != 0) KRATOS_ERROR(std::logic_error,"Epetra failure found","");

				
					//generate a new matrix pointer according to this graph
					TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );	
					pA.swap(pNewA);


					//generate new vector pointers according to the given map
					if( pb == NULL || TSparseSpace::Size(*pb) != BaseType::mEquationSystemSize)
					{
						TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(my_map) );
						pb.swap(pNewb);	
					}
					if( pDx == NULL || TSparseSpace::Size(*pDx) != BaseType::mEquationSystemSize)
					{
						TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(my_map) );
						pDx.swap(pNewDx);	
					}
					if( BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
					{
						TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(my_map) );
						BaseType::mpReactionsVector.swap(pNewReactionsVector); 
					}

					delete [] temp;
					delete [] assembling_temp;



					
				}
				else
				{
					if(TSparseSpace::Size1(*pA) == 0 || TSparseSpace::Size1(*pA) != BaseType::mEquationSystemSize || TSparseSpace::Size2(*pA) != BaseType::mEquationSystemSize)
					{
						KRATOS_ERROR(std::logic_error,"it should not come here resizing is not allowed this way!!!!!!!! ... ","");
					}
				}

				//


				//if needed resize the vector for the calculation of reactions
				if(BaseType::mCalculateReactionsFlag == true)
				{
				  //unsigned int ReactionsVectorSize = BaseType::mDofSet.size()-BaseType::mEquationSystemSize;

// 					BaseType::mpReactionsVector->resize(ReactionsVectorSize,false);

					KRATOS_ERROR(std::logic_error,"calculation of reactions not yet implemented with Trilinos","");
				}


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

		void BuildLHS_CompleteOnFreeRows(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A)
		{
			KRATOS_ERROR(std::logic_error,"method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver ","");
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
			//this->mReactionsVector = TSystemVectorType();

			if (this->GetEchoLevel()>0)
			{

				KRATOS_WATCH("TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit Clear Function called");
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
		Epetra_MpiComm& mrComm;
		int mguess_row_size;
		int mFirstMyId;
		int mLastMyId;

		unsigned int mdim;
		TVariableType const & mrVar_x;
		TVariableType const & mrVar_y;
		TVariableType const & mrVar_z;

		/*@} */
		/**@name Protected Operators*/
		/*@{ */

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
		  KRATOS_ERROR(std::logic_error, "This method is not implemented for Trilinos", "");
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

	}; /* Class TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_COMPONENTWISE_SLIP  defined */

