/* *********************************************************   
*          
*   Last Modified by:    $Author: anonymous $
*   Date:                $Date: 2009-01-15 14:50:24 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_QUASI_INCOMPRESSIBLE_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_QUASI_INCOMPRESSIBLE_BUILDER_AND_SOLVER


/* System includes */
#include <set>


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/geometry_utilities.h" 


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

	\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

	\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


	\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

	\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

	\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

	\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


	*/
	template
	<
	class TSparseSpace,
	class TDenseSpace , 
	class TLinearSolver,
	int TDim
	>
	class ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver 
		: public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		//typedef boost::shared_ptr< ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver );


		typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver> BaseType;

		typedef typename BaseType::TSchemeType TSchemeType;
		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;
		
		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
		typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

		typedef typename BaseType::NodesArrayType NodesContainerType;
		
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
		ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver(
			typename TLinearSolver::Pointer pNewLinearSystemSolver)
			: BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
		{
		}


		/** Destructor.
		*/
		virtual ~ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */


	
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

			KRATOS_ERROR(std::runtime_error, "For the quasi incompressible builder and solver this fct doesnt exist!", "");
			
			KRATOS_CATCH("")
		}

		void InitializeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			KRATOS_WATCH("EMPTY FUNCTION FOR THIS SOLVER")
			KRATOS_CATCH("")
		}

		void FinalizeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			KRATOS_WATCH("EMPTY FUNCTION FOR THIS SOLVER")
			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		//this is done in a purely nodal way taking advantage of the neighbour relatinoships
		//which are assumed to be calculated separately
		//HERE we store the displacements variables in a list
		void SetUpDofSet(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part
			) 
		{
			KRATOS_TRY

			//count dofs
			mnumber_of_active_nodes = 0;
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
					mnumber_of_active_nodes += 1;
				}		
			}

			//getting the dof position
			//unsigned int dof_position = r_model_part.NodesBegin()->GetDofPosition(PRESSURE);

			//fills the DofList
			BaseType::mDofSet.clear();
			BaseType::mDofSet.reserve( mnumber_of_active_nodes * TDim );
			if(TDim == 2)
			{
				for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
				{
					if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
					{
							BaseType::mDofSet.push_back( it->pGetDof(DISPLACEMENT_X) );
							BaseType::mDofSet.push_back( it->pGetDof(DISPLACEMENT_Y) );
					}		
				}
			}
			else if(TDim == 3)
			{
				for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
				{
					if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
					{
							BaseType::mDofSet.push_back( it->pGetDof(DISPLACEMENT_X) );
							BaseType::mDofSet.push_back( it->pGetDof(DISPLACEMENT_Y) );
							BaseType::mDofSet.push_back( it->pGetDof(DISPLACEMENT_Z) );
					}		
				}
			}
			//before it was like that:
			//this->mEquationSystemSize = rDofSet.size();
			this->mEquationSystemSize = BaseType::mDofSet.size();

			//throws an execption if there are no Degrees of freedom involved in the analysis
			if (BaseType::mDofSet.size()==0)
				KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

			BaseType::mDofSetIsInitialized = true;

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		//this function numbers the DOFS - from 1 onwards... (note, that only DISPLACEMENT DOFs are stored, PRESSURE is not!!!!)
		void SetUpSystem(
			ModelPart& r_model_part
			)
		{ 
			KRATOS_TRY
			
			//assing id to the nodes
			unsigned int index = 0;
			for(typename DofsArrayType::iterator i_dof = BaseType::mDofSet.begin() ; i_dof != BaseType::mDofSet.end() ; ++i_dof)
			{	
				//i_dof->EquationId() = index;
				i_dof->SetEquationId(index) ;
				index++;
			}
			KRATOS_CATCH("");
		}

		//**************************************************************************
		//**************************************************************************
		//
	
		
		void ResizeAndInitializeVectors(
			TSystemMatrixType& A,
			TSystemMatrixType& mD,
			TSystemVectorType& Dx,
			TSystemVectorType& b,
			TSystemMatrixType& mMconsistent,
			TSystemVectorType& mMdiagInv,
			ElementsArrayType& rElements,
			ConditionsArrayType& rConditions,
			ProcessInfo& CurrentProcessInfo
			)
		{
			KRATOS_TRY


			//resizing the system vectors and matrix
			if (A.size1() == 0 || this->GetReshapeMatrixFlag() == true) //if the matrix is not initialized
			{
				A.resize(this->mEquationSystemSize,this->mEquationSystemSize,false);
				//ConstructMatrixStructure(A);
			}
			if(Dx.size() != this->mEquationSystemSize)
				Dx.resize(this->mEquationSystemSize,false);
			if(b.size() != this->mEquationSystemSize)
				b.resize(this->mEquationSystemSize,false);
			
			if(BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
			{
				TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0) );
				BaseType::mpReactionsVector.swap(pNewReactionsVector); 
			}

			//resize auxiliaries
			unsigned int reduced_dim = this->mEquationSystemSize / TDim;
			if(mD.size1() != reduced_dim)
				mD.resize(reduced_dim,this->mEquationSystemSize,false);
			
			if(mMconsistent.size1() != reduced_dim)
				mMconsistent.resize(reduced_dim,reduced_dim,false);
			
			if(mMdiagInv.size() != reduced_dim )
				mMdiagInv.resize(reduced_dim,false);
			KRATOS_CATCH("")

		}
		
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
			TSparseSpace::SetToZero( *(BaseType::mpReactionsVector) );

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
			LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;


			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
				
				//calculate elemental contribution
				pScheme->CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);
				AssembleRHS(b,RHS_Contribution,EquationId);
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
				AssembleLHS(A,LHS_Contribution,EquationId);
				AssembleRHS(b,RHS_Contribution,EquationId);
			}

			KRATOS_CATCH("")

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

	public:
		/**@name Static Member Variables */
		/*@{ */


		/*@} */
		/**@name Member Variables */
		/*@{ */

		TSystemMatrixType mD;
		TSystemMatrixType mMconsistent;
		TSystemVectorType mMdiagInv;
		TSystemVectorType mpreconditioner;
		unsigned int mnumber_of_active_nodes;
		/*@} */
		/**@name Private Operators*/
		/*@{ */


		/*@} */
		/**@name Private Operations*/
		/*@{ */

		//**************************************************************************
		//
			//**************************************************************************
		void AssembleLHS(
			TSystemMatrixType& A,
			LocalSystemMatrixType& LHS_Contribution,
			Element::EquationIdVectorType& EquationId
			)
		{
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
							A(i_global,j_global) += LHS_Contribution(i_local,j_local);
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

			if (BaseType::mCalculateReactionsFlag==false) //if we don't need to calculate reactions
			{
				for (unsigned int i_local=0; i_local<local_size; i_local++)
				{
					unsigned int i_global=EquationId[i_local];
					if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs
					{	// ASSEMBLING THE SYSTEM VECTOR
						b[i_global] += RHS_Contribution[i_local];
					}
				}
			}
			else //when the calculation of reactions is needed
			{
			  TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
				for (unsigned int i_local=0; i_local<local_size; i_local++)
				{
					unsigned int i_global=EquationId[i_local];
					if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs
					{	// ASSEMBLING THE SYSTEM VECTOR
						b[i_global] += RHS_Contribution[i_local];
					}
					else //on "fixed" DOFs
					{	// Assembling the Vector of REACTIONS
						ReactionsVector[i_global-BaseType::mEquationSystemSize] -= RHS_Contribution[i_local];
					}
				}
			}
		}


		void ConstructMatrixStructure(
			TSystemMatrixType& A, ModelPart& r_model_part
			)
		{ 
			KRATOS_TRY

			std::vector<int>  indices;
			indices.reserve(1000);

			//count non zeros
			int total_nnz = 0;
			for (typename NodesArrayType::iterator node_iterator =r_model_part.NodesBegin(); node_iterator !=r_model_part.NodesEnd(); ++node_iterator)
			{
				
				//not adding thelonely nodes:
				if( (node_iterator->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{//we add one because we have to account for the contribution of the node itself
				total_nnz +=1+(node_iterator->GetValue(NEIGHBOUR_NODES)).size();
				}
			}
			
			//reserve space in the matrix
			A.reserve(total_nnz* TDim * TDim,false);
			
			unsigned int row_index;

			//fill the matrix row by row
			unsigned int dof_position = r_model_part.NodesBegin()->GetDofPosition(DISPLACEMENT_X);
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES);
				if( neighb_nodes.size() != 0 )
				{
					//first row in the block
					row_index = it->GetDof(DISPLACEMENT_X,dof_position).EquationId();
					
					//add id of the current node
					//NOTE: here and in the following we ASSUME that the ids of DISPLACEMENT_X _Y and _Z are sequential
					for(unsigned int kk = 0; kk<TDim; kk++)
					{
						indices.push_back(row_index + kk);
					}

					//filling and order the first neighbours list
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
						i != neighb_nodes.end(); i++)
					{
						unsigned int tmp = (i->GetDof(DISPLACEMENT_X,dof_position)).EquationId();
						for(unsigned int kk = 0; kk<TDim; kk++)
						{
							indices.push_back(tmp + kk);
						}
					}
					std::sort(indices.begin(),indices.end());

					//fill in the system matrix A
					for(unsigned int kk = 0; kk<TDim; kk++)
					{
						for(unsigned int j=0; j<indices.size(); j++)
						{
							A.push_back(row_index + kk,indices[j] , 0.00);
						}
					}

					//clean the indices (it is a work array)
					indices.erase(indices.begin(),indices.end());
				}
			}

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void ConstructMatrixStructure_Mconsistent(
					TSystemMatrixType& Mconsistent,  ModelPart& r_model_part)
		{
			KRATOS_TRY
		
			std::vector<int>  indices;
			indices.reserve(1000);

			//KRATOS_WATCH("contruct matrix structure Mconsistent 0")
			
			int total_nnz = 0;
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				//not to do include lonely nodes in the system
				if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
				//we add one because we have to account for the contribution of the node itself
				total_nnz += 1+(it->GetValue(NEIGHBOUR_NODES)).size();
				}
			}
			
			Mconsistent.reserve(total_nnz,false);	
			
			unsigned int row_index;
			//fill the matrix row by row
			unsigned int dof_position = r_model_part.NodesBegin()->GetDofPosition(DISPLACEMENT_X);

			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES);
				if( neighb_nodes.size() != 0 )
				{
					//first row in the block
					row_index = it->GetDof(DISPLACEMENT_X,dof_position).EquationId();
					
					//add id of the current node
					//NOTE: here and in the following we ASSUME that the ids of DISPLACEMENT_X _Y and _Z are sequential
					
					//we store in the array of indices the column numbers of the pressure index of the respective node, which coincides
					//with the index of DISP_X, divided by TDim (pressure is scalar - no need to store 2 more indices, as it was in 
					//the case of vector (displ)
					
					//CHECK THIS!!!
					//indices.push_back(row_index/3.0);
					indices.push_back(row_index/TDim);
					
					//filling and order the first neighbours list
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
						i != neighb_nodes.end(); i++)
					{
						unsigned int tmp = (i->GetDof(DISPLACEMENT_X,dof_position)).EquationId();
						indices.push_back(tmp/TDim);
							
					}
					std::sort(indices.begin(),indices.end());

					//fill M (the consistent mass matrix)-note that the "pressure index" is assumed to concide the	DISPLACEMENT_X	index divided by 3
					for(unsigned int j=0; j<indices.size(); j++)
					{
						Mconsistent.push_back(row_index/TDim, indices[j] , 0.00);
						//KRATOS_WATCH(Mconsistent)
					}
					//clean the indices (it is a work array)
					indices.erase(indices.begin(),indices.end());					
				}
			}

			KRATOS_CATCH("")
		}
		
		
		//**************************************************************************
		//**************************************************************************
		void ConstructMatrixStructure_DivergenceMatrixD(
					TSystemMatrixType& mD,  ModelPart& r_model_part)
		{
			KRATOS_TRY
			
			std::vector<int>  indices;
			indices.reserve(1000);

			//count non zeros
			int total_nnz = 0;
			for (typename NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				//not to add lonely nodes to the system
				if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
				//we add one because we have to account for the contribution of the node itself
				total_nnz += 1 + (it->GetValue(NEIGHBOUR_NODES)).size();
				}
			}
			
			mD.reserve(total_nnz* TDim,false);

			unsigned int row_index;
			//fill the matrix row by row
			unsigned int dof_position = r_model_part.NodesBegin()->GetDofPosition(DISPLACEMENT_X);
			
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES);
				if( neighb_nodes.size() != 0 )
				{
					//first row in the block
					row_index = it->GetDof(DISPLACEMENT_X,dof_position).EquationId();
					
					//add id of the current node
					//NOTE: here and in the following we ASSUME that the ids of DISPLACEMENT_X _Y and _Z are sequential
					for(unsigned int kk = 0; kk<TDim; kk++)
					{
						indices.push_back(row_index + kk);
					}

					//filling and order the first neighbours list
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
						i != neighb_nodes.end(); i++)
					{
						unsigned int tmp = (i->GetDof(DISPLACEMENT_X,dof_position)).EquationId();
						for(unsigned int kk = 0; kk<TDim; kk++)
						{
							indices.push_back(tmp + kk);
						}
					}
					std::sort(indices.begin(),indices.end());


					//fill D (the divergence matrix) - note that the "pressure index" is assumed to concide the	DISPLACEMENT_X	index divided by 3
					for(unsigned int j=0; j<indices.size(); j++)
					{
						mD.push_back(row_index/TDim, indices[j] , 0.00);
					}
		
					//clean the indices (it is a work array)
					indices.erase(indices.begin(),indices.end());
				}
			}
			KRATOS_CATCH("")
		
		}
	
		//**************************************************************************
		//**************************************************************************
		
		void BuildAuxiliaries(
			TSystemMatrixType& mD,
			ModelPart& r_model_part)
		{
			KRATOS_TRY
			KRATOS_WATCH("BUILDING AUXILIARY MATRIX D")
			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
			array_1d<double,TDim+1> N;
			array_1d<unsigned int ,TDim+1> local_indices;
			//array_1d<double,TDim+1> rhs_contribution;
			double Volume;
			double temp;
			

			//getting the dof position
			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);

			double aaa = 1.0/(TDim+1.0);
			//if the element is not having all the nodes IS_STRUCTURE, assemble it, otherwise do nothing
			for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin(); 
				i!=r_model_part.ElementsEnd(); i++)
			{	

				
				Geometry< Node<3> >& geom = i->GetGeometry();
				//counting the n-r of structure nodes  
				unsigned int str_nr=0;
				
				//for (int k = 0;k<TDim+1;k++)
				for (unsigned int k = 0;k<geom.size();k++)
				{
				str_nr+=(unsigned int)(i->GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE));
				}
				///////////////////////////////////////////////////////////////////////////////////////////////
				//if the element is not having all the nodes IS_STRUCTURE, assemble it, otherwise do nothing
				// that means, that the entries corresponding to the structural elements are zero
				///////////////////////////////////////////////////////////////////////////////////////////////
				if (geom.size()!=str_nr)
				{				
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
							
					//finiding local indices
					//for(int ii = 0; ii<TDim+1; ii++)
					for(unsigned int ii = 0; ii<geom.size(); ii++)
					{
						local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
					}
					//building matrix D (transpose of the gradient integrated by parts)
							
					temp = Volume*aaa;


					for(unsigned int row = 0; row<TDim+1; row++)
					{
						unsigned int row_index = local_indices[row] / (TDim); //ATTENTION! here i am doing a dangerous op
						//KRATOS_WATCH(row_index)
						for(unsigned int col = 0; col<TDim+1; col++)
						{											
							for(unsigned int kkk = 0; kkk<TDim; kkk++)
							{
								//check if the below is correct (copied it from Mass matrix)
								unsigned int col_index = local_indices[col]+kkk;
								//unsigned int col_index = col + kkk;
								mD(row_index,col_index) += temp * DN_DX(col,kkk);
							}
						}

					
					}
				}
				
			}

			KRATOS_CATCH (" ")
		}
		
		//**************************************************************************
		//**************************************************************************
		//
		//assembles consistent and lumped mass matrices
		void AssembleMassMatrices(TSystemMatrixType& Mconsistent, TSystemVectorType& mMdiagInv,  ModelPart& r_model_part)
		{
			//first we assemble the diagonal mass matrix
			KRATOS_TRY
			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
			array_1d<double,TDim+1> N;
			array_1d<unsigned int ,TDim+1> local_indices;
			//array_1d<double,TDim+1> rhs_contribution;
			double Volume;
			double temp;
			//getting the dof position
			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);

			double aaa = 1.0/(TDim+1.0);
			
			for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin(); 
				i!=r_model_part.ElementsEnd(); i++)
			{	
			
				Geometry< Node<3> >& geom = i->GetGeometry();
				//counting number of structural nodes
				int str_nr=0;
				//for (int k = 0;k<TDim+1;k++)
				for (unsigned int k = 0;k<geom.size();k++)
				{
				str_nr+=int(i->GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE));
				}
				//we do not do anything for the elements of the structure (all nodes are IS_STR)
				if (geom.size()!=str_nr)	
				{
				
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
				
					//finiding local indices
					//for(int ii = 0; ii<TDim+1; ii++)
					for(unsigned int ii = 0; ii<geom.size(); ii++)
					{
						local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
					}	
				
					temp = Volume*aaa;
					for(unsigned int row = 0; row<TDim+1; row++)
					{						
						unsigned int row_index = local_indices[row] / (TDim);
						mMdiagInv[row_index] += temp;						
					}
				}
				
			}
				
			//inverting the mass matrix
			for(unsigned int i = 0; i<TSparseSpace::Size(mMdiagInv); i++)
			{
				if (mMdiagInv[i]*mMdiagInv[i]>0.000000000000001)
					mMdiagInv[i] = 1.0/mMdiagInv[i];
				else{ //if (mMdiagInv[i]==0.0)
					//KRATOS_WATCH(mMdiagInv[i])					
					mMdiagInv[i] = 0.0;					
					}
			}


			//KRATOS_WATCH(mMdiagInv)
			//AND NOW WE BUILD THE CONSISTENT MASS MATRIX
				
			for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin(); 
				i!=r_model_part.ElementsEnd(); i++)
			{	
				
				Geometry< Node<3> >& geom = i->GetGeometry();
				unsigned int str_nr=0;
				for (unsigned int k = 0;k<i->GetGeometry().size();k++)
				{
				  str_nr+=(unsigned int)(i->GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE));
				}
				
				if (geom.size()!=str_nr)
				{
				
					
				
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
				
					//finiding local indices
					//for(int ii = 0; ii<TDim+1; ii++)
					for(unsigned int ii = 0; ii<geom.size(); ii++)
					{
						local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
					}
								
					temp = Volume*aaa;
					//element mass matrix has a shape:
					//			2 1 1
					//  A/12.0* 		1 2 1	in 2D
					//			1 1 2
					//			
					//			and
					//
					//			2 1 1 1
					//	V/20.0* 1 2 1 1		in 3D
					//		    1 1 2 1
					//		    1 1 1 2
				
					//nothing should be added in case of membrane
					for(unsigned int row = 0; row<TDim+1; row++)
						{
							unsigned int row_index = local_indices[row] / (TDim); //pressure is a scalar=>matrix size is Tdim times smaller than for vector
							for(unsigned int col = 0; col<TDim+1; col++)
							{
									unsigned int col_index = local_indices[col] /(TDim);
									if (row_index==col_index)
										{
										//Mconsistent(row_index,col_index) += temp * 2.0;
										if (TDim==2)
												Mconsistent(row_index,col_index) += 0.25*temp * 2.0;
										if (TDim==3)
												Mconsistent(row_index,col_index) += 0.2*temp * 2.0;
										}
									else
										{
								
										//Mconsistent(row_index,col_index) += temp ;
										if (TDim==2)
											Mconsistent(row_index,col_index) += 0.25*temp ;
										else if (TDim==3)
											Mconsistent(row_index,col_index) += 0.2*temp ;
									
										}
							}

					
						}
				}
			}

		
			KRATOS_CATCH("")
			
		}
		
		
		//output += trans(input)*input
		//
		void calc_prod_precond_vec( TSystemVectorType& vec,
							 TSystemVectorType& precond,
							 TSystemVectorType& result) 
		{
			KRATOS_TRY

			if ( precond.size()!=vec.size() ) 
				KRATOS_ERROR(std::logic_error,"preconditioner size is wrong","")
			if ( precond.size()!=result.size() ) 
				KRATOS_ERROR(std::logic_error,"preconditioner size is wrong","")
			TSparseSpace::SetToZero(result);

			//typedef  unsigned int size_type;
			//typedef  double value_type;
			//KRATOS_WATCH(precond)
			for (unsigned int i=0; i<precond.size();i++)
			{
			result[i]=precond[i]*vec[i];
			}
			KRATOS_CATCH("");
		}

		/*

		void calc_prod_vec_precond_vec( TSystemVectorType& vec,
							 TSystemVectorType& precond,
							 double& result) 
		{
			KRATOS_TRY

			if ( precond.size()!=vec.size() ) 
				KRATOS_ERROR(std::logic_error,"preconditioner size is wrong","")
			if ( precond.size()!=result.size() ) 
				KRATOS_ERROR(std::logic_error,"preconditioner size is wrong","")
			
			//typedef  unsigned int size_type;
			//typedef  double value_type;
			//KRATOS_WATCH(precond)
			for (int i=0; i<precond.size();i++)
			{
			result[i]=precond[i]*vec[i];
			}
			result = TSparseSpace::Dot(
			KRATOS_CATCH("");
		}
		*/
		
		void calc_GMinvD_prod(TSystemMatrixType& mD,
							 TSystemVectorType& Minv,
							 TSystemVectorType& x,
							 
							 TSystemVectorType& WorkArray,
							 TSystemVectorType& destination) 
		{
			KRATOS_TRY
			//typedef  unsigned int size_type;
			//typedef  double value_type;

			TSparseSpace::SetToZero(WorkArray);
			//TSparseSpace::SetToZero(destination);

			//WorkArray = D * x
			TSparseSpace::Mult(mD, x, WorkArray);
			
			//KRATOS_WATCH(WorkArray)

			//WorkArray = Minv * WorkArray
			for(unsigned int i=0; i<WorkArray.size(); i++)
			{WorkArray[i] *= Minv[i];}
			
			//destination = trans(D) * WorkArray
			
			
			
			//TSparseSpace::TransposeMult(D, x, WorkArray);
			TSparseSpace::TransposeMult(mD, WorkArray, destination);
			//KRATOS_WATCH(destination)

			KRATOS_CATCH("");
		}
		
		//*************************************************************************************************
		//*************************************************************************************************
		void ReturnDx(	 TSystemVectorType& Dx, TSystemVectorType& xi) 
		{
		KRATOS_TRY
		//TSparseSpace::SetToZero(Dx);
		
		if ( Dx.size()!=xi.size() ) 
			KRATOS_ERROR(std::logic_error,"Dx and xi sizes mismatch","")
		
		Dx=xi;
		KRATOS_CATCH("");
		}

		
		//*************************************************************************************************
		//*************************************************************************************************
		
		void CalculatePreconditionerDiagonalMatrix(const TSystemMatrixType& D,
							 const TSystemVectorType& Minv,
							 const TSystemMatrixType& A,
							 TSystemVectorType& preconditioner) 
		{
			KRATOS_TRY
			typedef  unsigned int size_type;
			typedef double value_type;
			
		

			TSparseSpace::SetToZero(preconditioner);
			
			
			
			if ( preconditioner.size()!=A.size1() ) 
				KRATOS_ERROR(std::logic_error,"preconditioner size is wrong","")
			
			//get diagonal of matrix A
			for(unsigned int i = 0; i<A.size1(); i++)
			{
				preconditioner[i] = A(i,i);
			}
						
			//TSparseSpace::SetToZero(preconditioner);
			
			//calculate and add diagonal of G*Minv*D
			//using that G*Minv*D(i,i) = D_k
			
			
			for (size_type k = 0; k < D.size1 (); ++ k) 
			{
				size_type begin = D.index1_data () [k];
				size_type end = D.index1_data () [k + 1];


				for (size_type i = begin; i < end; ++ i)
				{
					unsigned int index_i = D.index2_data () [i];
					value_type data_i = D.value_data()[i];
					preconditioner[index_i] += Minv[k]*data_i*data_i;
				}
			}
			
			
			//KRATOS_WATCH(preconditioner)
			//invert the preconditioner matrix
			for(unsigned int i = 0; i<A.size1(); i++)
			{
				
				
				
				if (fabs(preconditioner[i])>1e-20)
					//preconditioner[i] = 1.00/preconditioner[i];
					preconditioner[i] = 1.00/preconditioner[i];
				else 
					preconditioner[i] = 1.0;
				
				if (preconditioner[i]<0.0)
				{
				preconditioner[i]=1.0;
					//KRATOS_ERROR(std::logic_error,"NEGATIVE PRECONDITIONER","")
				}
				
				
							
			}
			KRATOS_CATCH("");

		}

		//*************************************************************************************************
		//*************************************************************************************************
		
		bool ConvergenceCheck (TSystemVectorType& residual, TSystemVectorType& b, const double& tolerance, const int& iter_number, const int& max_iter_number) 
		{
		//const DataType abs_toll = DataType(1e-15);
		//
		//absolute tolerance = 1e-15
		//
		if (iter_number>max_iter_number)
				KRATOS_ERROR(std::logic_error,"MAX NUMBER OF ITERATIONS EXCEEDED, UR CG DIDNT CONVERGE","")

		if (TSparseSpace::TwoNorm(residual)<1e-15)
			return true;
		
		
		else 
			{
			const double& ratio = TSparseSpace::TwoNorm(residual)/TSparseSpace::TwoNorm(b);
			//KRATOS_WATCH(ratio)
			return(  (ratio) < tolerance);
			}
		
		}

		//*************************************************************************************************
		//*************************************************************************************************
		//   Additional functions for CG
		/*
		void UpdateX(TSystemVectorType& xi, const double& alpha, TSystemVectorType& di)
		{
		TSparseSpace::UnaliasedAdd(xi, alpha, di);
		}
		
		void UpdateResidual()
		{
		
		}

		void UpdateSearchDirection()
		{
		
		}
		*/

		//*************************************************************************************************
		//*************************************************************************************************

		void ModifyForDirichlet (TSystemMatrixType& A, TSystemVectorType& b) 
		{
			KRATOS_TRY

			double large_number = 1e20;
			for(typename DofsArrayType::iterator i_dof = BaseType::mDofSet.begin() ; i_dof != BaseType::mDofSet.end() ; ++i_dof)
			{	
				if(i_dof->IsFixed() == true)
				{
					unsigned int eq_id = i_dof->EquationId();

					A(eq_id,eq_id) += large_number;
					//b[eq_id] = 0.0001;
				}
			}

			KRATOS_CATCH("");
		}
		void CalculateNodalPressureForce (TSystemMatrixType& mD,TSystemVectorType& mMdiagInv,ModelPart& r_model_part) 
		{
			KRATOS_TRY
			
			int i=0;
			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);
			const int size = TSparseSpace::Size(mMdiagInv);
			
			TSystemVectorType p(size);
			TSystemVectorType f_p(3*size);
			i=0;
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
				if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0)//  && in->FastGetSolutionStepValue(IS_FLUID)==1.0)
				{
				i=in->GetDof(DISPLACEMENT_X,dof_position).EquationId()/TDim;
				p[i]=in->FastGetSolutionStepValue(PRESSURE);
				}
				
			}
			
			TSparseSpace::TransposeMult(mD, p, f_p);
			

			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
				if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0)// && in->FastGetSolutionStepValue(IS_FLUID)==1.0)
				{
					
			
					in->FastGetSolutionStepValue(FORCE_X)=f_p[in->GetDof(DISPLACEMENT_X,dof_position).EquationId()];
					in->FastGetSolutionStepValue(FORCE_Y)=f_p[in->GetDof(DISPLACEMENT_Y,dof_position).EquationId()];
					in->FastGetSolutionStepValue(FORCE_Z)=f_p[in->GetDof(DISPLACEMENT_Z,dof_position).EquationId()];
				
				}
			}
		
		KRATOS_CATCH("");
		}
		/*
		void UpdatePressuresNew (TSystemMatrixType& mMconsistent, TSystemVectorType& mMdiagInv,ModelPart& r_model_part, double bulk_modulus, double density) 
		{
		KRATOS_TRY
			//getting the dof position
			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);
			
			//!!!! LATER ON - CHANGE THE WAY TO COMPUTE BULK MODULUS INSTEAD OF PASSING IT AS A PARAMETER
			//resetting  the pressures to zero
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{				
			in->FastGetSolutionStepValue(PRESSURE)=0.0;				
			}
			//for pressure vectors
			const int size = TSparseSpace::Size(mMdiagInv);
			
			TSystemVectorType p_n(size);
			//TSystemVectorType p_n1(size);
			TSystemVectorType history(size);
			
			//TSparseSpace::SetToZero(p_n1);
			TSparseSpace::SetToZero(p_n);
			TSparseSpace::SetToZero(history);
			

			//assuming that the bulk modulus is the same for all nodes in the model part
			//p_n is the history, d_a - change_of_nodal_area/current_nodal_area
			int i=0;
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
				if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0 )// && in->FastGetSolutionStepValue(IS_FLUID)==1.0)
				{
				i=in->GetDof(DISPLACEMENT_X,dof_position).EquationId()/TDim;
				p_n[i]=in->FastGetSolutionStepValue(PRESSURE,1);
				
				}

			
			} 
			//KRATOS_WATCH(p_n)
						
			//history (multiplied by the consistent mass matrix) and then by the inverse lumped mass matrix
			TSparseSpace::Mult(mMconsistent, p_n, history);			
			
			int aa=0;
			
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{	
			if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0)// && in->FastGetSolutionStepValue(IS_FLUID)==1.0)
				{
				aa=in->GetDof(DISPLACEMENT_X,dof_position).EquationId()/TDim;
					if (in->FastGetSolutionStepValue(IS_FLUID)==1.0)
						{
						in->FastGetSolutionStepValue(PRESSURE)=(mMdiagInv[aa]*history[aa])+bulk_modulus*density*(in->FastGetSolutionStepValue(NODAL_AREA) - in->FastGetSolutionStepValue(NODAL_AREA,1))/(in->FastGetSolutionStepValue(NODAL_AREA));	
						}
				}
					
			}
		
			
			KRATOS_CATCH("");
		}
		*/
//void UpdatePressuresNew (TSystemMatrixType& mMconsistent, TSystemVectorType& mMdiagInv,ModelPart& r_model_part, double bulk_modulus, double density)
void UpdatePressuresNew (TSystemMatrixType& mMconsistent, TSystemVectorType& mMdiagInv,ModelPart& r_model_part, double bulk_modulus, double density)  
		{
		KRATOS_TRY
			//getting the dof position
			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);
//			const double dt = r_model_part.GetProcessInfo()[DELTA_TIME];
			
			//!!!! LATER ON - CHANGE THE WAY TO COMPUTE BULK MODULUS INSTEAD OF PASSING IT AS A PARAMETER
			//resetting  the pressures to zero
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{				
			in->FastGetSolutionStepValue(PRESSURE)=0.0;				
			}
			//for pressure vectors
			const int size = TSparseSpace::Size(mMdiagInv);
			
			TSystemVectorType p_n(size);
			//TSystemMatrixType aux(size,size);
			//aux=ZeroMatrix(size,size);
			TSystemVectorType temp(size);
			TSystemVectorType history(size);
			
			//TSparseSpace::SetToZero(p_n1);
			TSparseSpace::SetToZero(p_n);
			TSparseSpace::SetToZero(history);


			
			
			

			//assuming that the bulk modulus is the same for all nodes in the model part
			//p_n is the history, d_a - change_of_nodal_area/current_nodal_area
			int i=0;
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
				if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0 )// && in->FastGetSolutionStepValue(IS_FLUID)==1.0)
				{
				i=in->GetDof(DISPLACEMENT_X,dof_position).EquationId()/TDim;
				p_n[i]=in->FastGetSolutionStepValue(PRESSURE,1);
				
				}

			
			} 
			//KRATOS_WATCH(p_n)
						
			//history (multiplied by the consistent mass matrix) and then by the inverse lumped mass matrix
			TSparseSpace::Mult(mMconsistent, p_n, history);			
			
			//KRATOS_WATCH(history)
			int aa=0;
			
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{	
			if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0)// && in->FastGetSolutionStepValue(IS_FLUID)==1.0)
				{
				aa=in->GetDof(DISPLACEMENT_X,dof_position).EquationId()/TDim;
					if (in->FastGetSolutionStepValue(IS_FLUID)==1.0)
						{
						//+temp[aa]/density
						in->FastGetSolutionStepValue(PRESSURE)=(mMdiagInv[aa]*history[aa])+bulk_modulus*density*(in->FastGetSolutionStepValue(NODAL_AREA) - in->FastGetSolutionStepValue(NODAL_AREA,1))/(in->FastGetSolutionStepValue(NODAL_AREA));
				
						//this one is without mass matrix difference stab, just the laplacian
						//in->FastGetSolutionStepValue(PRESSURE)=p_n[aa]+temp[aa]/density+bulk_modulus*density*(in->FastGetSolutionStepValue(NODAL_AREA) - in->FastGetSolutionStepValue(NODAL_AREA,1))/(in->FastGetSolutionStepValue(NODAL_AREA));		
						}
				}
					
			}
		
			
			KRATOS_CATCH("");
		}
		//this function updates pressure after the Dx is obtained at every step of N-R procedure
		void UpdatePressures (	TSystemMatrixType& mD, 	
		                             TSystemMatrixType& mMconsistent, TSystemVectorType& mMdiagInv,ModelPart& r_model_part, double bulk_modulus, double density) 
		{
			KRATOS_TRY
			//getting the dof position
//			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);
			const double dt = r_model_part.GetProcessInfo()[DELTA_TIME];
			
			//!!!! LATER ON - CHANGE THE WAY TO COMPUTE BULK MODULUS INSTEAD OF PASSING IT AS A PARAMETER
			
			//for pressure vectors
			const int size = TSparseSpace::Size(mMdiagInv);
			//for displacement vectors
			const int size_disp = TDim*TSparseSpace::Size(mMdiagInv);
			
			TSystemVectorType p_n(size);
			TSystemVectorType dp(size);
			TSystemVectorType p_n1(size);
			TSystemVectorType history(size);
			//TSystemVectorType temp1(size);
			//TSystemVectorType temp2(size);
			
			TSparseSpace::SetToZero(p_n);
			TSparseSpace::SetToZero(dp);
			TSparseSpace::SetToZero(p_n1);
			TSparseSpace::SetToZero(history);
			//TSparseSpace::SetToZero(temp1);
			//TSparseSpace::SetToZero(temp2);

			TSystemMatrixType aux(size,size);
			TSystemVectorType temp(size);


			TSystemVectorType displ(size_disp);
			/*
			TSystemMatrixType GlobLapl (size,size);
			TSystemMatrixType  LocLapl (TDim+1,TDim+1);

			for (typename ElementsArrayType::iterator im=r_model_part.ElementsBegin(); im!=r_model_part.ElementsEnd(); ++im)
			{
			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
			array_1d<double,TDim+1> N;
			array_1d<unsigned int ,TDim+1> local_indices;
			
			Geometry< Node<3> >& geom = im->GetGeometry();
			//calculating elemental values
			double Volume;
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
		
			array_1d<double,3> ms_vel_gauss = ZeroVector(3);

			const array_1d<double,3>& fv0 = geom[0].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& fv1 = geom[1].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& fv2 = geom[2].FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3> fv3 = ZeroVector(3);
			if (TDim==3)
				fv3 = geom[3].FastGetSolutionStepValue(VELOCITY);
				

			double nu = geom[0].FastGetSolutionStepValue(VISCOSITY)+
							geom[1].FastGetSolutionStepValue(VISCOSITY) +
							geom[2].FastGetSolutionStepValue(VISCOSITY);
	
			double density = geom[0].FastGetSolutionStepValue(DENSITY)+
							geom[1].FastGetSolutionStepValue(DENSITY) +
							geom[2].FastGetSolutionStepValue(DENSITY);		

			ms_vel_gauss=fv0+fv1+fv2;
			if (TDim==2)
				{
				nu*=0.33333333333;
				density*=0.33333333333;
				ms_vel_gauss*=0.33333333333;
				}



			if (TDim==3)
				{
				ms_vel_gauss+=fv3;
				nu+=geom[3].FastGetSolutionStepValue(VISCOSITY);
				density+=geom[3].FastGetSolutionStepValue(DENSITY);	
				ms_vel_gauss*=0.25;
				nu*=0.25;
				density*=0.25;
				}
		      

			//finiding local indices
			//for(int ii = 0; ii<TDim+1; ii++)
			for(unsigned int ii = 0; ii<geom.size(); ii++)
			{
				local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
			}
				

			//the structural elements should not contribute to the Laplacian
			int str_nr=0;
			for (unsigned int k = 0;k<geom.size();k++)
			{
			  str_nr+=(unsigned int)(geom[k].FastGetSolutionStepValue(IS_STRUCTURE));
			}
			int switch_var=0;	
			//set to zero the entries of the str. elements
			if (str_nr==TDim+1)
				switch_var=0;
			else
				switch_var =1;

			//ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
			//but with one integration N=0.333333333
			double norm_u;
			double h;
			if (TDim==2)
			{
			ms_vel_gauss[0] =  0.33333333333333*(fv0[0]+fv1[0]+fv2[0]);
			ms_vel_gauss[1] =  0.33333333333333*(fv0[1]+fv1[1]+fv2[1]);
			ms_vel_gauss[2] =  0.0;

			//calculating parameter tau (saved internally to each element)
			h = sqrt(2.00*Volume);
			norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
			norm_u = sqrt(norm_u);
			}
			if (TDim==3)
			{
			ms_vel_gauss[0] =  0.25*(fv0[0]+fv1[0]+fv2[0]+fv3[0]);
			ms_vel_gauss[1] =  0.25*(fv0[1]+fv1[1]+fv2[1]+fv3[1]);
			ms_vel_gauss[2] =  0.25*(fv0[2]+fv1[2]+fv2[2]+fv3[2]);

			//calculating parameter tau (saved internally to each element)
			h = sqrt(2.00*Volume);
			norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1] + ms_vel_gauss[2]*ms_vel_gauss[2];
			norm_u = sqrt(norm_u);
			}
			//- 4.0/(bulk_modulus*h*h)
			//double tau = 1.00 / ( 4.00*nu/(h*h) - bulk_modulus*dt/h+2.00*norm_u/h);	
			//double tau=(bulk_modulus)*dt*h/(norm_u+nu/h);
			//double tau=(bulk_modulus)*dt*dt;//h/(norm_u+*nu/h);						
			//my last proposal
			double tau = (bulk_modulus)*dt*nu/(norm_u*norm_u+(nu/dt));
			//Ric's proposal - doesnt work - checked with 2d-splash 
			//double tau = (bulk_modulus)*dt*1.0/((1.0/dt)+(nu/h*h));

			//SWITCHED OFF THE STABILIZATION!
			switch_var=0;			


			noalias(LocLapl)=switch_var*prod(DN_DX,trans(DN_DX));		
		
			
				
			
			
			

			for(unsigned int row = 0; row<TDim+1; row++)
				{
					unsigned int row_index = local_indices[row] / (TDim); 
					for(unsigned int col = 0; col<TDim+1; col++)
					{
					unsigned int col_index = local_indices[col] /(TDim);

						GlobLapl(row_index, col_index)+=tau*Volume*LocLapl(row,col);							
					}				
				}
			//end of the loop over elements			
			}
			
			for (int i=0;i<mMdiagInv.size(); i++)
			{
			aux(i,i)=GlobLapl(i,i)*mMdiagInv(i);
			}
			*/
			//assuming that the bulk modulus is the same for all nodes in the model part
			//
			//additionally here we update densities, simply by implyimg: ro_0xV_0=ro_1xV_1
			int i=0;
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
			if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
				//in pn we save old pressures
				if (i<size)
					p_n[i]=in->FastGetSolutionStepValue(PRESSURE,1);
				i++;
				//here we update densities
				//if (in->FastGetSolutionStepValue(NODAL_AREA)!=0.0)
				//	in->FastGetSolutionStepValue(DENSITY)=in->FastGetSolutionStepValue(DENSITY,1)*in->FastGetSolutionStepValue(NODAL_AREA,1)/in->FastGetSolutionStepValue(NODAL_AREA);
				
				}
			}
			//temp = prod(aux, p_n);
			//history (multiplied by the consistent mass matrix)
			TSparseSpace::Mult(mMconsistent, p_n, history);
			
			//now we compute the pressure increment
			//first we save in the p_n1 the current deltap = KDd //Dx is denoted by d
			//
			//we store displacements in one big vector 
			
			for(typename DofsArrayType::iterator i_dof = BaseType::mDofSet.begin() ; i_dof != BaseType::mDofSet.end() ; ++i_dof)

			{

					displ[i_dof->EquationId()]=i_dof->GetSolutionStepValue()-i_dof->GetSolutionStepValue(1);							
			}

			
			TSparseSpace::Mult(mD, displ, dp);
			//KRATOS_WATCH(bulk_modulus)
			
			dp*=bulk_modulus; 
			dp*=density;


			//now we add the history (multiplied by the consistent mass matrix)
			//adding: mMconsistent*p_n + KDdipsl
			
			
			//p_n1=(temp+dp);

			//and now we multiply the result with the inverse of the lumped mass matrix
			//we reutilize the auxilliary matrix temp
			
			for (int ii=0; ii<size;ii++)
			{
			//temp1[ii]=mMdiagInv[ii]*p_n1[ii];
			p_n1[ii]=mMdiagInv[ii]*(history[ii]+dp[ii]);
			}
			
			
			//this is just to check
			//for (int ii=0; ii<size;ii++)
			//{
			//temp2[ii]=(mMdiagInv[ii]*dp[ii])+p_n[ii];
			//}
			
			//resetting  the pressures to zero
			//
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
				if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
				in->FastGetSolutionStepValue(PRESSURE)=0.0;
				}
			}
			
			int aa=0;
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{	
			if( (in->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
				//not to add the "lonely" nodes , that are not part of the model (e.g. structure walls)
					if (aa<size)	
						in->FastGetSolutionStepValue(PRESSURE)=p_n1[aa];//+temp[aa]/density;
					   //in->FastGetSolutionStepValue(PRESSURE)=temp2[aa];
				aa++;
				}
			}

			//KRATOS_WATCH("PRESSURE UPDATE FUNCTION INSIDE BULDER AND SOLVER")
			/*
			for (typename NodesArrayType::iterator in=r_model_part.NodesBegin(); in!=r_model_part.NodesEnd(); ++in)
			{
			KRATOS_WATCH(in->FastGetSolutionStepValue(PRESSURE));
			}
			*/
			KRATOS_CATCH("");
		}

		//**************************************************************************
		//**************************************************************************
	/*	
	void SystemSolve(
			const TSystemMatrixType& A,
			const TSystemMatrixType& D,
			const TSystemVectorType& mMass_inverse,
			const TSystemVectorType& mpreconditioner,
			TSystemVectorType& x,
			const TSystemVectorType& b
			)
		{
			KRATOS_TRY

			const int size = TSparseSpaceType::Size(rX);
	
			unsigned int IterationsNumber = 0;
    
			TSystemVectorType r(size);
			TSystemVectorType q(size);

			
	
	PreconditionedMult(rA,rX,r);
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

	BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);
    
	VectorType p(r);
	VectorType q(size);
         
	double roh0 = TSparseSpaceType::Dot(r, r);
	double roh1 = roh0;
	double beta = 0;

	if(fabs(roh0) < 1.0e-30) //modification by Riccardo
//	if(roh0 == 0.00)
	  return false;
	    
	do
	  {
	    PreconditionedMult(rA,p,q);

	    double pq = TSparseSpaceType::Dot(p,q);

	    //if(pq == 0.00)
	    if(fabs(pq) <= 1.0e-30)
	      break;

	    double alpha = roh0 / pq;
        
	    TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
	    TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);
         
	    roh1 = TSparseSpaceType::Dot(r,r);

	    beta = (roh1 / roh0);
	    TSparseSpaceType::ScaleAndAdd(1.00, r, beta, p);
	      
	    roh0 = roh1;

	    BaseType::mResidualNorm = sqrt(roh1);
	    BaseType::mIterationsNumber++;
	  } while(BaseType::IterationNeeded() && (fabs(roh0) > 1.0e-30)
	
	  			
			KRATOS_CATCH("");
		}
	*/

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

	}; /* Class ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_QUASI_INCOMPRESSIBLE_BUILDER_AND_SOLVER  defined */


