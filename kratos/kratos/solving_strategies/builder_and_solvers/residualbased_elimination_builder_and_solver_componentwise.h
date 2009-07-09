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
*   Date:                $Date: 2008-11-10 14:23:33 $
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVERCOMPONENTWISE )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVERCOMPONENTWISE


/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif


/* External includes */
#include "boost/smart_ptr.hpp"


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

	This is a specialization of the standard buliding strategy to the case in which a single variable is to be used in the 
	building.

	the creation of the DofList and the construction of the system matrix is in this case much faster
	as the neighborhood relationships are considered to be known


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
	class TDenseSpace , 
	class TLinearSolver,
	class TVariableType
	>
	class ResidualBasedEliminationBuilderAndSolverComponentwise 
		: public ResidualBasedEliminationBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedEliminationBuilderAndSolverComponentwise );


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
		ResidualBasedEliminationBuilderAndSolverComponentwise(
			typename TLinearSolver::Pointer pNewLinearSystemSolver,TVariableType const& Var)
			: ResidualBasedEliminationBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
			, rVar(Var)
		{

			/* 			std::cout << "using the standard builder and solver " << std::endl; */

		}


		/** Destructor.
		*/
		virtual ~ResidualBasedEliminationBuilderAndSolverComponentwise(){}


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
			TSparseSpace::SetToZero( *(BaseType::mpReactionsVector) );

#ifndef _OPENMP

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
			LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			unsigned int pos = (r_model_part.Nodes().begin())->GetDofPosition(rVar);

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
				//calculate elemental contribution
				(*it)->InitializeNonLinearIteration(CurrentProcessInfo);
				(*it)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

				Geometry< Node<3> >& geom = (*it)->GetGeometry();
				if(EquationId.size() != geom.size()) EquationId.resize(geom.size(),false);

                                for(unsigned int i=0; i<geom.size(); i++)
					EquationId[i] = geom[i].GetDof(rVar,pos).EquationId();

                                //assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);
				AssembleRHS(b,RHS_Contribution,EquationId);
			}

			LHS_Contribution.resize(0,0,false);
			RHS_Contribution.resize(0,false);
			EquationId.resize(0,false);

			// assemble all conditions
			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			{
				//calculate elemental contribution
				(*it)->InitializeNonLinearIteration(CurrentProcessInfo);
				(*it)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

				Geometry< Node<3> >& geom = (*it)->GetGeometry();
				if(EquationId.size() != geom.size()) EquationId.resize(geom.size(),false);

                                for(unsigned int i=0; i<geom.size(); i++)
				{
					EquationId[i] = geom[i].GetDof(rVar,pos).EquationId();
				}

				//assemble the elemental contribution
				AssembleLHS(A,LHS_Contribution,EquationId);
				AssembleRHS(b,RHS_Contribution,EquationId);
			}
#else
                        ////////////////////////////////////////////////////////////////////////
                        //********************************************************************//
                        ////////////////////////////////////////////////////////////////////////
                        //OPENMP
                        ////////////////////////////////////////////////////////////////////////
                        //********************************************************************//
                        ////////////////////////////////////////////////////////////////////////

//create a partition of the element array
			int number_of_threads = omp_get_max_threads();

			//creating an array of lock variables of the size of the system matrix
			std::vector< omp_lock_t > lock_array(A.size1());
			


			int A_size = A.size1();
			for(int i = 0; i<A_size; i++)
			    omp_init_lock(&lock_array[i]);
			

                        vector<unsigned int> element_partition;
			CreatePartition(number_of_threads, pElements.size(), element_partition);
                        if (this->GetEchoLevel()>0)
			{
                            KRATOS_WATCH( number_of_threads );
                            KRATOS_WATCH( element_partition );
                        }

                        
                        double start_prod = omp_get_wtime();

			#pragma omp parallel for firstprivate(number_of_threads)
			for(int k=0; k<number_of_threads; k++)
			{
				//contributions to the system
				LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
				LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

				//vector containing the localization in the system of the different
				//terms
				Element::EquationIdVectorType EquationId;
				ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
				typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
				typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];

                                unsigned int pos = (r_model_part.Nodes().begin())->GetDofPosition(rVar);


				// assemble all elements
				for (typename ElementsArrayType::ptr_iterator it=it_begin; it!=it_end; ++it)
				{

                                    //calculate elemental contribution
                                    (*it)->InitializeNonLinearIteration(CurrentProcessInfo);
                                    (*it)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

                                    Geometry< Node<3> >& geom = (*it)->GetGeometry();
                                    if(EquationId.size() != geom.size()) EquationId.resize(geom.size(),false);

                                    for(unsigned int i=0; i<geom.size(); i++)
                                            EquationId[i] = geom[i].GetDof(rVar,pos).EquationId();

                                    //assemble the elemental contribution
				    Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);

/*                                    #pragma omp critical
                                    {
                                        AssembleLHS(A,LHS_Contribution,EquationId);
                                        AssembleRHS(b,RHS_Contribution,EquationId);
                                    }*/
                                }
                        }

			vector<unsigned int> condition_partition;
			CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

 			#pragma omp parallel for firstprivate(number_of_threads)
			for(int k=0; k<number_of_threads; k++)
			{
				//contributions to the system
				LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
				LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

				Condition::EquationIdVectorType EquationId;

				ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

				typename ConditionsArrayType::ptr_iterator it_begin=ConditionsArray.ptr_begin()+condition_partition[k];
				typename ConditionsArrayType::ptr_iterator it_end=ConditionsArray.ptr_begin()+condition_partition[k+1];

                                unsigned int pos = (r_model_part.Nodes().begin())->GetDofPosition(rVar);

				// A all elements
				for (typename ConditionsArrayType::ptr_iterator it=it_begin; it!=it_end; ++it)
				{

                                    //calculate elemental contribution
                                    (*it)->InitializeNonLinearIteration(CurrentProcessInfo);
                                    (*it)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

                                    Geometry< Node<3> >& geom = (*it)->GetGeometry();
                                    if(EquationId.size() != geom.size()) EquationId.resize(geom.size(),false);

                                    for(unsigned int i=0; i<geom.size(); i++)
                                    {
                                            EquationId[i] = geom[i].GetDof(rVar,pos).EquationId();
                                    }

                                    //assemble the elemental contribution
				    Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);
/*                                    #pragma omp critical
                                    {
                                        AssembleLHS(A,LHS_Contribution,EquationId);
                                        AssembleRHS(b,RHS_Contribution,EquationId);
                                    }*/
                                }
                        }
                        if (this->GetEchoLevel()>0)
			{
                            double stop_prod = omp_get_wtime();
                            std::cout << "parallel building time: " << stop_prod - start_prod << std::endl;
                        }

			for(int i = 0; i<A_size; i++)
			    omp_destroy_lock(&lock_array[i]);

#endif



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


			//fills a list of "active" nodes defined as nodes which have neighbours
			// AND no fixed pressure
			mActiveNodes.clear();
			mActiveNodes.reserve(r_model_part.Nodes().size() );
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
					mActiveNodes.push_back(*(it.base() ));
				}
			}

            //fills the DofList and give a unique progressive tag to each node
			BaseType::mDofSet.clear();
			BaseType::mDofSet.reserve(mActiveNodes.size() );

			for(WeakPointerVector< Node<3> >::iterator iii = mActiveNodes.begin(); iii!=mActiveNodes.end(); iii++)
			{
				BaseType::mDofSet.push_back( iii->pGetDof(rVar) );
			}

			//throws an execption if there are no Degrees of freedom involved in the analysis
			if (BaseType::mDofSet.size()==0)
				KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

			BaseType::mDofSetIsInitialized = true;

			KRATOS_CATCH("")
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

				if(pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
				{
					TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0,0) );
					pA.swap(pNewA);					
				}
				if(pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
				{
					TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0) );
					pDx.swap(pNewDx);
				}
				if(pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
				{
					TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0) );
					pb.swap(pNewb);
				}
				if(BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
				{
					TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0) );
					BaseType::mpReactionsVector.swap(pNewReactionsVector); 
				}

				TSystemMatrixType& A  = *pA;
				TSystemVectorType& Dx = *pDx;
				TSystemVectorType& b  = *pb;

				//resizing the system vectors and matrix
				if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
				{
					A.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,false);
					ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
				}
				else
				{
					if(A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
					{
						KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
						A.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,true);
						ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
					}
				}
				if(Dx.size() != BaseType::mEquationSystemSize)
					Dx.resize(BaseType::mEquationSystemSize,false);
				if(b.size() != BaseType::mEquationSystemSize)
					b.resize(BaseType::mEquationSystemSize,false);

				//


				//if needed resize the vector for the calculation of reactions
				if(BaseType::mCalculateReactionsFlag == true)
				{
					unsigned int ReactionsVectorSize = BaseType::mDofSet.size()-BaseType::mEquationSystemSize;
					if(BaseType::mpReactionsVector->size() != ReactionsVectorSize)
						BaseType::mpReactionsVector->resize(ReactionsVectorSize,false);
				}

				//swapping pointers
// 				pA.swap(pNewA);
// 				pDx.swap(pNewDx);
// 				pb.swap(pNewb);

				KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************
		void Clear() 
		{
			this->mDofSet = DofsArrayType();

			TSparseSpace::Clear( (this->mpReactionsVector) );
// 			*(this->mpReactionsVector) = TSystemVectorType();

			if (this->GetEchoLevel()>0)
			{
				KRATOS_WATCH("ResidualBasedEliminationBuilderAndSolver Clear Function called");
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
		void ConstructMatrixStructure(
					      TSystemMatrixType& A,
			ElementsContainerType& rElements,
			ConditionsArrayType& rConditions,
			ProcessInfo& CurrentProcessInfo)
		{ 
			KRATOS_TRY

				std::vector< std::vector<int> > index_list(BaseType::mEquationSystemSize);

			int total_size = 0;

			unsigned int pos = (mActiveNodes.begin())->GetDofPosition(rVar);
			//constructing the system matrix row by row
			int index_i;
			for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin(); 
				in!=mActiveNodes.end(); in++)
			{
				Node<3>::DofType& current_dof = in->GetDof(rVar,pos);
				if( current_dof.IsFixed() == false)
				{
					index_i = (current_dof).EquationId();
					WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);

					std::vector<int>& indices = index_list[index_i];
					indices.reserve(neighb_nodes.size());

					//filling the first neighbours list
					indices.push_back(index_i);
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
						i != neighb_nodes.end(); i++)
					{

						Node<3>::DofType& neighb_dof = i->GetDof(rVar,pos);
						if(neighb_dof.IsFixed() == false )
						{
							int index_j = (neighb_dof).EquationId();
							indices.push_back(index_j);
						}
					}

					//sorting the indices and elminating the duplicates
					std::sort(indices.begin(),indices.end());
					typename std::vector<int>::iterator new_end = std::unique(indices.begin(),indices.end());

//					indices.erase(new_end,indices.end());

					total_size += indices.size();
				}
			}

			A.reserve(total_size,false);

			//setting to zero the matrix (and the diagonal matrix)
			for(unsigned int i=0; i<BaseType::mEquationSystemSize; i++)
			{
				std::vector<int>& indices = index_list[i];
				for(unsigned int j=0; j<indices.size(); j++)
				{
					A.push_back(i,indices[j] , 0.00);
				}
			}	

			KRATOS_CATCH("")
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
		TVariableType const & rVar;
		WeakPointerVector<Node<3> > mActiveNodes;

		/*@} */
		/**@name Private Operators*/
		/*@{ */
                //******************************************************************************************
		//******************************************************************************************
		inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
		{
			partitions.resize(number_of_threads+1);
			int partition_size = number_of_rows / number_of_threads;
			partitions[0] = 0;
			partitions[number_of_threads] = number_of_rows;
			for(int i = 1; i<number_of_threads; i++)
			   partitions[i] = partitions[i-1] + partition_size ;
		}
		
		

		//**************************************************************************
#ifdef _OPENMP
		void Assemble(
			TSystemMatrixType& A,
			TSystemVectorType& b,
			const LocalSystemMatrixType& LHS_Contribution,
			const LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			std::vector< omp_lock_t >& lock_array
			)
		{
			unsigned int local_size = LHS_Contribution.size1();

			for (unsigned int i_local=0; i_local<local_size; i_local++)
			{
				unsigned int i_global=EquationId[i_local];
								
				if ( i_global < BaseType::mEquationSystemSize )
				{
				        omp_set_lock(&lock_array[i_global]);
					
					b[i_global] += RHS_Contribution(i_local);
					for (unsigned int j_local=0; j_local<local_size; j_local++)
					{
						unsigned int j_global=EquationId[j_local];
						if ( j_global < BaseType::mEquationSystemSize )
						{
							A(i_global,j_global) += LHS_Contribution(i_local,j_local);
						}
					}
					
					omp_unset_lock(&lock_array[i_global]);
					
					
				}
			}
		}
#endif


		/*@} */
		/**@name Private Operations*/
		/*@{ */






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

	}; /* Class ResidualBasedEliminationBuilderAndSolverComponentwise */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVERCOMPONENTWISE  defined */

