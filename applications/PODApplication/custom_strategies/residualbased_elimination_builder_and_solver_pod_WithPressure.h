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
*   Date:                $Date: 2008-11-19 16:12:53 $
*   Revision:            $Revision: 1.10 $
*
* ***********************************************************/


#if !defined(RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_POD_WITHPRESSURE )
#define  RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_POD_WITHPRESSURE



//  REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>
 
 
 #include <boost/numeric/ublas/matrix_proxy.hpp>
 #include <boost/numeric/ublas/matrix_expression.hpp>

 #include <iostream>
 #include <fstream>
 #include <vector>
 
 

/* System includes */
#include <set>
#include <omp.h>


/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "pod_application.h"

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
	class TDenseSpace , //= DenseSpace<double>,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
	class ResidualBasedEliminationBuilderAndSolverPOD_WithPressure 
		: public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		//typedef boost::shared_ptr< ResidualBasedEliminationBuilderAndSolverPOD_WithPressure<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedEliminationBuilderAndSolverPOD_WithPressure );


		typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver> BaseType;

		typedef typename BaseType::TSchemeType TSchemeType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		
		/// POD LHS AND POD RHS and POD X
		

		typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
		typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
 
		typedef typename BaseType::NodesArrayType NodesArrayType;
		typedef typename BaseType::ElementsArrayType ElementsArrayType;
		typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

		typedef typename BaseType::ElementsContainerType ElementsContainerType;
		typedef typename ModelPart::NodesContainerType NodesContainerType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		ResidualBasedEliminationBuilderAndSolverPOD_WithPressure(
			typename TLinearSolver::Pointer pNewLinearSystemSolver)
			: BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
		{
		  
		  

			/* 			std::cout << "using the standard builder and solver " << std::endl; */

		}


		/** Destructor.
		*/
		virtual ~ResidualBasedEliminationBuilderAndSolverPOD_WithPressure(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

		//**************************************************************************
		//**************************************************************************
		
		
		
		 /**
  * Invert a matrix via gauss-jordan algorithm (PARTIAL PIVOT)
  *
  * @param m The matrix to invert. Must be square.
  * @param singular If the matrix was found to be singular, then this
  *        is set to true, else set to false.
  * @return If singular is false, then the inverted matrix is returned.
  *         Otherwise it contains random values.
  */
 template<class T>
 //#define T double /// for debug
 boost::numeric::ublas::matrix<T>
 gjinverse(const boost::numeric::ublas::matrix<T> &m, bool &singular)
 {
     using namespace boost::numeric::ublas;

     const int size = m.size1();

     // Cannot invert if non-square matrix or 0x0 matrix.
     // Report it as singular in these cases, and return 
     // a 0x0 matrix.
     if (size != m.size2() || size == 0)
     {
         singular = true;
         matrix<T> A(0,0);
         return A;
     }

     // Handle 1x1 matrix edge case as general purpose 
     // inverter below requires 2x2 to function properly.
     if (size == 1)
     {
         matrix<T> A(1, 1);
         if (m(0,0) == 0.0)
         {
             singular = true;
             return A;
         }
         singular = false;
         A(0,0) = 1/m(0,0);
         return A;
     }

     // Create an augmented matrix A to invert. Assign the
     // matrix to be inverted to the left hand side and an
     // identity matrix to the right hand side.
     matrix<T> A(size, 2*size);
     matrix_range<matrix<T> > Aleft(A, 
                                    range(0, size), 
                                    range(0, size));
     Aleft = m;
     matrix_range<matrix<T> > Aright(A, 
                                     range(0, size), 
                                     range(size, 2*size));
     Aright = identity_matrix<T>(size);

     // Swap rows to eliminate zero diagonal elements.
     for (int k = 0; k < size; k++)
     {
         if ( A(k,k) == 0 ) // XXX: test for "small" instead
         {
             // Find a row(l) to swap with row(k)
             int l = -1;
             for (int i = k+1; i < size; i++) 
             {
                 if ( A(i,k) != 0 )
                 {
                     l = i; 
                     break;
                 }
             }

             // Swap the rows if found
             if ( l < 0 ) 
             {
                 std::cerr << "Error:" <<  __FUNCTION__ << ":"
                           << "Input matrix is singular, because cannot find"
                           << " a row to swap while eliminating zero-diagonal.";
                 singular = true;
                 return Aleft;
             }
             else 
             {
                 matrix_row<matrix<T> > rowk(A, k);
                 matrix_row<matrix<T> > rowl(A, l);
                 rowk.swap(rowl);

 #if defined(DEBUG) || !defined(NDEBUG)
                 std::cerr << __FUNCTION__ << ":"
                           << "Swapped row " << k << " with row " << l 
                           << ":" << A << "\n";
 #endif
             }
         }
     }

     // Doing partial pivot
     for (int k = 0; k < size; k++)
     {
         // normalize the current row
         for (int j = k+1; j < 2*size; j++)
             A(k,j) /= A(k,k);
         A(k,k) = 1;

         // normalize other rows
         for (int i = 0; i < size; i++)
         {
             if ( i != k )  // other rows  // FIX: PROBLEM HERE
             {
                 if ( A(i,k) != 0 )
                 {
                     for (int j = k+1; j < 2*size; j++)
                         A(i,j) -= A(k,j) * A(i,k);
                     A(i,k) = 0;
                 }
             }
         }

 #if defined(DEBUG) || !defined(NDEBUG)
         std::cerr << __FUNCTION__ << ":"
                   << "GJ row " << k << " : " << A << "\n";
 #endif
     }

     singular = false;
     return Aright;
     
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
			//getting the nodes from the model_part
		        NodesContainerType& rNodes = r_model_part.Nodes(); //gives access to the list of nodes

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
			std::cout << "this is single coeff POD For the 2 velocities and 1 Cooficent for the Pressure builder and solver" << std::endl;
			
			unsigned int number_of_pod_modes = (r_model_part.NodesBegin()->GetValue(POD_VELOCITY_X)).size();
			      
			/// the reduction pod local matrix
			mLHS_POD_Contribution.resize(2*number_of_pod_modes,2*number_of_pod_modes,false);
			mRHS_POD_Contribution.resize(2*number_of_pod_modes,false);
			mX_POD_Contribution.resize(2*number_of_pod_modes,false);
			
			noalias(mLHS_POD_Contribution) = ZeroMatrix(2*number_of_pod_modes,2*number_of_pod_modes);
			noalias(mRHS_POD_Contribution) = ZeroVector(2*number_of_pod_modes);
			noalias(mX_POD_Contribution) = ZeroVector(2*number_of_pod_modes);

			int number_of_threads = OpenMPUtils::GetNumThreads();

                        vector<unsigned int> element_partition;
			CreatePartition(number_of_threads, pElements.size(), element_partition);

			
			#pragma omp parallel for firstprivate(number_of_threads ) schedule(static,1)
			for(int k=0; k<number_of_threads; k++)
			{
				//contributions to the system
				LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
				LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
				
				Matrix local_LHS_POD(2*number_of_pod_modes,2*number_of_pod_modes);
				Vector local_RHS_POD(2*number_of_pod_modes);
				Vector local_X_POD(2*number_of_pod_modes);
				noalias(local_LHS_POD) = ZeroMatrix(2*number_of_pod_modes,2*number_of_pod_modes);
				noalias(local_RHS_POD) = ZeroVector(2*number_of_pod_modes);

				
				//vector containing the localization in the system of the different
				//terms
				Element::EquationIdVectorType EquationId;
				ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
				typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
				typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];
				
				// assemble all elements
				for (typename ElementsArrayType::ptr_iterator it=it_begin; it!=it_end; ++it)
				{
					//calculate elemental contribution
					pScheme->CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

					LocalSystemMatrixType V=LocalSystemMatrixType(2*number_of_pod_modes,LHS_Contribution.size2()); 
					
					for(unsigned int index=0; index<(*it)->GetGeometry().size(); index++)
					{
					    const Vector& vel_x = (*it)->GetGeometry()[index].GetValue(POD_VELOCITY_X);
					    const Vector& vel_y = (*it)->GetGeometry()[index].GetValue(POD_VELOCITY_Y);
					    const Vector& press = (*it)->GetGeometry()[index].GetValue(POD_PRESSURE);
					    unsigned int base = index*3;
					    
					  for(unsigned int i = 0; i<number_of_pod_modes; i++)
					  {
					    V(i,base) = vel_x[i];
					    V(i,base+1) = vel_y[i];
					    V(i,base+2) = 0.0;
					  }  
					  
					    for(unsigned int i =0;i<(number_of_pod_modes); i++)
					    {
					      V(i+number_of_pod_modes,base) =0.0; 
					      V(i+number_of_pod_modes,base+1) =0.0;
					      V(i+number_of_pod_modes,base+2) =press[i]; 
					    }
					  
					  
// 			           KRATOS_WATCH(V); 
				      }
 				  
				      ///working on the LHS
				      Matrix tmp = prod(LHS_Contribution,trans(V));
				      noalias(local_LHS_POD) += prod(V,tmp);
			      
				      /// working on the RHS
				      noalias(local_RHS_POD)+=prod(V,RHS_Contribution);
								
				      pScheme->CleanMemory(*it);
			         }
			         
			         #pragma omp critical
				{
				  noalias(mLHS_POD_Contribution) += local_LHS_POD;
				  noalias(mRHS_POD_Contribution) += local_RHS_POD;
				}
			}
			
			

			
			/// INversion LU         inverting ,LHS_POD_contribution
// 			using namespace boost::numeric::ublas;
 	                typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
// 			typedef boost::numeric::ublas::identity_matrix<std::size_t> identitymatrix;
// 			typedef boost::numeric::ublas::matrix<std::size_t> Matrix;

		  UblasSolve(mLHS_POD_Contribution,mRHS_POD_Contribution,mX_POD_Contribution);

// 		  KRATOS_WATCH(mX_POD_Contribution);
		  
			
			
                         //double EndTime = GetTickCount();

		      //std::cout << "total time " << EndTime - StartTime << std::endl;
		      //std::cout << "writing in the system matrix " << ccc << std::endl;
		      //std::cout << "calculating the elemental contrib " << ddd << std::endl;
		      
		      
			LHS_Contribution.resize(0,0,false);
			RHS_Contribution.resize(0,false);

			/// assemble all conditions
// 			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
// 			{
// 				//calculate elemental contribution
// 				pScheme->Condition_CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);
// 
// 				//assemble the elemental contribution
// 				AssembleLHS(A,LHS_Contribution,EquationId);  // A is the Global matarix after assembling!!
// 				AssembleRHS(b,RHS_Contribution,EquationId);
// 			}
			
			
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
			if( TSparseSpace::Size(b) != 0)
				norm_b = TSparseSpace::TwoNorm(b);
			else
				norm_b = 0.00;

			if(norm_b != 0.00)
				BaseType::mpLinearSystemSolver->Solve(A,Dx,b);
			else
				TSparseSpace::SetToZero(Dx);

			//prints informations about the current time
			if (this->GetEchoLevel()>1)
			{
				std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
			}

			KRATOS_CATCH("")

		}
		
//         template<class T>
//         bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
//  	using namespace boost::numeric::ublas;
//  	typedef permutation_matrix<std::size_t> pmatrix;
//  	// create a working copy of the input
//  	matrix<T> A(input);
//  	// create a permutation matrix for the LU-factorization
//  	pmatrix pm(A.size1());
// 
//  	// perform LU-factorization
//  	int res = lu_factorize(A,pm);
//         if( res != 0 ) return false;
// 
//  	// create identity matrix of "inverse"
//  	inverse.assign(ublas::identity_matrix<T>(A.size1()));
// 
//  	// backsubstitute to get the inverse
//  	lu_substitute(A, pm, inverse);
// 
//  	return true;
//        }


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

			Timer::Start("Build");
			Build(pScheme,r_model_part,A,b);
			Timer::Stop("Build");

			//does nothing...dirichlet conditions are naturally dealt with in defining the residual
			ApplyDirichletConditions(pScheme,r_model_part,A,Dx,b);

			if (this->GetEchoLevel()== 3)
			{
				std::cout << "before the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
			}

// 			boost::timer solve_time;
			Timer::Start("Solve");
			
			/// here we set the solution step value to zero, so to be sure that the update will do nothing
			TSparseSpace::SetToZero( Dx);
			
			for(NodesContainerType::iterator it = r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
			{
			   const Vector& vel_pod_x = (it)->GetValue(POD_VELOCITY_X);
			   const Vector& vel_pod_y = (it)->GetValue(POD_VELOCITY_Y);
			   const Vector& press_pod = (it)->GetValue(POD_PRESSURE);
			   array_1d<double,3>& vel = it->FastGetSolutionStepValue(VELOCITY);


			  
			  if((it)->IsFixed(VELOCITY_X) == false)
 			  {
			      for(unsigned int i = 0; i<vel_pod_x.size(); i++)
			         vel[0] += vel_pod_x[i]*mX_POD_Contribution[i];
 			  }
				    
 			  if((it)->IsFixed(VELOCITY_Y) == false)
 			  {
			    for(unsigned int i = 0; i<vel_pod_x.size(); i++)
				    { 
			         vel[1] += vel_pod_y[i]*mX_POD_Contribution[i];
// 			         KRATOS_WATCH(vel_pod_y[i]);
				    }
 			  }
 			  vel[2] = 0.0;
			  
 			  if((it)->IsFixed(PRESSURE) == false)
			  { 			   
			  
			  unsigned int base2=0; 
			  for(unsigned int i = vel_pod_x.size(); i<2*press_pod.size(); i++)
				    { 
			         it->FastGetSolutionStepValue(PRESSURE) += press_pod[base2++]*mX_POD_Contribution[i];
			    
				    }
			  }
// 			  KRATOS_WATCH(vel);
			  
			  
			  
			}
			
// 			KRATOS_WATCH("fgdjkgfhjdfhjkhsdjkfgh");
// 			KRATOS_WATCH(VELOCITY);
// 			KRATOS_WATCH(PRESSURE);
// 			KRATOS_WATCH(POD_VELOCITY_X);

// KRATOS_WATCH(r_model_part.Nodes()[587].FastGetSolutionStepValue(VELOCITY));
// 	

// 			SystemSolve(mLHS_POD_Contribution,mX_POD_Contribution,mRHS_POD_Contribution);

// 			SystemSolve(A,Dx,b);

			Timer::Stop("Solve");

// 			if(this->GetEchoLevel()>0)
// 			{
// 				std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
// 			}
			if (this->GetEchoLevel()== 3)
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
			TSparseSpace::SetToZero( *(BaseType::mpReactionsVector) );

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

KRATOS_WATCH("setting up the dofs");
				//Gets the array of elements from the modeler
				ElementsArrayType& pElements = r_model_part.Elements();

			Element::DofsVectorType ElementalDofList;

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			DofsArrayType Doftemp;
			BaseType::mDofSet = DofsArrayType();
			//mDofSet.clear();

			//double StartTime = GetTickCount();
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
				// gets list of Dof involved on every element
//aaa = GetTickCount();
				pScheme->GetElementalDofList(*it,ElementalDofList,CurrentProcessInfo);
//bbb += GetTickCount() - aaa;
/*KRATOS_WATCH((*it)->Id());
std::cout << "node ids" << std::endl;
for(unsigned int i=0; i<((*it)->GetGeometry()).size(); i++)
	std::cout << ((*it)->GetGeometry())[i].Id() << " ";
std::cout << std::endl;
for(unsigned int i=0; i<ElementalDofList.size(); i++)
	std::cout << (ElementalDofList[i]->Id()) << " ";
std::cout << std::endl;*/

//KRATOS_WATCH(ElementalDofList);

//ccc = GetTickCount();
			for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
				{
					Doftemp.push_back(*i);
					//mDofSet.push_back(*i);
				}
//ddd += GetTickCount() - ccc;
			}
	
//std::cout << "searching " << bbb << std::endl;
//std::cout << "inserting " << ddd << std::endl;

			//taking in account conditions
			ConditionsArrayType& pConditions = r_model_part.Conditions();
			for (typename ConditionsArrayType::ptr_iterator it=pConditions.ptr_begin(); it!=pConditions.ptr_end(); ++it)
			{
				// gets list of Dof involved on every element
				pScheme->GetConditionDofList(*it,ElementalDofList,CurrentProcessInfo);

//ccc = GetTickCount();
				for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
				{
					//mDofSet.push_back(*i);
					Doftemp.push_back(*i);
				}
//ddd += GetTickCount() - ccc;
			}
//std::cout << "searching " << bbb << std::endl;
//std::cout << "inserting " << ddd << std::endl;
/*for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
{
	KRATOS_WATCH(*dof_iterator);
}
std::cout << "DofTemp before Unique" << Doftemp.size() << std::endl;
*/
//ccc = GetTickCount();
			Doftemp.Unique();
//std::cout << "DofTemp after Unique" << Doftemp.size() << std::endl;
			BaseType::mDofSet = Doftemp;

//ddd = GetTickCount() - ccc;
//std::cout << "Unique " << ddd << std::endl;

			//throws an execption if there are no Degrees of freedom involved in the analysis
			if (BaseType::mDofSet.size()==0)
				KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

			BaseType::mDofSetIsInitialized = true;
KRATOS_WATCH("finished setting up the dofs");
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
			int free_id = 0;
			int fix_id = BaseType::mDofSet.size();

			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
				if (dof_iterator->IsFixed())
					dof_iterator->SetEquationId(--fix_id);
				else
					dof_iterator->SetEquationId(free_id++);

			BaseType::mEquationSystemSize = fix_id;

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
					ConstructMatrixStructure(A,rElements,rConditions,CurrentProcessInfo);
				}
				else
				{
					if(A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
					{
						KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
						A.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,true);
						ConstructMatrixStructure(A,rElements,rConditions,CurrentProcessInfo);
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
			//refresh RHS to have the correct reactions
			BuildRHS(pScheme,r_model_part,b);

			int i;
			int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);

			typename DofsArrayType::ptr_iterator it2;
// KRATOS_WATCH(*BaseType::mpReactionsVector);
			//updating variables
			TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
			for (it2=BaseType::mDofSet.ptr_begin();it2 != BaseType::mDofSet.ptr_end(); ++it2)
			{
				if ( (*it2)->IsFixed()  )
				{
					i=(*it2)->EquationId();
					i-=systemsize; 
/*KRATOS_WATCH((*it2)->GetSolutionStepReactionValue());
KRATOS_WATCH(ReactionsVector[i]);*/
					(*it2)->GetSolutionStepReactionValue() = ReactionsVector[i];
				}
			}
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

			if(this->mpReactionsVector != NULL)
				TSparseSpace::Clear( (this->mpReactionsVector) );
// 			this->mReactionsVector = TSystemVectorType();

			if (this->GetEchoLevel()>0)
			{

				KRATOS_WATCH("ResidualBasedEliminationBuilderAndSolverPOD_WithPressure Clear Function called");
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

			std::size_t equation_size = A.size1();
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
			int data_size = 0;
			for(std::size_t i = 0 ; i < indices.size() ; i++)
			{
				data_size += indices[i].size();
			}
			A.reserve(data_size,false);

            //filling with zero the matrix (creating the structure)
            Timer::Start("MatrixStructure");
#ifndef _OPENMP
            for(std::size_t i = 0 ; i < indices.size() ; i++)
            {
                std::vector<std::size_t>& row_indices = indices[i];
                std::sort(row_indices.begin(), row_indices.end());
                
                for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
                {
                    A.push_back(i,*it,0.00);
                }
                row_indices.clear(); 
            }
#else
            int number_of_threads = omp_get_max_threads();
            vector<unsigned int> matrix_partition;
            CreatePartition(number_of_threads, indices.size(), matrix_partition);
            KRATOS_WATCH( matrix_partition );
            for( int k=0; k<number_of_threads; k++ )
            {
                #pragma omp parallel
                if( omp_get_thread_num() == k )
                {
                    for( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; i++ )
                    {
                        std::vector<std::size_t>& row_indices = indices[i];
                        std::sort(row_indices.begin(), row_indices.end());
                        
                        for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
                        {
                            A.push_back(i,*it,0.00);
                        }
                        row_indices.clear(); 
                    }
                }
            }
#endif
            Timer::Stop("MatrixStructure");
        }

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
		
		
		int number_of_pod_modes;
		
		LocalSystemMatrixType mLHS_POD_Contribution ;
   		LocalSystemVectorType mRHS_POD_Contribution ;
		LocalSystemVectorType mX_POD_Contribution ;
		LocalSystemVectorType Solution ;
// 		LocalSystemMatrixType inverse;
		 

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
			unsigned int local_size = LHS_Contribution.size1();
			for (unsigned int i_local=0; i_local<local_size; i_local++)
			{
				unsigned int i_global=EquationId[i_local];
				if ( i_global < BaseType::mEquationSystemSize )
				{
					for (unsigned int j_local=0; j_local<local_size; j_local++)
					{
						int j_global=EquationId[j_local];

						A(i_global,j_global) += LHS_Contribution(i_local,j_local);
					}
				}
			}
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
				//note that computation of reactions is not performed here!
			}
		}
#endif

		bool UblasSolve (const Matrix& input, const Vector& rhs, Vector& x)
		{
		    using namespace boost::numeric::ublas;
		    typedef permutation_matrix<std::size_t> pmatrix;
		    // create a working copy of the input
		    Matrix A(input);
		    // create a permutation matrix for the LU-factorization
		    pmatrix pm(A.size1());

		    // perform LU-factorization
		    int res = lu_factorize(A,pm);
		    if( res != 0 ) return false;
		    
		    noalias(x) = rhs;

		    // backsubstitute to get the inverse
		    lu_substitute(A,pm,x);
		    
		    Vector aaa=rhs;
		    aaa -= prod(input,x);
			    
		    return true;
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

	}; /* Class ResidualBasedEliminationBuilderAndSolverPOD_WithPressure */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_POD_WITHPRESSURE  defined */

