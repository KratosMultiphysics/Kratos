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
*   Date:                $Date: 2008-11-11 14:03:41 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_TRILINOS_DISPLACEMENT_CRITERIA )
#define  KRATOS_TRILINOS_DISPLACEMENT_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "Epetra_MpiComm.h"

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
	template<class TSparseSpace,
	class TDenseSpace 
	>
	class TrilinosDisplacementCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< TrilinosDisplacementCriteria< TSparseSpace, TDenseSpace > > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION( TrilinosDisplacementCriteria );

		typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

		typedef TSparseSpace SparseSpaceType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		TrilinosDisplacementCriteria(
			TDataType NewRatioTolerance,
			TDataType AlwaysConvergedNorm,
			Epetra_MpiComm& rComm)
			: ConvergenceCriteria< TSparseSpace, TDenseSpace >(),mrComm(rComm)
		{
			mRatioTolerance = NewRatioTolerance;
			mAlwaysConvergedNorm = AlwaysConvergedNorm;

			//mActualizeRHSIsNeeded = false;
		}

		/** Destructor.
		*/
		virtual ~TrilinosDisplacementCriteria(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

		/*Criterias that need to be called after getting the solution */
		bool PostCriteria(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			const TSystemMatrixType& A,
			const TSystemVectorType& Dx,
			const TSystemVectorType& b
			)
		{
		  if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
			{
			  	TDataType mFinalCorrectionNorm = TSparseSpace::TwoNorm(Dx);

				TDataType ratio = 0.00;

				mReferenceDispNorm = CalculateReferenceNorm(rDofSet);

				ratio = mFinalCorrectionNorm/mReferenceDispNorm;

				double aaa = SparseSpaceType::Size(Dx);

				double AbsoluteNorm = (mFinalCorrectionNorm/sqrt(aaa));

				if(mrComm.MyPID() == 0) //print performed only by the first processor
					std::cout << "DISPLACEMENT CRITERIA :: obtained tol = " << ratio << ";  expected ratio = " << mRatioTolerance << "absolute tol = " << AbsoluteNorm << std::endl;

				if ( ratio <= mRatioTolerance || AbsoluteNorm<mAlwaysConvergedNorm ) // || (mFinalCorrectionNorm/x.size())<=1e-7) 
				{
					if(mrComm.MyPID() == 0) //print performed only by the first processor
						std::cout <<"convergence is achieved" <<std::endl;

					return true;
				}
				else
				{
					return false;
				}
			}
			else //in this case all the displacements are imposed!
			{
				return true;
			}
		}

		void Initialize(
			ModelPart& r_model_part
			) 
		{
			BaseType::mConvergenceCriteriaIsInitialized = true;
		}

		void InitializeSolutionStep(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			const TSystemMatrixType& A,
			const TSystemVectorType& Dx,
			const TSystemVectorType& b
			)
		{
		}

		void FinalizeSolutionStep(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			const TSystemMatrixType& A,
			const TSystemVectorType& Dx,
			const TSystemVectorType& b
			){}



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

	private:
		/**@name Static Member Variables */
		/*@{ */


		/*@} */
		/**@name Member Variables */
		/*@{ */
		TDataType mRatioTolerance;
		TDataType mAlwaysConvergedNorm;
		
		Epetra_MpiComm& mrComm;


		TDataType mReferenceDispNorm;
		/*@} */
		/**@name Private Operators*/
		/*@{ */

		TDataType CalculateReferenceNorm(DofsArrayType& rDofSet)
		{
			TDataType local_ReferenceDispNorm = TDataType();
			TDataType temp;

			for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
			{
				if(i_dof->IsFree())
				{
					temp = i_dof->GetSolutionStepValue();
					local_ReferenceDispNorm += temp*temp;
				}
			}
			
			//perform the sum between all of the nodes
			TDataType ReferenceDispNorm = TDataType();
			mrComm.SumAll(&local_ReferenceDispNorm, &ReferenceDispNorm, 1);


			ReferenceDispNorm = sqrt(ReferenceDispNorm);
			return ReferenceDispNorm;
		}

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

	}; /* Class ClassName */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_DISPLACEMENT_CRITERIA  defined */

