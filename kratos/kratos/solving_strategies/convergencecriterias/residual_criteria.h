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
*   Date:                $Date: 2007-03-06 10:30:34 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_RESIDUAL_CRITERIA )
#define  KRATOS_NEW_RESIDUAL_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"

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
	class TDenseSpace = DenseSpace<double>
	>
	class ResidualCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ResidualCriteria */
		typedef boost::shared_ptr< ResidualCriteria< TSparseSpace, TDenseSpace > > Pointer;		

		typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofSetType DofSetType;

		//    typedef typename BaseType::DofArrayType DofArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		//    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		//    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/


		ResidualCriteria(Model::Pointer pNewModel,
			TDataType NewRatioTolerance)
			: ConvergenceCriteria< TSparseSpace, TDenseSpace >(pNewModel)
		{
			mRatioTolerance = NewRatioTolerance;

			//if the norm is smaller than this, convergence is considered achieved
			mAlwaysConvergedNorm = 1e-7;

			//this criteria needs the recalculation of the residual!!
			SetActualizeRHSFlag(true);

			mInitialResidualIsSet=false;
		}

		/** Destructor.
		*/
		virtual ~ResidualCriteria(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

		/*Criterias that need to be called after getting the solution */
		bool PostCriteria(
			const String& ElementGroupName,
			DofSetType& rDofSet,
			const TSystemMatrixType& A,
			const TSystemVectorType& Dx,
			const TSystemVectorType& b,
			const ProcessInfo& CurrentProcessInfo
			)
		{

			if (Dx.size() != 0) //if we are solving for something
			{
				int It = CurrentProcessInfo.GetNonLinearIterationNumber();

				if (mInitialResidualIsSet == false)
				{
					mInitialResidualNorm = sqrt(std::inner_product(b.begin(),b.end(),b.begin(),TDataType()));
					mCurrentResidualNorm = mInitialResidualNorm;
					mInitialResidualIsSet = true;
				}
				else 
				{
					//std::cout << "B = " << b << std::endl;
					mCurrentResidualNorm = sqrt(std::inner_product(b.begin(),b.end(),b.begin(),TDataType()));
				}


				TDataType ratio;

				if(mInitialResidualNorm == 0.00) ratio = 0.00;

				else ratio = mCurrentResidualNorm/mInitialResidualNorm;

				double b_size = b.size();
				std::cout << "RESIDUAL CRITERIA :: ratio = " << ratio << "norm value = " << (mCurrentResidualNorm/sqrt(b_size)) << std::endl;
				if ( 
					ratio <= mRatioTolerance 
					|| 
					(mCurrentResidualNorm/sqrt(b_size))<mAlwaysConvergedNorm
					)  
				{
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


		virtual void InitializeSolutionStep(
			const String& ElementGroupName,
			DofSetType& rDofSet,
			const TSystemMatrixType& A,
			const TSystemVectorType& Dx,
			const TSystemVectorType& b,
			const ProcessInfo& CurrentProcessInfo
			)
		{mInitialResidualIsSet = false;}



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

		bool mInitialResidualIsSet;

		TDataType mRatioTolerance;

		TDataType mInitialResidualNorm;

		TDataType mCurrentResidualNorm;

		TDataType mGuaranteedNormValue;

		TDataType mAlwaysConvergedNorm;

		/*@} */
		/**@name Private Operators*/
		/*@{ */


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

#endif /* KRATOS_NEW_RESIDUAL_CRITERIA  defined */

