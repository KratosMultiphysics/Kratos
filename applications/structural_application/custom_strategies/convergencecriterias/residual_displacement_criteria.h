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
*   Last Modified by:    $Author:    Nelson$
*   Date:                $Date:      2009-05-22$
*   Revision:            $Revision:  1.0 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_RESIDUAL_DISPLACEMENT_CRITERIA)
#define  KRATOS_NEW_RESIDUAL_DISPLACEMENT_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"

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
	class Residual_Displacement_Criteria : 
						
						public  DisplacementCriteria <TSparseSpace, TDenseSpace>,
						public  ResidualCriteria <TSparseSpace, TDenseSpace>
						
					      
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< Residual_Displacement_Criteria< TSparseSpace, TDenseSpace > > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION(Residual_Displacement_Criteria);

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
		Residual_Displacement_Criteria( TDataType NewRatioTolerance,
			                        TDataType AlwaysConvergedNorm)
			:  
			   DisplacementCriteria<TSparseSpace, TDenseSpace > (NewRatioTolerance,AlwaysConvergedNorm),
			   ResidualCriteria < TSparseSpace, TDenseSpace >  (NewRatioTolerance,AlwaysConvergedNorm)


		{
		       KRATOS_TRY
		       KRATOS_CATCH("")			
		}

		/** Destructor.
		*/
		virtual ~Residual_Displacement_Criteria(){}


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

			mdisplacement_criteria =  this-> DisplacementCriteria< TSparseSpace, TDenseSpace >::PostCriteria(r_model_part,rDofSet,A,Dx,b);
  			mresidual_criteria     =  this-> ResidualCriteria< TSparseSpace, TDenseSpace > ::PostCriteria(r_model_part,rDofSet,A,Dx,b);
                         if(mdisplacement_criteria==true && mresidual_criteria==true)
  			  { 
  			    return true;
  			  } 
                          else
  			  {
  			    return false;
  			  } 

		}

		void Initialize(ModelPart& r_model_part) 
		{
	//		this->DisplacementCriteria::Initialize(r_model_part);
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
		
		bool mdisplacement_criteria;
		bool mresidual_criteria;

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

#endif /* KRATOS_NEW_DISPLACEMENT_CRITERIA  defined */
