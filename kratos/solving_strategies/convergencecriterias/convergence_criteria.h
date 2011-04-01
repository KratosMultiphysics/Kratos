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
 *   Date:                $Date: 2008-02-27 13:54:45 $
 *   Revision:            $Revision: 1.3 $
 *
 * ***********************************************************/


#if !defined(KRATOS_NEW_CONVERGENCE_CRITERIA )
#define  KRATOS_NEW_CONVERGENCE_CRITERIA


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
    class TDenseSpace //= DenseSpace<double>
    >
    class ConvergenceCriteria
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        typedef typename TSparseSpace::DataType TDataType;
        typedef typename TSparseSpace::MatrixType TSystemMatrixType;
        typedef typename TSparseSpace::VectorType TSystemVectorType;

        typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
        typedef typename TDenseSpace::VectorType LocalSystemVectorType;

        typedef Dof<TDataType> TDofType;
        typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> > DofsArrayType;
        /* 		typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType; */

        /** Counted pointer of ConvergenceCriteria */
        //typedef boost::shared_ptr< ConvergenceCriteria< TSparseSpace, TDenseSpace > > Pointer;
        KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriteria);
        /*@} */
        /**@name Life Cycle
         */
        /*@{ */

        /** Constructor.
         */
        ConvergenceCriteria()
        {
            mActualizeRHSIsNeeded = false;
            mConvergenceCriteriaIsInitialized = false;
        }

        /** Destructor.
         */
        virtual ~ConvergenceCriteria()
        {
        }


        /*@} */
        /**@name Operators
         */

        /*@{ */

        void SetActualizeRHSFlag(bool flag)
        {
            mActualizeRHSIsNeeded = flag;
        }

        bool GetActualizeRHSflag()
        {
            return mActualizeRHSIsNeeded;
        }

        /*Criterias that need to be called before getting the solution */
        virtual bool PreCriteria(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                const TSystemMatrixType& A,
                const TSystemVectorType& Dx,
                const TSystemVectorType& b
                )
        {
            return true;
        }

        /*Criterias that need to be called after getting the solution */
        virtual bool PostCriteria(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                const TSystemMatrixType& A,
                const TSystemVectorType& Dx,
                const TSystemVectorType& b
                )
        {
            return true;
        }

        virtual void Initialize(
                ModelPart& r_model_part
                )
        {
            mConvergenceCriteriaIsInitialized = true;
        }

        virtual void InitializeSolutionStep(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                const TSystemMatrixType& A,
                const TSystemVectorType& Dx,
                const TSystemVectorType& b
                )
        {
        }

        virtual void FinalizeSolutionStep(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                const TSystemMatrixType& A,
                const TSystemVectorType& Dx,
                const TSystemVectorType& b
                )
        {
        }
        
        /**
         * This function is designed to be called once to perform all the checks needed
         * on the input provided. Checks can be "expensive" as the function is designed
         * to catch user's errors.
         * @param r_model_part
         * @return 0 all ok
         */
        virtual int Check(ModelPart& r_model_part)
        {
            KRATOS_TRY

            return 0;
            KRATOS_CATCH("");
        }



        bool mActualizeRHSIsNeeded;
        bool mConvergenceCriteriaIsInitialized;


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

} /* namespace Kratos.*/

#endif /* KRATOS_NEW_CONVERGENCE_CRITERIA  defined */

