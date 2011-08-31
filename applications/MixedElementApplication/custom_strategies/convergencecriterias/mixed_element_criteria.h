/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
/* *********************************************************   
 *
 *   Last Modified by:    $Author: janosch $
 *   Date:                $Date: 2007-04-13 15:59:32 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/


#if !defined(KRATOS_MIXED_ELEMENT_CRITERIA )
#define  KRATOS_MIXED_ELEMENT_CRITERIA


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
    class TDenseSpace
    >
    class MixedElementConvergenceCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        typedef boost::shared_ptr< MixedElementConvergenceCriteria< TSparseSpace, TDenseSpace > > Pointer;

        typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

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
        MixedElementConvergenceCriteria(
                TDataType NewRatioTolerance,
                TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
        {
            mRatioTolerance = NewRatioTolerance;
            mAlwaysConvergedNorm = AlwaysConvergedNorm;

            //mActualizeRHSIsNeeded = false;
        }

        /** Destructor.
         */
        virtual ~MixedElementConvergenceCriteria()
        {
        }


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
            if (Dx.size() != 0) //if we are solving for something
            {

                double delta_disp_norm = 0.0;
                double disp_norm = 0.0;

                double delta_stress_norm = 0.0;
                double stress_norm = 0.0;

                for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
                {
                    if (i_dof->IsFree())
                    {
                        if (i_dof->GetVariable() == DISPLACEMENT_X || i_dof->GetVariable() == DISPLACEMENT_Y || i_dof->GetVariable() == DISPLACEMENT_Z)
                        {
                            delta_disp_norm += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                            disp_norm += i_dof->GetSolutionStepValue() * i_dof->GetSolutionStepValue();
                        } else
                        {
                            delta_stress_norm += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                            stress_norm += i_dof->GetSolutionStepValue() * i_dof->GetSolutionStepValue();
                        }

                    }
                }

                delta_disp_norm = sqrt(delta_disp_norm);
                disp_norm = sqrt(disp_norm);
                delta_stress_norm = sqrt(delta_stress_norm);
                stress_norm = sqrt(stress_norm);

                double disp_ratio = delta_disp_norm / disp_norm;
                double stress_ratio = delta_stress_norm / stress_norm;
                std::cout << "delta_disp_norm = " << delta_disp_norm << ";  disp_norm = " << disp_norm << std::endl;
               std::cout << "delta_stress_norm = " << delta_stress_norm << ";  stress_norm = " << stress_norm << std::endl;

                std::cout << "DISPLACEMENTS :: disp ratio = " << disp_ratio << ";  Expected ratio = " << mRatioTolerance << "Absolute tol reached = " << delta_disp_norm << std::endl;
                std::cout << "OTHER ::              ratio = " << stress_ratio << ";  Expected ratio = " << mRatioTolerance << "Absolute tol reached = " << delta_stress_norm << std::endl;

                if (
                        (disp_ratio < mRatioTolerance || delta_disp_norm < mAlwaysConvergedNorm)
                        &&
                        (delta_stress_norm < mRatioTolerance || delta_stress_norm < mAlwaysConvergedNorm)
                        )
                {
                    std::cout << "Congratulations the time step solution is converged..." << std::endl;
                    return true;
                } else
                {
                    return false;
                } 

            }

            return true;
        }

        void Initialize(
                ModelPart & r_model_part
                )
        {
        }

        void InitializeSolutionStep(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                const TSystemMatrixType& A,
                const TSystemVectorType& Dx,
                const TSystemVectorType & b
                )
        {
        }

        void FinalizeSolutionStep(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                const TSystemMatrixType& A,
                const TSystemVectorType& Dx,
                const TSystemVectorType & b
                )
        {
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

    private:
        /**@name Static Member Variables */
        /*@{ */


        /*@} */
        /**@name Member Variables */
        /*@{ */
        TDataType mRatioTolerance;
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

} /* namespace Kratos.*/

#endif /* KRATOS_MIXED_ELEMENT_CRITERIA  defined */

