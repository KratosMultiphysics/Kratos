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
 *   Date:                $Date: 2008-11-13 12:12:17 $
 *   Revision:            $Revision: 1.4 $
 *
 * ***********************************************************/


#if !defined(KRATOS_VARIABLE_UTILS )
#define  KRATOS_VARIABLE_UTILS


/* System includes */
#include "includes/define.h"
#include "includes/model_part.h"

/* External includes */


/* Project includes */


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
    class VariableUtils
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        /*@} */
        /**@name Life Cycle 
         */
        /*@{ */

        /** Constructor.
         */


        /** Destructor.
         */

        /*@} */
        /**@name Operators 
         */
        /*@{ */


        /*@} */
        /**@name Operations */
        /*@{ */

        //***********************************************************************
        //***********************************************************************

        void SaveVectorVar(const Variable< array_1d<double, 3 > >& OriginVariable,
                const Variable< array_1d<double, 3 > >& SavedVariable,
                ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
                    
            #pragma omp parallel for
            for(int k=0; k< static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                array_1d<double, 3 > & destination = i->GetValue(SavedVariable);
                noalias(destination) = i->FastGetSolutionStepValue(OriginVariable);
            }
            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************

        void SaveScalarVar(const Variable< double >& OriginVariable,
                Variable< double >& SavedVariable,
                ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k < static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                i->GetValue(SavedVariable) = i->FastGetSolutionStepValue(OriginVariable);
            }
            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************

        void CopyVectorVar(const Variable< array_1d<double, 3 > >& OriginVariable,
                Variable< array_1d<double, 3 > >& DestinationVariable,
                ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k < static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                noalias(i->FastGetSolutionStepValue(DestinationVariable)) = i->FastGetSolutionStepValue(OriginVariable);
            }
            KRATOS_CATCH("")
        }


        //***********************************************************************
        //***********************************************************************

        void CopyScalarVar(const Variable< double >& OriginVariable,
                Variable< double >& DestinationVariable,
                ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k< static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                i->FastGetSolutionStepValue(DestinationVariable) = i->FastGetSolutionStepValue(OriginVariable);
            }
            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************

        void SetToZero_VectorVar(const Variable< array_1d<double, 3 > >& Variable, ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k < static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                noalias(i->FastGetSolutionStepValue(Variable)) = ZeroVector(3);
            }
            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************

        void SetToZero_VelocityVectorVar(ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k< static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                noalias(i->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
//                (i)->Fix(VELOCITY_X);
//                (i)->Fix(VELOCITY_Y);
//                (i)->Fix(VELOCITY_Z);
                //i->FastGetSolutionStepValue(VELOCITY).Fix;
            }

            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************
        void SetToZero_VelocityVectorVar1(ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k < static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                noalias(i->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
//                (i)->Fix(VELOCITY_X);
//                (i)->Fix(VELOCITY_Y);
//                (i)->Fix(VELOCITY_Z);
                //i->FastGetSolutionStepValue(VELOCITY).Fix;
            }

            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************
        void SetToZero_ScalarVar(const Variable< double >& Variable, ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            #pragma omp parallel for
            for(int k=0; k < static_cast<int>(rNodes.size()) ;k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin()+k;
                i->FastGetSolutionStepValue(Variable) = 0.0;
            }
            KRATOS_CATCH("")
        }

        //***********************************************************************
        //***********************************************************************
        //returns a list of nodes filtered using the given var and value
        ModelPart::NodesContainerType SelectNodeList(const Variable< double >& Variable,
                const double& value,
                ModelPart::NodesContainerType& rOriginNodes)
        {
            KRATOS_TRY
            ModelPart::NodesContainerType selected_nodes;
            for (ModelPart::NodesContainerType::iterator i = rOriginNodes.begin(); i != rOriginNodes.end(); i++)
            {
                if (i->FastGetSolutionStepValue(Variable) == value)
                {
                    selected_nodes.push_back(*(i.base()));
                }
            }

            return selected_nodes;
            KRATOS_CATCH("")
        }

        /*@} */
        /**@name Acces */
        /*@{ */


        /*@} */
        /**@name Inquiry */
        /*@{ */


        /*@} */
        /**@name Friends */
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
        /**@name Private  Acces */
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

#endif /* KRATOS_VARIABLE_UTILS  defined */

