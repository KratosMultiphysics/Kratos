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
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_ACTIVATION_UTILITIES )
#define  KRATOS_ACTIVATION_UTILITIES


/* System includes */


/* External includes */


/* Project includes */
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"

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

/// Tool to evaluate the normals on nodes based on the normals of a set of surface conditions
class ActivationUtilities
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
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

    /// This function is used to activate and deactivate elements and conditions
    ///depending on the value of rVariable
    ///we need to distinguish 2 cases:
    /// CASE 1 - active_if_lower_than_reference=true
    ///       here elements and conditions will be active if any of their nodes
    ///       comply with the condition     rVariable<reference_value
    /// CASE 2 - active_if_lower_than_reference=false
    ///       here elements and conditions will be active if any of their nodes
    ///       comply with the condition     rVariable>reference_value    
    ///
    ///@param rmodel_part model_part containing the elements and conditions to be activated/deactivated
    ///@param rVariable nodal variable to be used in the comparison
    ///@param reference_value reference value used in the comparison
    ///@param active_if_lower_than_reference see description above
    void ActivateElementsAndConditions( ModelPart& rmodel_part,
                                        const Variable< double >& rVariable,
                                        const double reference_value,
                                        bool active_if_lower_than_reference)
    {
        KRATOS_TRY

        //KRATOS_WATCH(rVariable);
        ModelPart::ElementsContainerType::iterator el_begin = rmodel_part.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rmodel_part.ConditionsBegin();
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rmodel_part.Elements().size()); i++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + i;
            
            const Geometry< Node<3> >& geom = it->GetGeometry();
            it->Set(ACTIVE,false);
            
            for(unsigned int k=0; k<geom.size(); k++)
            {
                if( geom[k].FastGetSolutionStepValue(rVariable) < reference_value) 
                {
                    it->Set(ACTIVE,true);
                    break;
                }
            }
        }
            
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rmodel_part.Conditions().size()); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = cond_begin + i;
            
            const Geometry< Node<3> >& geom = it->GetGeometry();
            it->Set(ACTIVE,false);
            for(unsigned int k=0; k<geom.size(); k++)
            {
                if( geom[k].FastGetSolutionStepValue(rVariable) < reference_value) 
                {
                    it->Set(ACTIVE,true);
                    break;
                }
                
            }
        }
        
         if( active_if_lower_than_reference == false) //flip everything
        {
             #pragma omp parallel for
            for(int i=0; i<static_cast<int>(rmodel_part.Elements().size()); i++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + i;
                it->Flip(ACTIVE);
            }
            
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(rmodel_part.Conditions().size()); i++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + i;
                it->Flip(ACTIVE);
            }
        }

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

    //ActivationUtilities(void);

    //ActivationUtilities(ActivationUtilities& rSource);


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_ACTIVATION_UTILITIES  defined */

