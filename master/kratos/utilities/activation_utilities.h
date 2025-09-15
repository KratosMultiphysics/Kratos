//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//

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
            
            const Geometry< Node >& geom = it->GetGeometry();
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
            
            const Geometry< Node >& geom = it->GetGeometry();
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

