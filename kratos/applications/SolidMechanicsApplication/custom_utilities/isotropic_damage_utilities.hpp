//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ISOTROPIC_DAMAGE_UTILITIES )
#define  KRATOS_ISOTROPIC_DAMAGE_UTILITIES


/* System includes */
#include <cmath>

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


class IsotropicDamageUtilities
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
    IsotropicDamageUtilities() {};

    /** Destructor.
     */
    ~IsotropicDamageUtilities() {};

    /** Operators.
     */

    /**
     * @}
     */
    /**
     * Calculates perturbations in each direction of a given vector.
     * @param v1 the given vector used to obtain the perturbations.
     */
    
    Vector PerturbationVector( Vector v1, double min_tol = 1e-10, double max_tol = 1e-5 )
    {
        Vector PerturbationVector = ZeroVector(v1.size());
        
        //Maximum and minimum vector components
        double max_component = fabs(v1[0]) , min_component = fabs(v1[0]);

        for( unsigned int i=1; i<v1.size(); i++ )
        {
            if( fabs(v1[i]) < min_component )
            {
                min_component = fabs(v1[i]);
            }   
            else if( fabs(v1[i]) > max_component )
            {
                max_component = fabs(v1[i]);
            }
        }

        double aux = min_component*max_tol;

        if( aux < (max_component*min_tol) )
        {
            aux = max_component*min_tol;
        }
        
        //PerturbationVector
        for( unsigned int i=0; i<v1.size(); i++ )
        {
            if( fabs(v1[i]) > 1e-20 ) // different from zero
            {
                PerturbationVector[i] = v1[i]*max_tol;
            }
            else if( v1[i] >= 0.0 )
            {
                PerturbationVector[i] = aux;
            }
            else
            {
                PerturbationVector[i] = -aux;
            }
        }
        
        return PerturbationVector;
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

}; /* Class IsotropicDamageUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_ISOTROPIC_DAMAGE_UTILITIES defined */
