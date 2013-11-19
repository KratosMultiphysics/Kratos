//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPARISON_UTILITIES )
#define  KRATOS_COMPARISON_UTILITIES


/* System includes */
#include <cmath>

/* Project includes */
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

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


class ComparisonUtilities
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
    ComparisonUtilities() {};

    /** Destructor.
     */
    ~ComparisonUtilities() {};

    /** Operators.
     */

    //**************************************************************************
    //**************************************************************************



    //**************************************************************************
    //**************************************************************************

    double CalculateVonMises(const Vector& StressVector)
    {
        KRATOS_TRY

        double tolerance  = 1e-10;
        double zero       = 1e-10;
        double NormStress = 0.00;


        Matrix StressTensor    = MathUtils<double>::StressVectorToTensor(StressVector);

        CheckZeroDiagonalComponents (StressTensor);

        Vector PrincipalStress = ZeroVector(3);

        NormStress =SolidMechanicsMathUtilities<double>::NormTensor(StressTensor);

        Vector MainStresses;

        bool main_tensor = CheckPrincipalStresses( StressTensor );

        if(!main_tensor)
        {

            if(NormStress>1e-6)
            {
                MainStresses = SolidMechanicsMathUtilities<double>::EigenValues(StressTensor,tolerance,zero);
            }
            else
            {
                MainStresses = ZeroVector(3);
            }

        }
        else
        {
            MainStresses = ZeroVector(3);
            for(unsigned int i=0; i<StressTensor.size1(); i++)
                MainStresses[i]=StressTensor(i,i);
        }


        for(unsigned int i=0; i<MainStresses.size(); i++)
            PrincipalStress[i]=MainStresses[i];

        double SigmaEquivalent =  (0.5)*((PrincipalStress[0]-PrincipalStress[1])*(PrincipalStress[0]-PrincipalStress[1]) +
                                         (PrincipalStress[1]-PrincipalStress[2])*(PrincipalStress[1]-PrincipalStress[2]) +
                                         (PrincipalStress[2]-PrincipalStress[0])*(PrincipalStress[2]-PrincipalStress[0]));

        SigmaEquivalent = sqrt(SigmaEquivalent);

        //std::cout<<" SigmaEquivalent "<<SigmaEquivalent<<" StressVector "<<StressVector<<std::endl;

        return SigmaEquivalent;

        KRATOS_CATCH( "" )

    };



    void  CheckZeroDiagonalComponents (Matrix& StressTensor)
    {
        // No null diagonal terms are accepted in the eigenvalue calculation
        for(unsigned int i=0; i<StressTensor.size1(); i++)
        {
            if (fabs(StressTensor(i,i))<1e-10)
            {
                StressTensor(i,i) = 1e-10;
            }
        }
    }

    bool  CheckPrincipalStresses (Matrix& StressTensor)
    {
        // No null diagonal terms are accepted in the eigenvalue calculation
        bool main = true;
        for(unsigned int i=0; i<StressTensor.size1(); i++)
        {
            for(unsigned int j=0; j<StressTensor.size2(); j++)
            {
                if(i!=j)
                {
                    if (fabs(StressTensor(i,j))>1e-10)
                    {
                        main = false;
                    }
                }
            }
        }

        return main;
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

}; /* Class ComparisonUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_COMPARISON_UTILITIES defined */

