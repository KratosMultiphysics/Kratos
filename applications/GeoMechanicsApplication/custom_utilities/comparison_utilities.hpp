// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    JMCarbonell,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_COMPARISON_UTILITIES )
#define  KRATOS_GEO_COMPARISON_UTILITIES


/* System includes */

/* Project includes */
#include "custom_utilities/math_utilities.hpp"

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


class KRATOS_API(GEO_MECHANICS_APPLICATION) ComparisonUtilities
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

    double CalculateStressNorm(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        Matrix StressTensor(3,3); //3D stress tensor
        noalias(StressTensor) = ZeroMatrix(3,3);
        for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
        {
            for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
            {
                StressTensor(i,j) = LocalStressTensor(i,j);
            }
        }

        double StressNorm = ((StressTensor(0,0)*StressTensor(0,0))+(StressTensor(1,1)*StressTensor(1,1))+(StressTensor(2,2)*StressTensor(2,2))+
                             (StressTensor(0,1)*StressTensor(0,1))+(StressTensor(0,2)*StressTensor(0,2))+(StressTensor(1,2)*StressTensor(1,2))+
                             (StressTensor(1,0)*StressTensor(1,0))+(StressTensor(2,0)*StressTensor(2,0))+(StressTensor(2,1)*StressTensor(2,1)));

        StressNorm = sqrt(StressNorm);

        return StressNorm;

        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    double CalculateVonMises(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        Matrix StressTensor(3,3); //3D stress tensor
        noalias(StressTensor) = ZeroMatrix(3,3);
        for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
        {
            for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
            {
            StressTensor(i,j) = LocalStressTensor(i,j);
            }
        }

        //in general coordinates:
        double SigmaEquivalent =  (0.5)*((StressTensor(0,0)-StressTensor(1,1))*((StressTensor(0,0)-StressTensor(1,1)))+
                                         (StressTensor(1,1)-StressTensor(2,2))*((StressTensor(1,1)-StressTensor(2,2)))+
                                         (StressTensor(2,2)-StressTensor(0,0))*((StressTensor(2,2)-StressTensor(0,0)))+
                                       6*(StressTensor(0,1)*StressTensor(1,0)+StressTensor(1,2)*StressTensor(2,1)+StressTensor(2,0)*StressTensor(0,2)));

        if( SigmaEquivalent < 0 )
            SigmaEquivalent = 0;

        SigmaEquivalent = sqrt(SigmaEquivalent);

        return SigmaEquivalent;

        KRATOS_CATCH( "" )
    }



    void CheckZeroDiagonalComponents(Matrix& StressTensor)
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

    bool CheckPrincipalStresses(Matrix& StressTensor)
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

#endif /* KRATOS_GEO_COMPARISON_UTILITIES defined */

