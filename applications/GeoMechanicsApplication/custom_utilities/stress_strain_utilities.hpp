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

#pragma once

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


class KRATOS_API(GEO_MECHANICS_APPLICATION) StressStrainUtilities
{
public:
    /**@name Type Definitions */
    /*@{ */


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Operators.
     */

    static double CalculateStressNorm(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix LocalStressTensor = MathUtils<>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        double StressNorm = 0.;
        for(unsigned int i = 0; i < LocalStressTensor.size1(); ++i) {
            for(unsigned int j = 0; j < LocalStressTensor.size2(); ++j) {
                StressNorm += LocalStressTensor(i,j)*LocalStressTensor(i,j);
            }
        }
        return std::sqrt(StressNorm);

        KRATOS_CATCH( "" )
    }

    static double CalculateVonMisesStress(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        Matrix StressTensor(3,3); //3D stress tensor
        noalias(StressTensor) = ZeroMatrix(3,3);
        for (std::size_t i=0; i < LocalStressTensor.size1(); ++i) {
            for (std::size_t j=0; j < LocalStressTensor.size2(); ++j) {
            StressTensor(i,j) = LocalStressTensor(i,j);
            }
        }

        double SigmaEquivalent = 0.5*((StressTensor(0,0)-StressTensor(1,1))*(StressTensor(0,0)-StressTensor(1,1))+
                                      (StressTensor(1,1)-StressTensor(2,2))*(StressTensor(1,1)-StressTensor(2,2))+
                                      (StressTensor(2,2)-StressTensor(0,0))*(StressTensor(2,2)-StressTensor(0,0))+
                                      6.0*(StressTensor(0,1)*StressTensor(1,0)+
                                           StressTensor(1,2)*StressTensor(2,1)+
                                           StressTensor(2,0)*StressTensor(0,2) ));

        return std::sqrt(std::max(SigmaEquivalent, 0.));

        KRATOS_CATCH( "" )
    }

    static double CalculateTrace(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix StressTensor = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        double trace = 0.0;
        for (std::size_t i = 0; i < StressTensor.size1(); ++i) {
            trace += StressTensor(i,i);
        }

        return trace;

        KRATOS_CATCH( "" )
    }

    static double CalculateMeanStress(const Vector& StressVector)
    {
        KRATOS_TRY

        return CalculateTrace(StressVector) / (StressVector.size() == 3 ? 2.0 : 3.0);

        KRATOS_CATCH( "" )
    }

    static double CalculateVonMisesStrain(const Vector& StrainVector)
    {
        KRATOS_TRY

        return (2.0/3.0) * CalculateVonMisesStress(StrainVector);

        KRATOS_CATCH( "" )
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

}; /* Class StressStrainUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}