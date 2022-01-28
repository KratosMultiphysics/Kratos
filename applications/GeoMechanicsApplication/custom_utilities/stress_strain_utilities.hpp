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


class KRATOS_API(GEO_MECHANICS_APPLICATION) StressStrainUtilities
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
    StressStrainUtilities() {};

    /** Destructor.
     */
    ~StressStrainUtilities() {};

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
        for(unsigned int i=0; i<LocalStressTensor.size1(); ++i)
        {
            for(unsigned int j=0; j<LocalStressTensor.size2(); ++j)
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

    double CalculateVonMisesStress(const Vector& StressVector)
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

        //in general coordinates:
        double SigmaEquivalent =  (0.5)*((StressTensor(0,0)-StressTensor(1,1))*((StressTensor(0,0)-StressTensor(1,1)))+
                                         (StressTensor(1,1)-StressTensor(2,2))*((StressTensor(1,1)-StressTensor(2,2)))+
                                         (StressTensor(2,2)-StressTensor(0,0))*((StressTensor(2,2)-StressTensor(0,0)))+
                                       6*(StressTensor(0,1)*StressTensor(1,0)+StressTensor(1,2)*StressTensor(2,1)+StressTensor(2,0)*StressTensor(0,2)));

        if ( SigmaEquivalent < 0 )
             SigmaEquivalent = 0;

        SigmaEquivalent = sqrt(SigmaEquivalent);

        return SigmaEquivalent;

        KRATOS_CATCH( "" )
    }

    double CalculateTrace(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix StressTensor = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        double trace = 0.0;
        for (std::size_t i=0; i < StressTensor.size1(); ++i) {
            trace += StressTensor(i,i);
        }

        return trace;

        KRATOS_CATCH( "" )
    }

    double CalculateMeanStress(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix StressTensor = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        double trace = 0.0;
        for (std::size_t i=0; i < StressTensor.size1(); ++i) {
            trace += StressTensor(i,i);
        }

        return (trace / StressTensor.size1());

        KRATOS_CATCH( "" )
    }


    double CalculateVonMisesStrain(const Vector& StrainVector)
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

} /* namespace Kratos.*/

#endif /* KRATOS_GEO_COMPARISON_UTILITIES defined */

