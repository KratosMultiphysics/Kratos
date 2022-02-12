//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPARISON_UTILITIES )
#define  KRATOS_COMPARISON_UTILITIES


/* System includes */

/* Project includes */
#include "solid_mechanics_math_utilities.h"


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

        double StressNorm =  ((StressTensor(0,0)*StressTensor(0,0))+(StressTensor(1,1)*StressTensor(1,1))+(StressTensor(2,2)*StressTensor(2,2))+
			      (StressTensor(0,1)*StressTensor(0,1))+(StressTensor(0,2)*StressTensor(0,2))+(StressTensor(1,2)*StressTensor(1,2))+
			      (StressTensor(1,0)*StressTensor(1,0))+(StressTensor(2,0)*StressTensor(2,0))+(StressTensor(2,1)*StressTensor(2,1)));

	StressNorm = sqrt(StressNorm);

        return StressNorm;

        KRATOS_CATCH( "" )

    };

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

	//in principal stresses:

        // double tolerance  = 1e-10;
        // double zero       = 1e-10;
        // double NormStress = 0.00;

        // CheckZeroDiagonalComponents (StressTensor);

        // Vector PrincipalStress(3);
	// noalias(PrincipalStress) = ZeroVector(3);

        // NormStress =SolidMechanicsMathUtilities<double>::NormTensor(StressTensor);

        // Vector MainStresses(3);
        // noalias(MainStresses) = ZeroVector(3);	

        // bool main_tensor = CheckPrincipalStresses( StressTensor );

        // if(!main_tensor)
        // {

        //     if(NormStress>1e-6)
        //     {
        //         MainStresses = SolidMechanicsMathUtilities<double>::EigenValues(StressTensor,tolerance,zero);
        //     }
        //     else
        //     {
        //         noalias(MainStresses) = ZeroVector(3);
        //     }

        // }
        // else
        // {
        //     noalias(MainStresses) = ZeroVector(3);
        //     for(unsigned int i=0; i<StressTensor.size1(); i++)
        //         MainStresses[i]=StressTensor(i,i);
        // }


        // for(unsigned int i=0; i<MainStresses.size(); i++)
        //     PrincipalStress[i]=MainStresses[i];
	

        // double SigmaEquivalent =  (0.5)*((PrincipalStress[0]-PrincipalStress[1])*(PrincipalStress[0]-PrincipalStress[1]) +
        //                                  (PrincipalStress[1]-PrincipalStress[2])*(PrincipalStress[1]-PrincipalStress[2]) +
        //                                  (PrincipalStress[2]-PrincipalStress[0])*(PrincipalStress[2]-PrincipalStress[0]));


        // SigmaEquivalent = sqrt(SigmaEquivalent);

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


