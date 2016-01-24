//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_POROMECHANICS_MATH_UTILITIES )
#define  KRATOS_POROMECHANICS_MATH_UTILITIES

/* System includes */
#include <cmath>
#include "utilities/math_utils.h"

namespace Kratos
{

class PoromechanicsMathUtilities
{
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:
    
    //Constructor
    PoromechanicsMathUtilities() {}
    
    //------------------------------------------------------------------------------------
    
    //Destructor
    ~PoromechanicsMathUtilities() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateVonMises(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        Matrix StressTensor = ZeroMatrix(3); //3D stress tensor
        for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
        {
            for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
            {
                StressTensor(i,j) = LocalStressTensor(i,j);
            }
        }

        double SigmaEquivalent = 0.5 * ( (StressTensor(0,0)-StressTensor(1,1))*(StressTensor(0,0)-StressTensor(1,1)) +
                                         (StressTensor(1,1)-StressTensor(2,2))*(StressTensor(1,1)-StressTensor(2,2)) +
                                         (StressTensor(2,2)-StressTensor(0,0))*(StressTensor(2,2)-StressTensor(0,0)) +
                                         6 * ( StressTensor(0,1)*StressTensor(0,1)+StressTensor(1,2)*StressTensor(1,2)+StressTensor(2,0)*StressTensor(2,0) ) );

        if( SigmaEquivalent < 0 ) SigmaEquivalent = 0.0;

        SigmaEquivalent = sqrt(SigmaEquivalent);

        return SigmaEquivalent;

        KRATOS_CATCH( "" )
    }

}; /* Class PoromechanicsMathUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_POROMECHANICS_MATH_UTILITIES defined */
