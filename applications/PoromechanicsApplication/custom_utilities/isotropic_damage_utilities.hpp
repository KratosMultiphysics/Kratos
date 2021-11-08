//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ISOTROPIC_DAMAGE_UTILITIES_H_INCLUDED)
#define  KRATOS_ISOTROPIC_DAMAGE_UTILITIES_H_INCLUDED


/* System includes */
#include <cmath>

namespace Kratos
{

class IsotropicDamageUtilities
{
public:
    /**
     * @name type definitions
     * @{
     */

    /**
     * @}
     */
    /**
     * Calculates perturbations in each direction of a given vector.
     * @param InputVector the given vector used to obtain the perturbations.
     */

    static inline void ComputePerturbationVector( Vector& rPerturbationVector, const Vector& InputVector )
    {
        const unsigned int VSize = InputVector.size();
        if(rPerturbationVector.size() != VSize)
            rPerturbationVector.resize(VSize,false);

        const double MinTol = 1.0e-10;
        const double MaxTol = 1.0e-5;

        //Maximum and minimum vector components
        double max_component = fabs(InputVector[0]) , min_component = fabs(InputVector[0]);

        for( unsigned int i=1; i<VSize; i++ )
        {
            if( fabs(InputVector[i]) < min_component )
            {
                min_component = fabs(InputVector[i]);
            }
            else if( fabs(InputVector[i]) > max_component )
            {
                max_component = fabs(InputVector[i]);
            }
        }

        double aux = min_component*MaxTol;

        if( aux < (max_component*MinTol) )
        {
            aux = max_component*MinTol;
        }

        //PerturbationVector
        for( unsigned int i=0; i<VSize; i++ )
        {
            if( fabs(InputVector[i]) > 1.0e-20 ) // different from zero
            {
                rPerturbationVector[i] = InputVector[i]*MaxTol;
            }
            else if( InputVector[i] >= 0.0 )
            {
                rPerturbationVector[i] = aux;
            }
            else
            {
                rPerturbationVector[i] = -aux;
            }
        }
    }


}; /* Class IsotropicDamageUtilities */
} /* namespace Kratos.*/

#endif /* KRATOS_ISOTROPIC_DAMAGE_UTILITIES_H_INCLUDED defined */
