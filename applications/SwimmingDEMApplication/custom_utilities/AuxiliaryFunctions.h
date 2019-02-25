/* 
 * File:   GeometryFunctions.h
 * Author: msantasusana
 *
 * Created on 21 de mayo de 2012, 19:40
 */

#ifndef _SWIMMING_DEM_AUXILIARY_FUNCTIONS_H
#define	_SWIMMING_DEM_AUXILIARY_FUNCTIONS_H


#include <cmath>


namespace Kratos
{
    namespace AuxiliaryFunctions
    {
        static inline int swimming_DEM_floor(double x)
        {
          int i = (int)x; /* truncate */
          return i - ( i > x ); /* convert trunc to floor */
        }
    }
}
#endif	/* _SWIMMING_DEM_AUXILIARY_FUNCTIONS_H */

