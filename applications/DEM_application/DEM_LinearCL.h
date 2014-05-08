
#if !defined(DEM_LINEARCL_H_INCLUDED)
#define  DEM_LINEARCL_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"

namespace Kratos
{

class DEM_LinearCL:DEMConstitutiveLaw
{
public:

    void  CalculateContactForces(double LocalElasticContactForce[3],double indentation,SphericParticle *neighbour_iterator);
}


} /* namespace Kratos.*/
#endif /* DEM_LINEARCL_H_INCLUDED  defined */
