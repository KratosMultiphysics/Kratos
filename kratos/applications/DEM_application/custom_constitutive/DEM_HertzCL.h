
#if !defined(DEM_HERTZCL_H_INCLUDED)
#define  DEM_HERTZCL_H_INCLUDED

/* Project includes */
#include "DEM_constitutive_law.h"

namespace Kratos
{

class DEM_HertzCL:DEMDiscontinuumConstitutiveLaw
{
public:

    void  CalculateContactForces(double LocalElasticContactForce[3],double indentation,SphericParticle *neighbour_iterator);
};


} /* namespace Kratos.*/
#endif /* DEM_HERTZCL_H_INCLUDED  defined */
