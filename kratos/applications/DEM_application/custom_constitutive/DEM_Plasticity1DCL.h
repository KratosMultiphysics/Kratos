
#if !defined(DEM_PLASTICITY1D_H_INCLUDED)
#define  DEM_PLASTICITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"

namespace Kratos
{

class DEM_Plasticity1D:DEMConstitutiveLaw
{
public:

    void  CalculateContactForces(double LocalElasticContactForce[3],double indentation,SphericParticle *neighbour_iterator);
}
  

void PlasticityAndDamage1D(double LocalElasticContactForce[3], double kn, double equiv_young, double indentation, double corrected_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps);
  
      

} /* namespace Kratos.*/
#endif /* DEM_PLASTICITY_H_INCLUDED  defined */
