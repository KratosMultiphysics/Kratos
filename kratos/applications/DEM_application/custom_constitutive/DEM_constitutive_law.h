
#if !defined(DEM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"

namespace Kratos
{

/**
 * Base class of constitutive laws.
 */
class DEMConstitutiveLaw
{

public:

   virtual CalculateContactForces(double LocalElasticContactForce[3],double indentation,SphericParticle *neighbour_iterator);

   virtual PlasticityAndDamage(double LocalElasticContactForce[3], double kn_el, double equiv_young, double indentation, double      calculation_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

