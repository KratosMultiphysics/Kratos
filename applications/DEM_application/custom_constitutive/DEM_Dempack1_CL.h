
#if !defined(DEM_DEMPACK1_H_INCLUDED)
#define  DEM_DEMPACK1_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"

namespace Kratos
{

class DEM_Dempack1:DEMContinuumConstitutiveLaw
{
public:

      virtual void CalculateContactForces(double mRadius,double mSqrtOfRealMass, double other_radius,double otherSqrtMass,double distance,double initial_delta,ProcessInfo& rCurrentProcessInfo,PropertiesProxy *myProperties,PropertiesProxy *neighbourProperties,int mapping_new_ini,int mapping_new_cont, int i_neighbour_count,double LocalElasticContactForce[3]);


private:

        //Non-linear
         double mN1;
         double mN2;
         double mN3;
         double mC1;
         double mC2;
         double mC3;



    virtual void PlasticityAndDamage(double LocalElasticContactForce[3], double kn, double equiv_young, double indentation, double corrected_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps);

};
  
  

} /* namespace Kratos.*/
#endif /* DEM_DEMPACK1_H_INCLUDED  defined */
