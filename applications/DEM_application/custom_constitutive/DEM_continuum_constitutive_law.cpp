#include "DEM_application.h"
#include "DEM_continuum_constitutive_law.h"

namespace Kratos
{
  
    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw()
    {
        std::cout << " DEMContinuumConstitutiveLaw constructor.." << std::endl;

    } // Class DEMContinuumConstitutiveLaw

    int DEMContinuumConstitutiveLaw::Initialize()

    {
        DEMContinuumConstitutiveLaw test;
        std::cout << " DEMContinuumConstitutiveLaw print on screen..." << std::endl;
        return 0;
    }


    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const
    {
      std::cout << " Assigning DEMContinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
      pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER,this->Clone());
    }


    DEMContinuumConstitutiveLaw::Pointer DEMContinuumConstitutiveLaw::Clone() const
    {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEMContinuumConstitutiveLaw(*this));
        return p_clone;
    }
    
    DEMContinuumConstitutiveLaw::~DEMContinuumConstitutiveLaw()
    { std::cout << "Law destructor..." ; }
    
    
     void DEMContinuumConstitutiveLaw::PlasticityAndDamage(double LocalElasticContactForce[3], double kn_el, double equiv_young, double indentation, double      calculation_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps)
     {}



    void DEMContinuumConstitutiveLaw::CalculateContactForces(double LocalElasticContactForce[3],  double indentation){
    //  linear compression & tension


      } //Calculate Contact Forces with Linear CL

}
