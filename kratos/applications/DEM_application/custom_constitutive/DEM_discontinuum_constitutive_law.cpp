#include "DEM_application.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos
{
  
    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw()
    {
        std::cout << " DEMDiscontinuumConstitutiveLaw constructor.." << std::endl;

    } // Class DEMDiscontinuumConstitutiveLaw

    int DEMDiscontinuumConstitutiveLaw::Initialize()

    {
        DEMDiscontinuumConstitutiveLaw test;
        std::cout << " DEMDiscontinuumConstitutiveLaw print on screen..." << std::endl;
        return 0;
    }


    void DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const
    {
      std::cout << " Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
      pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER,this->Clone());
    }


    DEMDiscontinuumConstitutiveLaw::Pointer DEMDiscontinuumConstitutiveLaw::Clone() const
    {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }
    
    DEMDiscontinuumConstitutiveLaw::~DEMDiscontinuumConstitutiveLaw()
    { std::cout << "Law destructor..." ; }
    
    
     void DEMDiscontinuumConstitutiveLaw::PlasticityAndDamage(double LocalElasticContactForce[3], double kn_el, double equiv_young, double indentation, double      calculation_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps)
     {}



    void DEMDiscontinuumConstitutiveLaw::CalculateContactForces(double LocalElasticContactForce[3],  double indentation){
    //  linear compression & tension


      } //Calculate Contact Forces with Linear CL

}
