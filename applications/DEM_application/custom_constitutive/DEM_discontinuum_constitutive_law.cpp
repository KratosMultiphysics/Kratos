#include "DEM_application.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos
{
  
    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw()
    {
        //std::cout << " DEMDiscontinuumConstitutiveLaw constructor..." << std::endl;

    } // Class DEMDiscontinuumConstitutiveLaw
    
    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw( const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw){
        //std::cout << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    }

//    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw( const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw){
//        //std::cout << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
//    }


    void DEMDiscontinuumConstitutiveLaw::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {}


    void DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const
    {
      //std::cout << " Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
      pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER,this->Clone());
    }


    DEMDiscontinuumConstitutiveLaw::Pointer DEMDiscontinuumConstitutiveLaw::Clone() const
    {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }
    
    DEMDiscontinuumConstitutiveLaw::~DEMDiscontinuumConstitutiveLaw()
    { 
        //std::cout << "Law destructor..." ; 
    }
    
    
    void DEMDiscontinuumConstitutiveLaw::CalculateContactForces(double LocalElasticContactForce[3],
                                                                double indentation,
                                                                double kn_el,
                                                                double LocalDeltDisp[3],
                                                                double kt_el,
                                                                int& neighbour_failure_id,
                                                                double equiv_tg_of_fri_ang){}


     void DEMDiscontinuumConstitutiveLaw::PlasticityAndDamage(double LocalElasticContactForce[3],
                                                              double kn_el,
                                                              double equiv_young,
                                                              double indentation,
                                                              double calculation_area,
                                                              double radius_sum_i,
                                                              double& failure_criterion_state,
                                                              double& acumulated_damage,
                                                              int i_neighbour_count,
                                                              int mapping_new_cont,
                                                              int mapping_new_ini,
                                                              int time_steps){}



     void DEMDiscontinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                double indentation,
                                                                double equiv_visco_damp_coeff_normal,
                                                                double equiv_visco_damp_coeff_tangential,
                                                                bool sliding,
                                                                int mDampType)
     {

         //*** The compbrobation is component-wise since localContactForce and RelVel have in principle no relationship.
         // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
         // But in oposite direction the visco damping can't overpass the force...

         if (mDampType > 0){

             if (mDampType == 11 || mDampType == 10){
                 ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal * LocalRelVel[2];
}

             if (sliding == false && (mDampType == 1 || mDampType == 11)){ //only applied when no sliding to help to the regularized friction law or the spring convergence
                 ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
                 ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
             }
         }
     }


 //   void DEMDiscontinuumConstitutiveLaw::CalculateContactForces(double LocalElasticContactForce[3],  double indentation){

      } // KRATOS
