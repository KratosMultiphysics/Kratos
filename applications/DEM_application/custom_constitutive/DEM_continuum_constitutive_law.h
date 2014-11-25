
#if !defined(DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "containers/flags.h"

#include "custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"
namespace Kratos {

        /**
     * Base class of constitutive laws.
     */
    class /*__declspec( dllexport )*/ DEMContinuumConstitutiveLaw : public Flags {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMContinuumConstitutiveLaw);

        DEMContinuumConstitutiveLaw();

        DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        virtual ~DEMContinuumConstitutiveLaw();

        virtual DEMContinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateContactForces(double mRadius,
                double mSqrtOfRealMass,
                double other_radius,
                double otherSqrtMass,
                double distance,
                double initial_delta,
                int& neighbour_failure_id,
                ProcessInfo& rCurrentProcessInfo,
                PropertiesProxy *myProperties,
                PropertiesProxy *neighbourProperties,
                int mapping_new_ini,
                int mapping_new_cont,
                unsigned int i_neighbour_count,
                double LocalElasticContactForce[3],
                double ViscoDampingLocalContactForce[3],
                double LocalDeltDisp[3],
                Vector mcont_ini_neigh_area,
                array_1d<double, 4 > &mHistory_mapping_new_cont,
                double mDempack_damping,
                int mDampType,
                int mIniNeighbourFailureId_mapping_new_ini,
                double LocalCoordSystem[3][3],
                double RelVel[3]); //FF


        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding,
                int mDampType);




    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

    KRATOS_DEFINE_VARIABLE(DEMContinuumConstitutiveLaw::Pointer, DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER)

    inline void EvaluateFailureCriteria(
            const double contact_sigma,
            const double contact_tau,
            double& failure_criterion_state,
            bool& sliding,
            const int FailureCriterionOption,
            const double TauZero,
            const double TanContactInternalFriccion,
            const double SinContactInternalFriccion,
            const double CosContactInternalFriccion,
            int& NeighbourFailureId_i_neighbour_count, 
            int& IniNeighbourFailureId_mapping_new_ini, 
            const double TensionLimit) {

        //(1) MOHR-COULOMB FAILURE: (we don't consider rotational spring!!!!! here) need to be thought.
        if (FailureCriterionOption == 1) {//MOHR-COULOMB                    
            double sigma_max, sigma_min;

            if (contact_sigma >= 0) {
                sigma_max = contact_sigma;
                sigma_min = 0;
            } else {
                sigma_max = 0;
                sigma_min = contact_sigma;
            }

            //change into principal stresses
            double centre = 0.5 * (sigma_max + sigma_min);
            double radius = sqrt((sigma_max - centre)*(sigma_max - centre) + contact_tau * contact_tau);

            double sigma_I = centre + radius;
            double sigma_II = centre - radius;

            // Check:
            double distance_to_failure = (TauZero / (TanContactInternalFriccion) + centre) * SinContactInternalFriccion;
            failure_criterion_state = radius / distance_to_failure;

            if (sigma_I - sigma_II >= 2 * TauZero * CosContactInternalFriccion + (sigma_I + sigma_II) * SinContactInternalFriccion) {
                //breaks
                NeighbourFailureId_i_neighbour_count = 5;                       //mohr coulomb   
                IniNeighbourFailureId_mapping_new_ini = 5;
                failure_criterion_state = 1.0;
                sliding = true;
            }
        } //MOHR-COULOMB

        ///(2) UNCOUPLED FRACTURE
        if (FailureCriterionOption == 2) {//UNCOUPLED FRACTURE and DEMPACK        
            //double mTauZero = 0.5*sqrt(mCompressionLimit*mTensionLimit); 

            if (contact_sigma >= 0) {
                double tau_strength = TauZero + TanContactInternalFriccion*contact_sigma;
                failure_criterion_state = contact_tau / tau_strength;

                if (contact_tau > tau_strength) {
                    NeighbourFailureId_i_neighbour_count = 2; //shear in compression
                    IniNeighbourFailureId_mapping_new_ini = 2;
                    failure_criterion_state = 1.0;
                    sliding = true;
                }
            }//positive sigmas

            else { //negative sigmas            
                double tau_strength = TauZero;
                failure_criterion_state = GeometryFunctions::max(contact_tau / tau_strength, -contact_sigma / TensionLimit);

                if (contact_tau > tau_strength) {
                    NeighbourFailureId_i_neighbour_count = 3; //shear failure tension
                    IniNeighbourFailureId_mapping_new_ini = 3;
                    sliding = true;
                    failure_criterion_state = 1.0;                                     
                }
                
            } //negative values of sigma              

        } //UNCOUPLED FRACTURE
        if (failure_criterion_state > 1.0) failure_criterion_state = 1.0;

    }

inline void PlasticityAndDamage1D ( double LocalElasticContactForce[3], 
                                    const double kn_el, 
                                    double equiv_young, 
                                    double indentation, 
                                    double calculation_area, 
                                    double radius_sum_i, 
                                    double& failure_criterion_state, 
                                    double& acumulated_damage, 
                                    int i_neighbour_count, 
                                    int mapping_new_cont, 
                                    int mapping_new_ini,
                                    const double mN1,
                                    const double mN2,
                                    const double mN3,
                                    const double mYoungPlastic,
                                    const double mPlasticityLimit,
                                    const double mC1,
                                    const double mC2,
                                    const double mC3,
                                    const double mTensionLimit,
                                    const double mDamageMaxDisplacementFactor,
                                    array_1d <double, 4> &mHistory_mapping_new_cont,
                                    int &mNeighbourFailureId_i_neighbour_count,
                                    int &mIniNeighbourFailureId_mapping_new_ini,
                                    int time_steps) {   

      double kn_b = kn_el / mN1;
      double kn_c = kn_el / mN2;
      double kn_d = kn_el / mN3;
      double kp_el = mYoungPlastic/equiv_young * kn_el;
      double Yields_el = mPlasticityLimit * calculation_area;
        
      double Ncstr1_el = mC1 * calculation_area;
      double Ncstr2_el = mC2 * calculation_area;
      double Ncstr3_el = mC3 * calculation_area;
      double Ntstr_el  = mTensionLimit * calculation_area;
      double u_max = mHistory_mapping_new_cont[0];
      
      double& fn = LocalElasticContactForce[2];                                 //[2] means 'normal' contact force
                
      
      if( indentation >= 0.0 ){ //COMPRESSION

          fn = kn_el * indentation;         
          double u_ela1 = Ncstr1_el/kn_el;;
          double u_ela2 = u_ela1 + (Ncstr2_el-Ncstr1_el)/(kn_b);
          double u_ela3 = u_ela2 + (Ncstr3_el-Ncstr2_el)/(kn_c);

          if ( ( indentation > u_max ) || ( time_steps <= 1) ){                 //maximum historical intentation OR first step  MSIMSI 0

            mHistory_mapping_new_cont[0]  = indentation;                        // Guarda el treshold del màxim desplaçament
    
            if (indentation > u_ela3) {                                         //4rt tram
                                                                                
              fn = Ncstr3_el + ( indentation - u_ela3 )*kn_d;
              mHistory_mapping_new_cont[3] = kn_d/kn_el;
            }
            else if (indentation > u_ela2) {   
                                                                                //3r tram
              fn = Ncstr2_el + ( indentation - u_ela2 )*kn_c;
              mHistory_mapping_new_cont[3] = kn_c/kn_el; 
            }
            else{
                
              if( indentation > u_ela1){                                        //2n tram
              
                fn = Ncstr1_el + (indentation - u_ela1)*kn_b;
                mHistory_mapping_new_cont[3] = kn_b/kn_el;
              } 
            }
          mHistory_mapping_new_cont[1] = fn;                                    //actualitzar la força màxima a compressió.
          }
          
          else{                                                                 //Per sota del màxim.
   
              if(mHistory_mapping_new_cont[1] > 0.0) {                          //Màxim en compressió. 
              
                  double u_plas;                                                //MSIMSI 2 akesta operació de saber quant val la u_plastica es fa cada pas de temps i en realitat es fixe sempre.
                  if(Yields_el <= Ncstr1_el){                                   //si el punt de plastificació està en la primera rama elastica.
                  
                      u_plas = Yields_el/kn_el;
                  }
                  else{
                    
                    if(Yields_el <= Ncstr2_el)                                  //si està en la segona...
                    {
                        u_plas = u_ela1 + (Yields_el-Ncstr1_el)/(kn_b);
                    }
                     else if(Yields_el <= Ncstr3_el){                           //si està en la tercera...
                    
                        u_plas = u_ela2 + (Yields_el-Ncstr2_el)/(kn_c);
                    }                    
                    else {                                                      //en la quarta                   
                      u_plas = u_ela3 + (Yields_el-Ncstr3_el)/(kn_d);
                    }                   
                  }                 
                  if ( u_plas < u_max ) {                                       //si nosaltres estem per sota del maxim pero ja estem plastificant 
                  
                    fn = mHistory_mapping_new_cont[1] - kp_el*(u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                    mHistory_mapping_new_cont[3] = kp_el/kn_el;                    
                  }
                  else {                                                        // Esta en zona descarga elastica, ens despreocupem de la plasticitat                 
                    if ( indentation > u_ela3) {                                //en la 4a ramma                   
                      fn = Ncstr3_el + (indentation - u_ela3)*kn_d;                      
                    }                    
                    else if ( indentation > u_ela2) {                           //en la 3a ramma                    
                      fn = Ncstr2_el + (indentation - u_ela2)*kn_c;                      
                    }
                    
                    else{                    
                      if(indentation > u_ela1){                                 //en la 2a rama                     
                        fn = Ncstr1_el + (indentation-u_ela1)*kn_b;
                      }           
                    }           
                  }
              }                                                                 //si tenim precàrrega en compressió.              
          }                                                                     //Per sota del màxim.
        }                                                                       //Compression
      else {                                                                    //tension      
        fn = kn_el * indentation; 
        
        double u1 = Ntstr_el / kn_el;
        double u2 = u1*(1+ mDamageMaxDisplacementFactor);
 
          if(fabs(indentation) > u2) {                                          // FULL DAMAGE 
          
            mNeighbourFailureId_i_neighbour_count = 4;                          //tension failure
            mIniNeighbourFailureId_mapping_new_ini = 4;
            acumulated_damage = 1.0;
            fn = 0.0;
          }
          else{
            if (fabs(indentation) > u1)  {            
              double u_frac = (fabs(indentation) - u1)/(u2 - u1);
                                                                                //failure_criterion_state = fabs(indentation)/u2;
              acumulated_damage = u_frac;
                            
              if (u_frac > mHistory_mapping_new_cont[2])  {              
                mHistory_mapping_new_cont[2] = u_frac;
              }              
            }           
            fn = indentation * kn_el*(1.0 -  mHistory_mapping_new_cont[2]);     // normal adhesive force (gap +)            
          }     
        }                                                                       //Tension    
    }   
}                                                                               /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

