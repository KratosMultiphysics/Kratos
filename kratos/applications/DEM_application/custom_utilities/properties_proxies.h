//
// Authors: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

#if !defined(PROPERTIES_PROXIES_H_INCLUDED)
#define PROPERTIES_PROXIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
// Project includes
#include "includes/model_part.h"

#include "../kratos/includes/define.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../DEM_application_variables.h"

namespace Kratos {

    class PropertiesProxy {        
            
    public:
        
        PropertiesProxy(); //only used by serializer!
      
        unsigned int GetId();
        void    SetId(int id);
       
        double  GetYoung();
        double* pGetYoung();
        void    SetYoungFromProperties(double* young);
      
        double  GetPoisson();                                                     
        double* pGetPoisson();                                                   
        void    SetPoissonFromProperties(double* poisson);                                     
            
        double  GetRollingFriction();                                                 
        double* pGetRollingFriction();                                            
        void    SetRollingFrictionFromProperties(double* rolling_friction);        
      
        double  GetTgOfFrictionAngle();                                          
        double* pGetTgOfFrictionAngle();                                          
        void    SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle);
      
        double  GetCoefficientOfRestitution();                                             
        double* pGetCoefficientOfRestitution();                                            
        void    SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution);

        double  GetLnOfRestitCoeff();                                             
        double* pGetLnOfRestitCoeff();                                            
        void    SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff);     
      
        double  GetDensity();                                                    
        double* pGetDensity();                                                   
        void    SetDensityFromProperties(double* density);                       
      
        int     GetParticleMaterial();                                           
        int*    pGetParticleMaterial();                                          
        void    SetParticleMaterialFromProperties(int* particle_material);        
        
        double  GetParticleCohesion();                                            
        double* pGetParticleCohesion();                                           
        void    SetParticleCohesionFromProperties(double* particle_cohesion);
        
        double  GetParticleKNormal();
        double* pGetParticleKNormal();
        void    SetParticleKNormalFromProperties(double* particle_k_normal);

        double  GetParticleKTangential();
        double* pGetParticleKTangential();
        void    SetParticleKTangentialFromProperties(double* particle_k_tangential);
        
        //Conical damage    
        double  GetParticleContactRadius();
        double* pGetParticleContactRadius();
        void    SetParticleContactRadiusFromProperties(double* particle_contact_radius);
    
        double  GetParticleMaxStress();
        double* pGetParticleMaxStress();
        void    SetParticleMaxStressFromProperties(double* particle_max_stress);
    
        double  GetParticleAlpha();
        double* pGetParticleAlpha();
        void    SetParticleAlphaFromProperties(double* particle_alpha);
    
        double  GetParticleGamma();
        double* pGetParticleGamma();
        void    SetParticleGammaFromProperties(double* particle_gamma);
        
        double  GetContactSigmaMin();
        double* pGetContactSigmaMin();
        void    SetContactSigmaMinFromProperties(double* contact_sigma_min);
        
        double  GetContactTauZero();
        double* pGetContactTauZero();
        void    SetContactTauZeroFromProperties(double* contact_tau_zero);
        
        double  GetContactInternalFricc();
        double* pGetContactInternalFricc();
        void    SetContactInternalFriccFromProperties(double* contact_internal_fricc);
        
        PropertiesProxy operator= (PropertiesProxy props);
                       
    private:
        
        unsigned int mId;
        double* mYoung;
        double* mPoisson;
        double* mRollingFriction;
        double* mTgOfFrictionAngle;
        double* mCoefficientOfRestitution;
        double* mLnOfRestitCoeff;
        double* mDensity;
        int*    mParticleMaterial;
        double* mParticleCohesion;
        double* mParticleKNormal;
        double* mParticleKTangential;
        //Conical damage    
        double* mParticleContactRadius;
        double* mParticleMaxStress;
        double* mParticleAlpha;
        double* mParticleGamma;        
        double* mContactSigmaMin; 
        double* mContactTauZero;
        double* mContactInternalFricc;
                
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const;

        virtual void load(Serializer& rSerializer);       
    }; // class PropertiesProxy
    
        void AddPropertiesProxiesFromModelPartProperties(std::vector<PropertiesProxy>& vector_of_proxies,
                                                         ModelPart& rModelPart,
                                                         int& properties_counter);  
    
        void CreatePropertiesProxies(std::vector<PropertiesProxy>& vector_of_proxies,
                                     ModelPart& balls_mp,
                                     ModelPart& inlet_mp,
                                     ModelPart& clusters_mp);
    
        void CreatePropertiesProxies(std::vector<PropertiesProxy>& vector_of_proxies,
                                     ModelPart& r_model_part);
    
    
} // namespace Kratos

#endif // PROPERTIES_PROXIES_H_INCLUDED defined
