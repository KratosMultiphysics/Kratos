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

//#include "../kratos/includes/define.h"
//#include "../custom_elements/discrete_element.h"
//#include "../DEM_application_variables.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) PropertiesProxy {
            
    public:
        
        KRATOS_CLASS_POINTER_DEFINITION(PropertiesProxy);
        
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
        
        double  GetRollingFrictionWithWalls();
        double* pGetRollingFrictionWithWalls();
        void    SetRollingFrictionWithWallsFromProperties(double* rolling_friction_with_walls);
      
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
        
        // Dependent Friction    
        double  GetParticleContactRadius();
        double* pGetParticleContactRadius();
        void    SetParticleContactRadiusFromProperties(double* particle_contact_radius);
    
        double  GetParticleMaxStress();
        double* pGetParticleMaxStress();
        void    SetParticleMaxStressFromProperties(double* particle_max_stress);
    
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
        double* mRollingFrictionWithWalls;
        double* mTgOfFrictionAngle;
        double* mCoefficientOfRestitution;
        double* mLnOfRestitCoeff;
        double* mDensity;
        int*    mParticleMaterial;
        double* mParticleCohesion;
        double* mParticleKNormal;
        double* mParticleKTangential;
        //Dependent Friction
        double* mParticleContactRadius;
        double* mParticleMaxStress;
        double* mParticleGamma;        
        double* mContactSigmaMin; 
        double* mContactTauZero;
        double* mContactInternalFricc;
                
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const;

        virtual void load(Serializer& rSerializer);       
    }; // class PropertiesProxy

    
    inline std::ostream & operator<<( std::ostream& rOut, const PropertiesProxy& rTheProxies){
            rOut << "";
            return rOut;
    }
    

    class KRATOS_API(DEM_APPLICATION) PropertiesProxiesManager {
        
    public:
        KRATOS_CLASS_POINTER_DEFINITION(PropertiesProxiesManager);
        
        void AddPropertiesProxiesFromModelPartProperties(std::vector<PropertiesProxy>& vector_of_proxies,
                                                         ModelPart& rModelPart,
                                                         int& properties_counter);  
    
        void CreatePropertiesProxies(ModelPart& balls_mp,
                                     ModelPart& inlet_mp,
                                     ModelPart& clusters_mp);
    
        void CreatePropertiesProxies(ModelPart& r_model_part);
        
        std::vector<PropertiesProxy>& GetPropertiesProxies(ModelPart& r_model_part);
    }; // class PropertiesProxiesManager
    
    
} // namespace Kratos

#endif // PROPERTIES_PROXIES_H_INCLUDED defined
