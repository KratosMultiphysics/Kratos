//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2014-07-28 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
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
        PropertiesProxy(){} //only used by serializer!
      
        unsigned int     GetId()                                                          { return mId;                                 }
        void    SetId(int id)                                                    { mId = id;                                   }
       
        double  GetYoung()                                                       { return *mYoung;                             }
        double* pGetYoung()                                                      { return  mYoung;                             }
        void    SetYoungFromProperties(double* young)                            { mYoung = young;                             }  
      
        double  GetPoisson()                                                     { return *mPoisson;                           }
        double* pGetPoisson()                                                    { return  mPoisson;                           }
        void    SetPoissonFromProperties(double* poisson)                        { mPoisson = poisson;                         }                    
            
        double  GetRollingFriction()                                             { return *mRollingFriction;                   }      
        double* pGetRollingFriction()                                            { return  mRollingFriction;                   }  
        void    SetRollingFrictionFromProperties(double* rolling_friction)       { mRollingFriction = rolling_friction;        }      
      
        double  GetTgOfFrictionAngle()                                           { return *mTgOfFrictionAngle;                 } 
        double* pGetTgOfFrictionAngle()                                          { return  mTgOfFrictionAngle;                 } 
        void    SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { mTgOfFrictionAngle = tg_of_friction_angle;  }
      
        double  GetLnOfRestitCoeff()                                             { return *mLnOfRestitCoeff;                   } 
        double* pGetLnOfRestitCoeff()                                            { return  mLnOfRestitCoeff;                   } 
        void    SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { mLnOfRestitCoeff = ln_of_restit_coeff;      }  
      
        double  GetDensity()                                                     { return *mDensity;                           }
        double* pGetDensity()                                                    { return  mDensity;                           }
        void    SetDensityFromProperties(double* density)                        { mDensity = density;                         }  
      
        int     GetParticleMaterial()                                            { return *mParticleMaterial;                  }
        int*    pGetParticleMaterial()                                           { return  mParticleMaterial;                  }
        void    SetParticleMaterialFromProperties(int* particle_material)        { mParticleMaterial = particle_material;      }  
      
        PropertiesProxy operator= (PropertiesProxy props) {
          
            mId                = props.GetId();
            mYoung             = props.pGetYoung();
            mPoisson           = props.pGetPoisson();
            mRollingFriction   = props.pGetRollingFriction();
            mTgOfFrictionAngle = props.pGetTgOfFrictionAngle();
            mLnOfRestitCoeff   = props.pGetLnOfRestitCoeff();
            mDensity           = props.pGetDensity();
            mParticleMaterial  = props.pGetParticleMaterial();
                  
            return *this;
        } 
                
           
    private:
        unsigned int     mId;
        double* mYoung;
        double* mPoisson;
        double* mRollingFriction;
        double* mTgOfFrictionAngle;
        double* mLnOfRestitCoeff;
        double* mDensity;
        int*    mParticleMaterial;
      
      
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {

            rSerializer.save("mId",mId);
            /*rSerializer.save("mYoung",mYoung);
            rSerializer.save("mPoisson",mPoisson);
            rSerializer.save("mRollingFriction",mRollingFriction);
            rSerializer.save("mTgOfFrictionAngle",mTgOfFrictionAngle);
            rSerializer.save("mLnOfRestitCoeff",mLnOfRestitCoeff);
            rSerializer.save("mDensity",mDensity);
            rSerializer.save("mParticleMaterial",mParticleMaterial);*/
        }

        virtual void load(Serializer& rSerializer) {
            
            rSerializer.load("mId",mId);
            /*rSerializer.load("mYoung",mYoung);
            rSerializer.load("mPoisson",mPoisson);
            rSerializer.load("mRollingFriction",mRollingFriction);
            rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
            rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);
            rSerializer.load("mDensity",mDensity);
            rSerializer.load("mParticleMaterial",mParticleMaterial);  */        
        }
      
    };
    
    inline void AddPropertiesProxiesFromModelPartProperties(std::vector<PropertiesProxy>& vector_of_proxies, ModelPart& rModelPart, int& properties_counter){
        
          typedef PointerVectorSet<Properties, IndexedObject>               PropertiesContainerType;
          typedef PropertiesContainerType::iterator                PropertiesIterator;
                    
          for (PropertiesIterator props_it = rModelPart.GetMesh(0).PropertiesBegin(); props_it!= rModelPart.GetMesh(0).PropertiesEnd();   props_it++ ) {
              
              vector_of_proxies[properties_counter].SetId( props_it->GetId() );

              double* aux_pointer = &( props_it->GetValue(YOUNG_MODULUS) );
              vector_of_proxies[properties_counter].SetYoungFromProperties( aux_pointer );
              
              aux_pointer = &( props_it->GetValue(POISSON_RATIO) );
              vector_of_proxies[properties_counter].SetPoissonFromProperties(aux_pointer);
                                          
              aux_pointer = &( props_it->GetValue(ROLLING_FRICTION) );
              vector_of_proxies[properties_counter].SetRollingFrictionFromProperties(aux_pointer);
              
              //MA: I commented out the following part because the Initialization of the Proxies is done before the initialization of the inlet (where ProcessInfo for that ModelPart is set)
              /*if ( rModelPart.GetProcessInfo()[ROLLING_FRICTION_OPTION] )  {
                aux_pointer = &( props_it->GetValue(ROLLING_FRICTION) );
                vector_of_proxies[properties_counter].SetRollingFrictionFromProperties(aux_pointer);
              }
              else {
                vector_of_proxies[properties_counter].SetRollingFrictionFromProperties(NULL);
              }*/
              
              aux_pointer = &( props_it->GetValue(PARTICLE_FRICTION) );
              vector_of_proxies[properties_counter].SetTgOfFrictionAngleFromProperties(aux_pointer);
              
              aux_pointer = &( props_it->GetValue(LN_OF_RESTITUTION_COEFF) );
              vector_of_proxies[properties_counter].SetLnOfRestitCoeffFromProperties(aux_pointer);
              
              aux_pointer = &( props_it->GetValue(PARTICLE_DENSITY) );
              vector_of_proxies[properties_counter].SetDensityFromProperties(aux_pointer);
              
              int* int_aux_pointer = &( props_it->GetValue(PARTICLE_MATERIAL) );
              vector_of_proxies[properties_counter].SetParticleMaterialFromProperties(int_aux_pointer);
                                         
              properties_counter++;
                            
          }      
     }    
    
    inline void CreatePropertiesProxies(std::vector<PropertiesProxy>& vector_of_proxies, ModelPart& balls_mp, ModelPart& inlet_mp, ModelPart& clusters_mp){
          KRATOS_TRY
          
          vector_of_proxies.clear();    
          vector_of_proxies.resize( balls_mp.NumberOfProperties() + inlet_mp.NumberOfProperties() + clusters_mp.NumberOfProperties() );
          int properties_counter = 0;
          AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, balls_mp,    properties_counter);          
          AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, inlet_mp,    properties_counter);           
          AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, clusters_mp, properties_counter);                    
          
          return;          
          KRATOS_CATCH("")
    }
    
    inline void CreatePropertiesProxies(std::vector<PropertiesProxy>& vector_of_proxies, ModelPart& r_model_part){
          KRATOS_TRY
          
          vector_of_proxies.clear();    
          vector_of_proxies.resize( r_model_part.NumberOfProperties() );
          int properties_counter = 0;
          AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, r_model_part,    properties_counter);          
          
          return;          
          KRATOS_CATCH("")
    }
    
} // namespace Kratos.

#endif // PROPERTIES_PROXIES_H_INCLUDED defined
