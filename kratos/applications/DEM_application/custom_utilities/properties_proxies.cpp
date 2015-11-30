//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2014-07-28 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes

#include "DEM_application.h"
#include "properties_proxies.h"

namespace Kratos {

    PropertiesProxy::PropertiesProxy(){} //only used by serializer!
    
    unsigned int PropertiesProxy::GetId()                                                     { return mId;                                 }
    void    PropertiesProxy::SetId(int id)                                                    { mId = id;                                   }
       
    double  PropertiesProxy::GetYoung()                                                       { return *mYoung;                             }
    double* PropertiesProxy::pGetYoung()                                                      { return  mYoung;                             }
    void    PropertiesProxy::SetYoungFromProperties(double* young)                            { mYoung = young;                             }  
      
    double  PropertiesProxy::GetPoisson()                                                     { return *mPoisson;                           }
    double* PropertiesProxy::pGetPoisson()                                                    { return  mPoisson;                           }
    void    PropertiesProxy::SetPoissonFromProperties(double* poisson)                        { mPoisson = poisson;                         }                    
            
    double  PropertiesProxy::GetRollingFriction()                                             { return *mRollingFriction;                   }      
    double* PropertiesProxy::pGetRollingFriction()                                            { return  mRollingFriction;                   }  
    void    PropertiesProxy::SetRollingFrictionFromProperties(double* rolling_friction)       { mRollingFriction = rolling_friction;        }      
      
    double  PropertiesProxy::GetTgOfFrictionAngle()                                           { return *mTgOfFrictionAngle;                 } 
    double* PropertiesProxy::pGetTgOfFrictionAngle()                                          { return  mTgOfFrictionAngle;                 } 
    void    PropertiesProxy::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { mTgOfFrictionAngle = tg_of_friction_angle;  }
      
    double  PropertiesProxy::GetCoefficientOfRestitution()                                             { return *mCoefficientOfRestitution;                   } 
    double* PropertiesProxy::pGetCoefficientOfRestitution()                                            { return  mCoefficientOfRestitution;                   } 
    void    PropertiesProxy::SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution)     { mCoefficientOfRestitution = coefficient_of_restitution;      }  
    
    double  PropertiesProxy::GetLnOfRestitCoeff()                                             { return *mLnOfRestitCoeff;                   } 
    double* PropertiesProxy::pGetLnOfRestitCoeff()                                            { return  mLnOfRestitCoeff;                   } 
    void    PropertiesProxy::SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { mLnOfRestitCoeff = ln_of_restit_coeff;      }
 
    double  PropertiesProxy::GetDensity()                                                     { return *mDensity;                           }
    double* PropertiesProxy::pGetDensity()                                                    { return  mDensity;                           }
    void    PropertiesProxy::SetDensityFromProperties(double* density)                        { mDensity = density;                         }  
      
    int     PropertiesProxy::GetParticleMaterial()                                            { return *mParticleMaterial;                  }
    int*    PropertiesProxy::pGetParticleMaterial()                                           { return  mParticleMaterial;                  }
    void    PropertiesProxy::SetParticleMaterialFromProperties(int* particle_material)        { mParticleMaterial = particle_material;      }  
        
    double  PropertiesProxy::GetParticleCohesion()                                            { return *mParticleCohesion;                  }
    double* PropertiesProxy::pGetParticleCohesion()                                           { return  mParticleCohesion;                  }
    void    PropertiesProxy::SetParticleCohesionFromProperties(double* particle_cohesion)     { mParticleCohesion = particle_cohesion;      }
        
    PropertiesProxy PropertiesProxy::operator= (PropertiesProxy props) {
          
        mId                       = props.GetId();
        mYoung                    = props.pGetYoung();
        mPoisson                  = props.pGetPoisson();
        mRollingFriction          = props.pGetRollingFriction();
        mTgOfFrictionAngle        = props.pGetTgOfFrictionAngle();
        mCoefficientOfRestitution = props.pGetCoefficientOfRestitution();
        mLnOfRestitCoeff          = props.pGetLnOfRestitCoeff();
        mDensity                  = props.pGetDensity();
        mParticleMaterial         = props.pGetParticleMaterial();
        mParticleCohesion         = props.pGetParticleCohesion();
                       
        return *this;
    } 
                        
    void PropertiesProxy::save(Serializer& rSerializer) const {

        rSerializer.save("mId",mId);
        /*rSerializer.save("mYoung",mYoung);
        rSerializer.save("mPoisson",mPoisson);
        rSerializer.save("mRollingFriction",mRollingFriction);
        rSerializer.save("mTgOfFrictionAngle",mTgOfFrictionAngle);
        rSerializer.save("mLnOfRestitCoeff",mLnOfRestitCoeff);
        rSerializer.save("mDensity",mDensity);
        rSerializer.save("mParticleMaterial",mParticleMaterial);*/
        rSerializer.save("mParticleCohesion",mParticleCohesion);
    }

    void PropertiesProxy::load(Serializer& rSerializer) {
            
        rSerializer.load("mId",mId);
        /*rSerializer.load("mYoung",mYoung);
        rSerializer.load("mPoisson",mPoisson);
        rSerializer.load("mRollingFriction",mRollingFriction);
        rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
        rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);
        rSerializer.load("mDensity",mDensity);
        rSerializer.load("mParticleMaterial",mParticleMaterial);
        rSerializer.load("mParticleCohesion",mParticleCohesion);
        rSerializer.load("mParticleMaterial",mParticleMaterial);*/        
    }
      
    void AddPropertiesProxiesFromModelPartProperties(std::vector<PropertiesProxy>& vector_of_proxies,
                                                                      ModelPart& rModelPart,
                                                                      int& properties_counter) {
        
        typedef PointerVectorSet<Properties, IndexedObject>   PropertiesContainerType;
        typedef PropertiesContainerType::iterator                  PropertiesIterator;
                    
        for (PropertiesIterator props_it = rModelPart.GetMesh(0).PropertiesBegin(); props_it!= rModelPart.GetMesh(0).PropertiesEnd(); props_it++ ) {
              
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
              
            aux_pointer = &( props_it->GetValue(COEFFICIENT_OF_RESTITUTION) );
            vector_of_proxies[properties_counter].SetCoefficientOfRestitutionFromProperties(aux_pointer);

            //aux_pointer = &( props_it->GetValue(LN_OF_RESTITUTION_COEFF) );
            //vector_of_proxies[properties_counter].SetLnOfRestitCoeffFromProperties(aux_pointer);
             
            aux_pointer = &( props_it->GetValue(PARTICLE_DENSITY) );
            vector_of_proxies[properties_counter].SetDensityFromProperties(aux_pointer);
              
            int* int_aux_pointer = &( props_it->GetValue(PARTICLE_MATERIAL) );
            vector_of_proxies[properties_counter].SetParticleMaterialFromProperties(int_aux_pointer);
              
            aux_pointer = &( props_it->GetValue(PARTICLE_COHESION) );
            vector_of_proxies[properties_counter].SetParticleCohesionFromProperties(aux_pointer);
                                         
            properties_counter++;                         
        }      
    }    
    
    void CreatePropertiesProxies(std::vector<PropertiesProxy>& vector_of_proxies,
                                                  ModelPart& balls_mp,
                                                  ModelPart& inlet_mp,
                                                  ModelPart& clusters_mp) {
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
    
    void CreatePropertiesProxies(std::vector<PropertiesProxy>& vector_of_proxies, ModelPart& r_model_part) {

        KRATOS_TRY
          
        vector_of_proxies.clear();    
        vector_of_proxies.resize( r_model_part.NumberOfProperties() );
        int properties_counter = 0;
        AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, r_model_part, properties_counter);          
          
        return;          
    
        KRATOS_CATCH("")
    }
    
} // Namespace Kratos
