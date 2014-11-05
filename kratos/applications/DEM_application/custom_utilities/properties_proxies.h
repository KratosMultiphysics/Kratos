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

#include "../kratos/includes/define.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_constitutive/DEM_continuum_constitutive_law.h"

namespace Kratos {

    class PropertiesProxy {
            
    public:
      
        int     GetId()                                                          { return mId;                                 }
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
        
        void CreateAllPropertiesProxies(ModelPart& rBallsModelPart, ModelPart& rInletModelPart, ModelPart& rClustersModelPart){
            
            
            
        }
           
    private:
        int     mId;
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
            rSerializer.save("mYoung",mYoung);
            rSerializer.save("mPoisson",mPoisson);
            rSerializer.save("mRollingFriction",mRollingFriction);
            rSerializer.save("mTgOfFrictionAngle",mTgOfFrictionAngle);
            rSerializer.save("mLnOfRestitCoeff",mLnOfRestitCoeff);
            rSerializer.save("mDensity",mDensity);
            rSerializer.save("mParticleMaterial",mParticleMaterial);
        }

        virtual void load(Serializer& rSerializer) {
            
            rSerializer.load("mId",mId);
            rSerializer.load("mYoung",mYoung);
            rSerializer.load("mPoisson",mPoisson);
            rSerializer.load("mRollingFriction",mRollingFriction);
            rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
            rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);
            rSerializer.load("mDensity",mDensity);
            rSerializer.load("mParticleMaterial",mParticleMaterial);          
        }
      
    };
    
} // namespace Kratos.

#endif // PROPERTIES_PROXIES_H_INCLUDED defined