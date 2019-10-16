//
// Author:
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes

#include "properties_proxies.h"
#include "../DEM_application_variables.h"


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

    double  PropertiesProxy::GetRollingFrictionWithWalls()                                    { return *mRollingFrictionWithWalls;          }
    double* PropertiesProxy::pGetRollingFrictionWithWalls()                                   { return  mRollingFrictionWithWalls;          }
    void    PropertiesProxy::SetRollingFrictionWithWallsFromProperties(double* rolling_friction_with_walls) { mRollingFrictionWithWalls = rolling_friction_with_walls; }

    double  PropertiesProxy::GetTgOfFrictionAngle()                                           { return *mTgOfFrictionAngle;                 }
    double* PropertiesProxy::pGetTgOfFrictionAngle()                                          { return  mTgOfFrictionAngle;                 }
    void    PropertiesProxy::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { mTgOfFrictionAngle = tg_of_friction_angle;  }

    double  PropertiesProxy::GetCoefficientOfRestitution()                                    { return *mCoefficientOfRestitution;          }
    double* PropertiesProxy::pGetCoefficientOfRestitution()                                   { return  mCoefficientOfRestitution;          }
    void    PropertiesProxy::SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution) { mCoefficientOfRestitution = coefficient_of_restitution;      }

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

    double  PropertiesProxy::GetAmountOfCohesionFromStress()                                  { return *mAmountOfCohesionFromStress;        }
    double* PropertiesProxy::pGetAmountOfCohesionFromStress()                                 { return  mAmountOfCohesionFromStress;        }
    void    PropertiesProxy::SetAmountOfCohesionFromStressFromProperties(double* amount_of_cohesion_from_stress) { mAmountOfCohesionFromStress = amount_of_cohesion_from_stress; }

    double  PropertiesProxy::GetParticleKNormal()                                             { return *mParticleKNormal;                   }
    double* PropertiesProxy::pGetParticleKNormal()                                            { return  mParticleKNormal;                   }
    void    PropertiesProxy::SetParticleKNormalFromProperties(double* particle_k_normal)      { mParticleKNormal = particle_k_normal;       }

    double  PropertiesProxy::GetParticleKTangential()                                         { return *mParticleKTangential;               }
    double* PropertiesProxy::pGetParticleKTangential()                                        { return  mParticleKTangential;               }
    void    PropertiesProxy::SetParticleKTangentialFromProperties(double* particle_k_tangential) { mParticleKTangential = particle_k_tangential;}

    double  PropertiesProxy::GetParticleConicalDamageContactRadius()                          { return *mParticleConicalDamageContactRadius; }
    double* PropertiesProxy::pGetParticleConicalDamageContactRadius()                         { return  mParticleConicalDamageContactRadius; }
    void    PropertiesProxy::SetParticleConicalDamageContactRadiusFromProperties(double* particle_conical_damage_contact_radius) { mParticleConicalDamageContactRadius = particle_conical_damage_contact_radius; }

    double  PropertiesProxy::GetParticleConicalDamageMaxStress()                               { return *mParticleConicalDamageMaxStress;    }
    double* PropertiesProxy::pGetParticleConicalDamageMaxStress()                              { return  mParticleConicalDamageMaxStress;    }
    void    PropertiesProxy::SetParticleConicalDamageMaxStressFromProperties(double* particle_conical_damage_max_stress) { mParticleConicalDamageMaxStress = particle_conical_damage_max_stress; }

    double  PropertiesProxy::GetParticleConicalDamageGamma()                                   { return *mParticleConicalDamageGamma;        }
    double* PropertiesProxy::pGetParticleConicalDamageGamma()                                  { return  mParticleConicalDamageGamma;        }
    void    PropertiesProxy::SetParticleConicalDamageGammaFromProperties(double* particle_conical_damage_gamma) { mParticleConicalDamageGamma = particle_conical_damage_gamma; }

    double  PropertiesProxy::GetLevelOfFouling()                                              { return *mLevelOfFouling;                     }
    double* PropertiesProxy::pGetLevelOfFouling()                                             { return  mLevelOfFouling;                     }
    void    PropertiesProxy::SetLevelOfFoulingFromProperties(double* level_of_fouling)        { mLevelOfFouling = level_of_fouling;          }

    double  PropertiesProxy::GetContactTauZero()                                              { return *mContactTauZero;                     }
    double* PropertiesProxy::pGetContactTauZero()                                             { return  mContactTauZero;                     }
    void    PropertiesProxy::SetContactTauZeroFromProperties(double* contact_tau_zero)        { mContactTauZero = contact_tau_zero;          }

    double  PropertiesProxy::GetContactInternalFricc()                                        { return *mContactInternalFricc;               }
    double* PropertiesProxy::pGetContactInternalFricc()                                       { return  mContactInternalFricc;               }
    void    PropertiesProxy::SetContactInternalFriccFromProperties(double* contact_internal_fricc) { mContactInternalFricc = contact_internal_fricc; }

    PropertiesProxy PropertiesProxy::operator= (PropertiesProxy props) {

        mId                         = props.GetId();
        mYoung                      = props.pGetYoung();
        mPoisson                    = props.pGetPoisson();
        mRollingFriction            = props.pGetRollingFriction();
        mRollingFrictionWithWalls   = props.pGetRollingFrictionWithWalls();
        mTgOfFrictionAngle          = props.pGetTgOfFrictionAngle();
        mCoefficientOfRestitution   = props.pGetCoefficientOfRestitution();
        mLnOfRestitCoeff            = props.pGetLnOfRestitCoeff();
        mDensity                    = props.pGetDensity();
        mParticleMaterial           = props.pGetParticleMaterial();
        mParticleCohesion           = props.pGetParticleCohesion();
        mAmountOfCohesionFromStress = props.pGetAmountOfCohesionFromStress();
        mParticleKNormal            = props.pGetParticleKNormal();
        mParticleKTangential        = props.pGetParticleKTangential();

        return *this;
    }

    void PropertiesProxy::save(Serializer& rSerializer) const {
    }

    void PropertiesProxy::load(Serializer& rSerializer) {
    }

    void PropertiesProxiesManager::AddPropertiesProxiesFromModelPartProperties(std::vector<PropertiesProxy>& vector_of_proxies,
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

            aux_pointer = &( props_it->GetValue(ROLLING_FRICTION_WITH_WALLS) );
            vector_of_proxies[properties_counter].SetRollingFrictionWithWallsFromProperties(aux_pointer);

            //MA: I commented out the following part because the Initialization of the Proxies is done before the initialization of the inlet (where ProcessInfo for that ModelPart is set)
            /*if ( rModelPart.GetProcessInfo()[ROLLING_FRICTION_OPTION] )  {
            aux_pointer = &( props_it->GetValue(ROLLING_FRICTION) );
            vector_of_proxies[properties_counter].SetRollingFrictionFromProperties(aux_pointer);
            }
            else {
            vector_of_proxies[properties_counter].SetRollingFrictionFromProperties(NULL);
            }*/

            aux_pointer = &( props_it->GetValue(FRICTION) );
            vector_of_proxies[properties_counter].SetTgOfFrictionAngleFromProperties(aux_pointer);

            aux_pointer = &( props_it->GetValue(COEFFICIENT_OF_RESTITUTION) );
            vector_of_proxies[properties_counter].SetCoefficientOfRestitutionFromProperties(aux_pointer);

            aux_pointer = &( props_it->GetValue(PARTICLE_DENSITY) );
            vector_of_proxies[properties_counter].SetDensityFromProperties(aux_pointer);

            int* int_aux_pointer = &( props_it->GetValue(PARTICLE_MATERIAL) );
            vector_of_proxies[properties_counter].SetParticleMaterialFromProperties(int_aux_pointer);

            aux_pointer = &( props_it->GetValue(PARTICLE_COHESION) );
            vector_of_proxies[properties_counter].SetParticleCohesionFromProperties(aux_pointer);

            aux_pointer = &( props_it->GetValue(AMOUNT_OF_COHESION_FROM_STRESS) );
            vector_of_proxies[properties_counter].SetAmountOfCohesionFromStressFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(K_NORMAL));
            vector_of_proxies[properties_counter].SetParticleKNormalFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(K_TANGENTIAL));
            vector_of_proxies[properties_counter].SetParticleKTangentialFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(CONICAL_DAMAGE_CONTACT_RADIUS));
            vector_of_proxies[properties_counter].SetParticleConicalDamageContactRadiusFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(CONICAL_DAMAGE_MAX_STRESS));
            vector_of_proxies[properties_counter].SetParticleConicalDamageMaxStressFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(CONICAL_DAMAGE_GAMMA));
            vector_of_proxies[properties_counter].SetParticleConicalDamageGammaFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(LEVEL_OF_FOULING));
            vector_of_proxies[properties_counter].SetLevelOfFoulingFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(CONTACT_TAU_ZERO));
            vector_of_proxies[properties_counter].SetContactTauZeroFromProperties(aux_pointer);

            aux_pointer = &(props_it->GetValue(CONTACT_INTERNAL_FRICC));
            vector_of_proxies[properties_counter].SetContactInternalFriccFromProperties(aux_pointer);

            properties_counter++;
        }
    }

    void PropertiesProxiesManager::CreatePropertiesProxies(
                                                  ModelPart& balls_mp,
                                                  ModelPart& inlet_mp,
                                                  ModelPart& clusters_mp) {
        KRATOS_TRY

        balls_mp[VECTOR_OF_PROPERTIES_PROXIES] = std::vector<PropertiesProxy>();

        std::vector<PropertiesProxy>& vector_of_proxies = balls_mp[VECTOR_OF_PROPERTIES_PROXIES];
        vector_of_proxies.clear();
        vector_of_proxies.resize( balls_mp.NumberOfProperties() + inlet_mp.NumberOfProperties() + clusters_mp.NumberOfProperties() );
        int properties_counter = 0;
        AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, balls_mp,    properties_counter);
        AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, inlet_mp,    properties_counter);
        AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, clusters_mp, properties_counter);

        return;

        KRATOS_CATCH("")
    }

    void PropertiesProxiesManager::CreatePropertiesProxies(ModelPart& r_model_part) {

        KRATOS_TRY

        r_model_part[VECTOR_OF_PROPERTIES_PROXIES] = std::vector<PropertiesProxy>();

        std::vector<PropertiesProxy>& vector_of_proxies = r_model_part[VECTOR_OF_PROPERTIES_PROXIES];
        vector_of_proxies.clear();
        vector_of_proxies.resize( r_model_part.NumberOfProperties() );
        int properties_counter = 0;
        AddPropertiesProxiesFromModelPartProperties(vector_of_proxies, r_model_part, properties_counter);

        return;

        KRATOS_CATCH("")
    }

    std::vector<PropertiesProxy>& PropertiesProxiesManager::GetPropertiesProxies(ModelPart& r_model_part) {
        return r_model_part[VECTOR_OF_PROPERTIES_PROXIES];
    }

} // Namespace Kratos
