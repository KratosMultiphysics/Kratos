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

    double  PropertiesProxy::GetDensity()                                                     { return *mDensity;                           }
    double* PropertiesProxy::pGetDensity()                                                    { return  mDensity;                           }
    void    PropertiesProxy::SetDensityFromProperties(double* density)                        { mDensity = density;                         }

    int     PropertiesProxy::GetParticleMaterial()                                            { return *mParticleMaterial;                  }
    int*    PropertiesProxy::pGetParticleMaterial()                                           { return  mParticleMaterial;                  }
    void    PropertiesProxy::SetParticleMaterialFromProperties(int* particle_material)        { mParticleMaterial = particle_material;      }

    PropertiesProxy PropertiesProxy::operator= (PropertiesProxy props) {

        mId                         = props.GetId();
        mYoung                      = props.pGetYoung();
        mPoisson                    = props.pGetPoisson();
        mDensity                    = props.pGetDensity();
        mParticleMaterial           = props.pGetParticleMaterial();   

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

            aux_pointer = &( props_it->GetValue(PARTICLE_DENSITY) );
            vector_of_proxies[properties_counter].SetDensityFromProperties(aux_pointer);

            int* int_aux_pointer = &( props_it->GetValue(PARTICLE_MATERIAL) );
            vector_of_proxies[properties_counter].SetParticleMaterialFromProperties(int_aux_pointer);        

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
