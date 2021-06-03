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

        double  GetDensity();
        double* pGetDensity();
        void    SetDensityFromProperties(double* density);

        int     GetParticleMaterial();
        int*    pGetParticleMaterial();
        void    SetParticleMaterialFromProperties(int* particle_material);       

        PropertiesProxy operator= (PropertiesProxy props);

    private:

        unsigned int mId;
        double* mYoung;
        double* mPoisson;
        double* mDensity;
        int*    mParticleMaterial;

        friend class Serializer;

        void save(Serializer& rSerializer) const;

        void load(Serializer& rSerializer);
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
