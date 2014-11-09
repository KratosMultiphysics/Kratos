#ifndef KRATOS_BECHMARK_UTILITIES_H
#define KRATOS_BECHMARK_UTILITIES_H

// Project includes
#include "includes/model_part.h"
#include "../applications/DEM_application/custom_utilities/properties_proxies.h"

namespace Kratos
{
class BenchmarkUtils
{
public:

typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
typedef ModelPart::NodesContainerType::iterator     NodeIterator;
typedef PointerVectorSet<Properties, IndexedObject>               PropertiesContainerType;
typedef typename PropertiesContainerType::iterator                PropertiesIterator;

KRATOS_CLASS_POINTER_DEFINITION(BenchmarkUtils);

BenchmarkUtils(){}
/// Calculator

virtual ~BenchmarkUtils(){}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


void ComputeHydrodynamicForces(ModelPart& r_model_part)
{
    mFastProperties.clear();
    mFastProperties.resize(r_model_part.NumberOfProperties());
    int properties_counter = 0;
    AddPropertiesProxiesFromModelPartProperties(r_model_part, properties_counter);

    for (ModelPart::ElementIterator i_elem = r_model_part.ElementsBegin(); i_elem != r_model_part.ElementsEnd(); ++i_elem){
        Kratos::SphericParticle& particle = static_cast<Kratos::SphericParticle&>(*i_elem);
        particle.SetFastProperties(mFastProperties);
        particle.Initialize();
        particle.MemberDeclarationFirstStep(r_model_part.GetProcessInfo());
        particle.GetGeometry()[0].Set(INSIDE, true);
        particle.Set(BLOCKED, false);
        array_1d<double, 3> force;
        array_1d<double, 3> moment;
        array_1d<double, 3> gravity = r_model_part.GetProcessInfo()[GRAVITY];
        particle.ComputeAdditionalForces(force, moment, r_model_part.GetProcessInfo(), gravity);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

std::vector<PropertiesProxy> mFastProperties;

void AddPropertiesProxiesFromModelPartProperties(ModelPart& rModelPart, int& properties_counter){

    PropertiesProxy aux_props;

    for (PropertiesIterator props_it = rModelPart.GetMesh(0).PropertiesBegin(); props_it!= rModelPart.GetMesh(0).PropertiesEnd();   props_it++ ) {

        aux_props.SetId( props_it->GetId() );

        double* aux_pointer = &( props_it->GetValue(YOUNG_MODULUS) );
        aux_props.SetYoungFromProperties( aux_pointer );

        aux_pointer = &( props_it->GetValue(POISSON_RATIO) );
        aux_props.SetPoissonFromProperties(aux_pointer);

        aux_pointer = &( props_it->GetValue(ROLLING_FRICTION) );
        aux_props.SetRollingFrictionFromProperties(aux_pointer);

        aux_pointer = &( props_it->GetValue(PARTICLE_FRICTION) );
        aux_props.SetTgOfFrictionAngleFromProperties(aux_pointer);

        aux_pointer = &( props_it->GetValue(LN_OF_RESTITUTION_COEFF) );
        aux_props.SetLnOfRestitCoeffFromProperties(aux_pointer);

        aux_pointer = &( props_it->GetValue(PARTICLE_DENSITY) );
        aux_props.SetDensityFromProperties(aux_pointer);

        int* int_aux_pointer = &( props_it->GetValue(PARTICLE_MATERIAL) );
        aux_props.SetParticleMaterialFromProperties(int_aux_pointer);
        mFastProperties[properties_counter].SetParticleMaterialFromProperties(int_aux_pointer);

        properties_counter++;

    }

}
};
}
#endif // KRATOS_BECHMARK_UTILITIES_H
