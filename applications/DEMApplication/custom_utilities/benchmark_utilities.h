#ifndef KRATOS_BENCHMARK_UTILITIES_H
#define KRATOS_BENCHMARK_UTILITIES_H

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
typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;
typedef PropertiesContainerType::iterator           PropertiesIterator;

KRATOS_CLASS_POINTER_DEFINITION(BenchmarkUtils);

BenchmarkUtils(){}
/// Calculator

virtual ~BenchmarkUtils(){}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


void ComputeHydrodynamicForces(ModelPart& r_model_part)
{
    PropertiesProxiesManager().CreatePropertiesProxies(r_model_part);

    for (ModelPart::ElementIterator i_elem = r_model_part.ElementsBegin(); i_elem != r_model_part.ElementsEnd(); ++i_elem){
        Kratos::SphericParticle& particle = static_cast<Kratos::SphericParticle&>(*i_elem);
        particle.SetFastProperties(mFastProperties);
        particle.Initialize();
        particle.MemberDeclarationFirstStep(r_model_part.GetProcessInfo());
        particle.GetGeometry()[0].Set(INSIDE, true);
        particle.Set(BLOCKED, false);
        array_1d<double, 3> force;
        array_1d<double, 3> moment;
        array_1d<double, 3>& gravity = r_model_part.GetProcessInfo()[GRAVITY];
        particle.ComputeAdditionalForces(force, moment, r_model_part.GetProcessInfo(), gravity);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

std::vector<PropertiesProxy> mFastProperties;
};

}
#endif // KRATOS_BENCHMARK_UTILITIES_H
