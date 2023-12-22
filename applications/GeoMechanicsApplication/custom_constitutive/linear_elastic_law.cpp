// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "linear_elastic_law.h"

namespace Kratos
{

bool GeoLinearElasticLaw::RequiresInitializeMaterialResponse()
{
    // Avoid that `ConstitutiveLaw::InitializeMaterialResponseCauchy` throws
    // an exception
    return false;
}

ConstitutiveLaw::StrainMeasure GeoLinearElasticLaw::GetStrainMeasure()
{
    return StrainMeasure_Infinitesimal;
}

ConstitutiveLaw::StressMeasure GeoLinearElasticLaw::GetStressMeasure()
{
    return StressMeasure_Cauchy;
}

void GeoLinearElasticLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
}

void GeoLinearElasticLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
}

} // namespace Kratos
