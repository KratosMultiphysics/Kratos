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

#pragma once

#include "includes/constitutive_law.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoLinearElasticLaw : public ConstitutiveLaw
{
public:
    bool RequiresInitializeMaterialResponse() override;

    StrainMeasure GetStrainMeasure() override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
