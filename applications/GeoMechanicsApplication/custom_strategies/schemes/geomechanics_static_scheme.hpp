// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//
#pragma once

#include "includes/define.h"
#include "includes/model_part.h"

#include "custom_strategies/schemes/geomechanics_time_integration_scheme.hpp"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeoMechanicsStaticScheme : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsStaticScheme);

    GeoMechanicsStaticScheme()
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>({}, {})
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        // Since this is a static scheme, the derivatives are not updated
    }
};

} // namespace Kratos
