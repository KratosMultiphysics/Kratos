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

#include "includes/kratos_export_api.h"
#include "integration_scheme.h"

#include <cstddef>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LobattoIntegrationScheme : public IntegrationScheme
{
public:
    [[nodiscard]] std::size_t                GetNumberOfIntegrationPoints() const;
    [[nodiscard]] IntegrationPointVectorType GetIntegrationPoints() const;
};

} // namespace Kratos
