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

#include "lobatto_integration_scheme.h"

namespace Kratos
{

LobattoIntegrationScheme::LobattoIntegrationScheme(Geo::IntegrationPointVectorType IntegrationPoints)
    : mIntegrationPoints(std::move(IntegrationPoints))
{
}

std::size_t LobattoIntegrationScheme::GetNumberOfIntegrationPoints() const
{
    return mIntegrationPoints.size();
}

const Geo::IntegrationPointVectorType& LobattoIntegrationScheme::GetIntegrationPoints() const
{
    return mIntegrationPoints;
}

} // namespace Kratos
