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

LobattoIntegrationScheme::LobattoIntegrationScheme(std::size_t NumberOfPoints)
{
    KRATOS_ERROR_IF_NOT(NumberOfPoints == 2) << "Can't construct Lobatto integration scheme: got "
                                             << NumberOfPoints << " point(s)" << std::endl;

    mIntegrationPoints.emplace_back(-1.0, 1.0);
    mIntegrationPoints.emplace_back(1.0, 1.0);
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
