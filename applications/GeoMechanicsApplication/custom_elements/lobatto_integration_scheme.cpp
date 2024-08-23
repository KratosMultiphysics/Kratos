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
    CreateIntegrationPoints(NumberOfPoints);
}

std::size_t LobattoIntegrationScheme::GetNumberOfIntegrationPoints() const
{
    return mIntegrationPoints.size();
}

const Geo::IntegrationPointVectorType& LobattoIntegrationScheme::GetIntegrationPoints() const
{
    return mIntegrationPoints;
}

void LobattoIntegrationScheme::CreateIntegrationPoints(std::size_t NumberOfPoints)
{
    // A table with the positions of the points and the corresponding weights can be found, for
    // instance, here: https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Lobatto_rules
    switch (NumberOfPoints) {
    case 2:
        mIntegrationPoints.emplace_back(-1.0, 1.0);
        mIntegrationPoints.emplace_back(1.0, 1.0);
        break;

    default:
        KRATOS_ERROR << "Can't construct Lobatto integration scheme: got " << NumberOfPoints
                     << " point(s)" << std::endl;
    }
}

} // namespace Kratos
