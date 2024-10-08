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

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LobattoIntegrationScheme : public IntegrationScheme
{
public:
    explicit LobattoIntegrationScheme(std::size_t NumberOfPoints);
    ~LobattoIntegrationScheme() override = default;

    [[nodiscard]] std::size_t GetNumberOfIntegrationPoints() const override;
    [[nodiscard]] const Geo::IntegrationPointVectorType& GetIntegrationPoints() const override;

private:
    static Geo::IntegrationPointVectorType CreateIntegrationPoints(std::size_t NumberOfPoints);

    Geo::IntegrationPointVectorType mIntegrationPoints;
};

} // namespace Kratos
