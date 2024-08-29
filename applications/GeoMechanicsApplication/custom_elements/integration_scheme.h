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

namespace Kratos
{

class IntegrationScheme
{
public:
    virtual ~IntegrationScheme() = default;

    [[nodiscard]] virtual std::size_t GetNumberOfIntegrationPoints() const                    = 0;
    [[nodiscard]] virtual const Geo::IntegrationPointVectorType& GetIntegrationPoints() const = 0;
};

} // namespace Kratos
