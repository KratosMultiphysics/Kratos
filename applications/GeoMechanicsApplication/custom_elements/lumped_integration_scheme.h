// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

#include "includes/kratos_export_api.h"
#include "integration_scheme.h"

namespace Kratos
{
/**
 * @Class LumpedIntegrationScheme
 * @brief Integration scheme used for planar interface elements in GeoMechanicsApplication
 * @details Integration point locations coincide with node-pairs, numbering follows the node-pairs too. Integrations weights as in lumped mass matrix from diagonal scaling of a surface element ( see geometry.h )
 * @author Wijtze Pieter Kikstra
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) LumpedIntegrationScheme : public IntegrationScheme
{
public:
    explicit LumpedIntegrationScheme(std::size_t NumberOfPoints);
    ~LumpedIntegrationScheme() override = default;

    [[nodiscard]] std::size_t GetNumberOfIntegrationPoints() const override;
    [[nodiscard]] const Geo::IntegrationPointVectorType& GetIntegrationPoints() const override;

private:
    static Geo::IntegrationPointVectorType CreateIntegrationPoints(std::size_t NumberOfPoints);

    Geo::IntegrationPointVectorType mIntegrationPoints;
};

} // namespace Kratos
