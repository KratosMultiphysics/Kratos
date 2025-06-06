// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes
#include "includes/element.h"

#include <optional>
#include <string>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) CheckUtilities
{
public:
    static void CheckDomainSize(double                            DomainSize,
                                std::size_t                       Id,
                                const std::optional<std::string>& PrintName = std::nullopt);

    static void CheckHasSolutionStepsDataFor(const Geometry<Node>& rGeometry, const VariableData& rVariable);
    static void CheckUtilities::CheckHasDofsFor(const Geometry<Node>& rGeometry, const Variable<double>& rVariable);
    static void CheckUtilities::CheckProperty(size_t                          Id,
                                              const Properties&               rProperties,
                                              const Kratos::Variable<double>& rVariable,
                                              std::optional<double> MaxValue = std::nullopt);
    static void CheckUtilities::CheckProperty(size_t                               Id,
                                              const Properties&                    rProperties,
                                              const Kratos::Variable<std::string>& rVariable,
                                              const std::string&                   rName);
    static void CheckForNonZeroZCoordinateIn2D(size_t Dimension, const Geometry<Node>& rGeometry);

}; /* Class CheckUtilities*/
} /* namespace Kratos.*/
