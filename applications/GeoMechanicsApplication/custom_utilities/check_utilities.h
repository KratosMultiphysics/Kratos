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
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "geometries/geometry.h"

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
								
	static void CheckNodalVariables(const Geometry<Node>&            rGeometry,
                                    const std::vector<VariableData>& rVariables);
    static void CheckNodalDof(const Geometry<Node>&            rGeometry,
                              const std::vector<VariableData>& rVariables);
}; /* Class CheckUtilities*/
} /* namespace Kratos.*/
