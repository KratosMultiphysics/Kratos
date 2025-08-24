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
#include "geo_aliases.h"

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
								
	static void CheckHasNodalSolutionStepData(const Geometry<Node>&             rGeometry, 
        const Geo::ConstVariableDataRefs& rVariableRefs);
    static void CheckHasDofs(const Geometry<Node>& rGeometry, 
        const Geo::ConstVariableDataRefs& rVariableRefs);

private:
    static std::string PrintVectorContent(const std::vector<size_t>& rVector);
}; /* Class CheckUtilities*/
} /* namespace Kratos.*/
