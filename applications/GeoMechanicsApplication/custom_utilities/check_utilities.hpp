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
#include "includes/exception.h"

#include <optional>
#include <string>

namespace Kratos
{

class CheckUtilities
{
public:
    static void CheckDomainSize(double                            DomainSize,
                                std::size_t                       Id,
                                const std::optional<std::string>& PrintName = std::nullopt)
    {
        constexpr auto min_domain_size = 1.0e-15;
        KRATOS_ERROR_IF(DomainSize < min_domain_size)
            << PrintName.value_or("DomainSize") << " (" << DomainSize << ") is smaller than "
            << min_domain_size << " for element " << Id << std::endl;
    }

}; /* Class CheckUtilities*/
} /* namespace Kratos.*/
