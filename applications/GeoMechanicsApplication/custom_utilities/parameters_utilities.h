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
//                   Marjan Fathian,
//                   Gennady Markelov
//

#pragma once

#include "includes/kratos_parameters.h"

namespace Kratos
{

class ParametersUtilities
{
public:
    static Parameters CopyRequiredParameters(const Parameters& rSourceParameters,
                                             const std::vector<std::string>& rNamesOfParametersToCopy);

    static Parameters CopyOptionalParameters(const Parameters& rSourceParameters,
                                             const std::vector<std::string>& rNamesOfParametersToCopy);

    static void AppendParameterNameIfExists(const std::string& rParameterName,
                                            const Parameters& rSourceParameters,
                                            std::vector<std::string>& rResult);

    static bool HasTableAttached(const Parameters& rSettings);
    static bool HasTableAttached(const Parameters& rSettings, int component);
};

} // namespace Kratos
