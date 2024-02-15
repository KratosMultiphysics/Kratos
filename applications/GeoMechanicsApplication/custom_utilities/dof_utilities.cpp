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

#include <algorithm>

#include "dof_utilities.h"

namespace Kratos
{

std::vector<std::size_t> ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs)
{
    std::vector<std::size_t> result;
    std::transform(rDofs.begin(), rDofs.end(), std::back_inserter(result),
                   [](const auto p_dof) { return p_dof->EquationId(); });
    return result;
}

} // namespace Kratos
