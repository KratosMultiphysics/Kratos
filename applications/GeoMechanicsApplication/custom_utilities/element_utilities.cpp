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

#include "element_utilities.hpp"

namespace Kratos
{

Element::EquationIdVectorType GeoElementUtilities::ExtractEquationIdsFrom(const Element::DofsVectorType& rDofs)
{
    Element::EquationIdVectorType result;
    std::transform(rDofs.begin(), rDofs.end(), std::back_inserter(result),
                   [](const auto pDof){return pDof->EquationId();});
    return result;
}

} // namespace Kratos
