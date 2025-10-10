// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Gennady Markelov
//

#include "custom_elements/Pw_element.h"

namespace Kratos
{

template class PwElement<2, 2>;
template class PwElement<2, 3>;
template class PwElement<2, 4>;
template class PwElement<2, 5>;
template class PwElement<3, 2>;
template class PwElement<3, 3>;
template class PwElement<3, 4>;

} // Namespace Kratos
