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

#include "custom_elements/transient_Pw_line_element.h"

namespace Kratos
{

template class TransientPwLineElement<2, 2>;
template class TransientPwLineElement<2, 3>;
template class TransientPwLineElement<2, 4>;
template class TransientPwLineElement<2, 5>;
template class TransientPwLineElement<3, 2>;
template class TransientPwLineElement<3, 3>;
template class TransientPwLineElement<3, 4>;

} // Namespace Kratos
