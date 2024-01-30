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

#include "custom_elements/transient_thermal_element.h"

namespace Kratos
{

template class TransientThermalElement<2, 3>;
template class TransientThermalElement<2, 4>;
template class TransientThermalElement<2, 6>;
template class TransientThermalElement<2, 8>;
template class TransientThermalElement<2, 9>;
template class TransientThermalElement<2, 10>;
template class TransientThermalElement<2, 15>;
template class TransientThermalElement<3, 4>;
template class TransientThermalElement<3, 8>;
template class TransientThermalElement<3, 10>;
template class TransientThermalElement<3, 20>;
template class TransientThermalElement<3, 27>;

} // Namespace Kratos
