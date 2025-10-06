// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.hpp"

namespace Kratos
{
template class UpdatedLagrangianUPwDiffOrderElement<2, 6>;
template class UpdatedLagrangianUPwDiffOrderElement<2, 8>;
template class UpdatedLagrangianUPwDiffOrderElement<2, 9>;
template class UpdatedLagrangianUPwDiffOrderElement<2, 10>;
template class UpdatedLagrangianUPwDiffOrderElement<2, 15>;
template class UpdatedLagrangianUPwDiffOrderElement<3, 10>;
template class UpdatedLagrangianUPwDiffOrderElement<3, 20>;
template class UpdatedLagrangianUPwDiffOrderElement<3, 27>;
} // Namespace Kratos
