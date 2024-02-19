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

#pragma once

#include <vector>

#include "includes/dof.h"

namespace Kratos
{

std::vector<std::size_t> ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs);

Vector ExtractSolutionStepValues(const std::vector<Dof<double>*>& rDofs, int Step);
Vector ExtractFirstDerivatives(const std::vector<Dof<double>*>& rDofs, int Step);

Vector ExtractSolutionStepValuesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int Step);
Vector ExtractFirstTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int Step);
Vector ExtractSecondTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int Step);

} // namespace Kratos
