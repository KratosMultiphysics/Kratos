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

#include "stress_state_policy.h"

namespace Kratos
{

class StressStatePolicyFactory
{
public:
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Create() const;
};

} // namespace Kratos
