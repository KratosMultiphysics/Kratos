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

#include "scoped_output_file_access.h"
#include "strategy_wrapper.hpp"

namespace Kratos
{

ScopedOutputFileAccess::ScopedOutputFileAccess(StrategyWrapper& rStrategyWrapper)
    : mrStrategyWrapper{rStrategyWrapper}
{
    mrStrategyWrapper.InitializeOutput();
}

ScopedOutputFileAccess::~ScopedOutputFileAccess()
{
    mrStrategyWrapper.FinalizeOutput();
}

} // namespace Kratos