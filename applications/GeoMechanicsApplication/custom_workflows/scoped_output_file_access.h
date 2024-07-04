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

#include "includes/kratos_export_api.h"

namespace Kratos
{

class StrategyWrapper;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ScopedOutputFileAccess
{
public:
    explicit ScopedOutputFileAccess(StrategyWrapper& rStrategyWrapper);
    ~ScopedOutputFileAccess();
    ScopedOutputFileAccess(const ScopedOutputFileAccess&)            = delete;
    ScopedOutputFileAccess& operator=(const ScopedOutputFileAccess&) = delete;
    ScopedOutputFileAccess(ScopedOutputFileAccess&&)                 = delete;
    ScopedOutputFileAccess& operator=(ScopedOutputFileAccess&&)      = delete;

private:
    StrategyWrapper& mrStrategyWrapper;
};

} // namespace Kratos
