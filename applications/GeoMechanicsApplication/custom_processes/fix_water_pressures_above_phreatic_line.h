// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Jonathan Nuttall
//

#pragma once

#include "includes/table.h"
#include "processes/process.h"

namespace Kratos
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) FixWaterPressuresAbovePhreaticLineProcess : public Process
{
public:
    FixWaterPressuresAbovePhreaticLineProcess(ModelPart& rMainModelPart, const Parameters& rSettings);
    void ExecuteInitializeSolutionStep() override;

private:
    ModelPart&    mrModelPart;
    Table<double> mPhreaticLineTable;
};

} // namespace Kratos
