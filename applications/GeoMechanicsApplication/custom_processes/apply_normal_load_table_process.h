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
//

#pragma once

#include "processes/process.h"

#include <memory>
#include <vector>

namespace Kratos
{

class ApplyNormalLoadTableProcess : public Process
{
public:

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

private:
    std::vector<std::unique_ptr<Process>> mProcesses;
};

}
