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

//#include <memory>
//#include <vector>

namespace Kratos
{

class ModelPart;
class Parameters;


class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyNormalLoadTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyNormalLoadTableProcess);

    ApplyNormalLoadTableProcess(ModelPart&        rModelPart,
                                const Parameters& rProcessSettings);

    ~ApplyNormalLoadTableProcess() override;

    ApplyNormalLoadTableProcess(const ApplyNormalLoadTableProcess&) = delete;
    ApplyNormalLoadTableProcess& operator=(const ApplyNormalLoadTableProcess&) = delete;

    void ExecuteInitialize() override;
    void ExecuteInitializeSolutionStep() override;

private:
    ModelPart& mrModelPart;
    std::vector<std::unique_ptr<Process>> mProcesses;
};

}
