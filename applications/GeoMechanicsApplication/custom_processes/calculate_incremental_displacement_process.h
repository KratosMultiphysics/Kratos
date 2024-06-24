// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//
#pragma once

#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) CalculateIncrementalDisplacementProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CalculateIncrementalDisplacementProcess);

    CalculateIncrementalDisplacementProcess(ModelPart& rModelPart, const Parameters& rSettings);

    ~CalculateIncrementalDisplacementProcess() override;

    CalculateIncrementalDisplacementProcess(const CalculateIncrementalDisplacementProcess&)            = delete;
    CalculateIncrementalDisplacementProcess& operator=(const CalculateIncrementalDisplacementProcess&) = delete;

    void Execute() override;

private:
    ModelPart& mrModelPart;
};

} // namespace Kratos