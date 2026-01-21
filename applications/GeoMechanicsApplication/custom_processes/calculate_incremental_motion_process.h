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

#include <string>

namespace Kratos
{
class ModelPart;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) CalculateIncrementalMotionProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CalculateIncrementalMotionProcess);

    CalculateIncrementalMotionProcess(ModelPart& rModelPart, const Parameters& rSettings);

    CalculateIncrementalMotionProcess(const CalculateIncrementalMotionProcess&)            = delete;
    CalculateIncrementalMotionProcess& operator=(const CalculateIncrementalMotionProcess&) = delete;

    void                      Execute() override;
    [[nodiscard]] std::string Info() const override;

private:
    ModelPart&  mrModelPart;
    std::string mResultsVariableName;
    std::string mBaseVariableName;
};

} // namespace Kratos
