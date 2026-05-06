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

/// <summary>
/// Process to calculate the total motion of nodes in a model part. Where total motion is either
/// TOTAL_DISPLACEMENT or TOTAL_ROTATION.
/// </summary>
class KRATOS_API(GEO_MECHANICS_APPLICATION) CalculateTotalMotionProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CalculateTotalMotionProcess);

    CalculateTotalMotionProcess(ModelPart& rModelPart, const Parameters& rSettings);

    CalculateTotalMotionProcess(const CalculateTotalMotionProcess&)            = delete;
    CalculateTotalMotionProcess& operator=(const CalculateTotalMotionProcess&) = delete;

    void                      Execute() override;
    [[nodiscard]] std::string Info() const override;

private:
    ModelPart&  mrModelPart;
    std::string mResultsVariableName;
    std::string mIncrementalVariableName;
};

} // namespace Kratos
