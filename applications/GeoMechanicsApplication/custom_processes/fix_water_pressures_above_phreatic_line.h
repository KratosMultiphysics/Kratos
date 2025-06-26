//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
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
    bool          mMoveMesh = false;
};

} // namespace Kratos
