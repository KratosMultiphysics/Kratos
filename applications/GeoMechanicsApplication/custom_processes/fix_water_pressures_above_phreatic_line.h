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

class ModelPart;
class Parameters;

///
/// @brief This process fixes the water pressure DoF to zero for all nodes in the modelpart above the phreatic line.
/// The position of the phreatic line is defined by a table of x and y coordinates (but is constant in time).
/// The current position of the nodes is determined by the initial position plus the total displacement.
/// This means the process requires the TOTAL_DISPLACEMENT variable to be defined and correct at each time step.
class KRATOS_API(GEO_MECHANICS_APPLICATION) FixWaterPressuresAbovePhreaticLineProcess : public Process
{
public:
    FixWaterPressuresAbovePhreaticLineProcess(ModelPart& rMainModelPart, const Parameters& rSettings);
    void                      ExecuteInitializeSolutionStep() override;
    [[nodiscard]] std::string Info() const override;

private:
    ModelPart&    mrModelPart;
    Table<double> mPhreaticLineTable;
};

} // namespace Kratos
