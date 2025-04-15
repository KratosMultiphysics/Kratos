// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

#include "processes/process.h"
#include "apply_constant_phreatic_multi_line_pressure_process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) FixModelPartAbovePhreaticLineProcess : public ApplyConstantPhreaticMultiLinePressureProcess
{
    public:
    KRATOS_CLASS_POINTER_DEFINITION(FixModelPartAbovePhreaticLineProcess);

    FixModelPartAbovePhreaticLineProcess(ModelPart& rModelPart, const Parameters& rSettings);

    ~FixModelPartAbovePhreaticLineProcess() override;

    FixModelPartAbovePhreaticLineProcess(const FixModelPartAbovePhreaticLineProcess&)            = delete;
    FixModelPartAbovePhreaticLineProcess& operator=(const FixModelPartAbovePhreaticLineProcess&) = delete;

    void ExecuteInitializeSolutionStep() override;
    void ExecuteInitialize() override;
    double CalculateDistanceToPhreaticLine(const Node& rNode) const;
    int FindIndex(const Node& rNode) const;

    private:
    bool mMoveMeshFlag;
};

} // namespace Kratos
