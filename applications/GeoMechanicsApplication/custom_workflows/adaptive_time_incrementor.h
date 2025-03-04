// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "includes/kratos_export_api.h"
#include "time_incrementor.h"
#include <string>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) AdaptiveTimeIncrementor : public TimeIncrementor
{
public:
    AdaptiveTimeIncrementor(double                         StartTime,
                            double                         EndTime,
                            double                         StartIncrement,
                            std::pair<std::string, double> MinAllowableDeltaTime,
                            std::size_t                    MaxNumOfCycles     = 10,
                            double                         ReductionFactor    = 0.5,
                            double                         IncreaseFactor     = 2.0,
                            double                         MaxTimeStepFactor  = 1000.0,
                            std::size_t                    MinNumOfIterations = 3,
                            std::size_t                    MaxNumOfIterations = 15);

    [[nodiscard]] bool WantNextStep(const TimeStepEndState& rPreviousState) const override;
    [[nodiscard]] bool WantRetryStep(std::size_t CycleNumber, const TimeStepEndState& rPreviousState) const override;
    [[nodiscard]] double GetIncrement() const override;
    void                 PostTimeStepExecution(const TimeStepEndState& rResultantState) override;

private:
    double                         mEndTime;
    double                         mDeltaTime;
    std::pair<std::string, double> mMinAllowableDeltaTime;
    std::size_t                    mMaxNumOfCycles;
    double                         mReductionFactor;
    double                         mIncreaseFactor;
    double                         mMaxDeltaTime;
    std::size_t                    mMinNumOfIterations;
    std::size_t                    mMaxNumOfIterations;
};

} // namespace Kratos
