//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_TEMPORAL_METHOD_H_INCLUDED)
#define KRATOS_TEMPORAL_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"

// Application includes

// Application method includes

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace TemporalMethods
{
class TemporalMethod
{
public:
    /// Pointer definition of RansApplyFlagProcess
    KRATOS_CLASS_POINTER_DEFINITION(TemporalMethod);

    TemporalMethod(ModelPart& rModelPart, const int EchoLevel)
        : mrModelPart(rModelPart), mEchoLevel(EchoLevel)
    {
    }

    virtual ~TemporalMethod() = default;

    virtual void InitializeStatisticsVariables()
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base class InitializeStatisticsVariables. "
                        "Please implement it in derrived class.\n";

        KRATOS_CATCH("");
    }

    virtual void InitializeStatisticsMethod(double IntegrationStartTime)
    {
        mIntegrationStartTime = IntegrationStartTime;
        const ProcessInfo& r_process_info = this->GetModelPart().GetProcessInfo();
        if (!r_process_info[IS_RESTARTED])
        {
            this->InitializeStatisticsVariables();
        }
    }

    virtual void CalculateStatistics()
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base class CalculateStatistics. "
                        "Please implement it in derrived class.\n";

        KRATOS_CATCH("");
    }

    ModelPart& GetModelPart() const
    {
        return mrModelPart;
    }

    double GetTotalTime() const
    {
        KRATOS_TRY

        const ProcessInfo& r_process_info = this->GetModelPart().GetProcessInfo();
        const double current_time = r_process_info[TIME];
        const double total_time = current_time - mIntegrationStartTime;

        KRATOS_ERROR_IF(total_time < 0.0)
            << "Total integration time should be greater than or equal to "
               "zero. [ "
               "total_time  = "
            << total_time << ", TIME = " << current_time << " ].\n";

        return total_time;

        KRATOS_CATCH("");
    }

    double GetDeltaTime() const
    {
        const ProcessInfo& r_process_info = this->GetModelPart().GetProcessInfo();
        return r_process_info[DELTA_TIME];
    }

    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

private:
    ModelPart& mrModelPart;
    int mEchoLevel;
    double mIntegrationStartTime;
};
} // namespace TemporalMethods
} // namespace Kratos
#endif // KRATOS_TEMPORAL_METHOD_H_INCLUDED
