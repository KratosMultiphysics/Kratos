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

    TemporalMethod(ModelPart& rModelPart)
        : mrModelPart(rModelPart), mTotalTime(0.0)
    {
    }

    virtual void InitializeStatisticsVariables()
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base class InitializeStatisticsVariables. "
                        "Please implement it in derrived class.\n";

        KRATOS_CATCH("");
    }

    virtual void InitializeStatisticsMethod()
    {
        mTotalTime = 0.0;
        this->InitializeStatisticsVariables();
    }

    virtual void CalculateStatistics(const double DeltaTime)
    {
        this->FinalizeStatisticsTimeStep(DeltaTime);
    }

    virtual void FinalizeStatisticsTimeStep(const double DeltaTime)
    {
        mTotalTime += DeltaTime;
    }

    ModelPart& GetModelPart() const
    {
        return mrModelPart;
    }

    double GetTotalTime() const
    {
        return mTotalTime;
    }

private:
    ModelPart& mrModelPart;
    double mTotalTime;
};
} // namespace TemporalMethods
} // namespace Kratos
#endif // KRATOS_TEMPORAL_METHOD_H_INCLUDED