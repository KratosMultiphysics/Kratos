//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <sstream>

// Project includes

// Include base h
#include "fixed_interval_recurring_event_utility.h"

namespace Kratos {

template <class TDataType>
FixedIntervalRecurringEventUtility<TDataType>::FixedIntervalRecurringEventUtility(
    const TDataType InitialValue,
    const TDataType Interval)
    : mInterval(Interval),
      mNextEventValue(TDataType{})
{
    ScheduleNextEvent(InitialValue);
}

template <class TDataType>
bool FixedIntervalRecurringEventUtility<TDataType>::IsEventExpected(const TDataType CurrentValue) const
{
    return CurrentValue >= mNextEventValue;
}

template <class TDataType>
void FixedIntervalRecurringEventUtility<TDataType>::ScheduleNextEvent(const TDataType CurrentValue)
{
    if (mInterval > TDataType{}) {
        while (mNextEventValue <= CurrentValue) {
            mNextEventValue += mInterval;
        }
    }
}

template <class TDataType>
std::string FixedIntervalRecurringEventUtility<TDataType>::Info() const
{
    std::stringstream msg;

    msg << "FixedIntervalRecurringEventUtility [ Interval = "
        << mInterval << ", Next event expected at = "
        << mNextEventValue << " ]\n";

    return msg.str();
}

// template instantiations
template class FixedIntervalRecurringEventUtility<int>;
template class FixedIntervalRecurringEventUtility<double>;

} // namespace Kratos
