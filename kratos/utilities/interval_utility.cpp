//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "utilities/interval_utility.h"

// System includes
#include <cmath>


namespace Kratos
{
namespace Detail
{

template <class TValue>
IntervalUtility<TValue>::IntervalUtility(Parameters settings)
{
    KRATOS_TRY

    if (!settings.Has("interval")) {
        // No interval provided => set it to ["Begin", "End"]
        settings.AddEmptyArray("interval");
        Parameters interval = settings["interval"];
        interval.Append("Begin");
        interval.Append("End");
    }

    auto interval = settings["interval"];
    KRATOS_ERROR_IF_NOT(settings.IsArray() && settings.size() == 2) << "Expecting 'interval' as an array of size 2, but got:\n" << interval.PrettyPrintJsonString();

    // Replace "Begin" with the minimum representable value
    if(interval[1].Is<std::string>()) {
        if(interval[0].Get<std::string>() == "Begin") {
            interval[0].Set(std::numeric_limits<TValue>::min());
        } else {
            KRATOS_ERROR << "the first value of 'interval' can be \"Begin\" or a number, 'interval' currently:\n" << interval.PrettyPrintJsonString();
        }
    }

    // Replace "End" with the maximum representable value
    /// @todo setting "End" to the maximum representable value will lead to incorrect behaviour
    /// when the test value is also the maximum representable value. This is an unlikely case for
    /// @a double and @a int but still a logical inconsistency.
    if(interval[1].Is<std::string>()) {
        if(interval[1].Get<std::string>() == "End") {
            interval[1].Set(std::numeric_limits<TValue>::max());
        } else {
            KRATOS_ERROR << "the second value of 'interval' can be \"End\" or a number, 'interval' currently:\n" << interval.PrettyPrintJsonString();
        }
    }

    mBegin = interval[0].Get<TValue>();
    mEnd = interval[1].Get<TValue>();

    KRATOS_CATCH("");
}

template <class TValue>
TValue IntervalUtility<TValue>::GetIntervalBegin() const noexcept
{
    return mBegin;
}

template <class TValue>
TValue IntervalUtility<TValue>::GetIntervalEnd() const noexcept
{
    return mEnd;
}

template <>
bool IntervalUtility<double>::IsInInterval(double Value) const noexcept
{
    const double eps = std::max(1e-14 * mBegin, 1e-30);
    return Value > mBegin - eps && Value < mEnd + eps;
}

template <>
bool IntervalUtility<int>::IsInInterval(int Value) const noexcept
{
    return mBegin < Value && Value <= mEnd;
}

template class IntervalUtility<double>;

template class IntervalUtility<int>;

} // namespace Detail
} // namespace Kratos
