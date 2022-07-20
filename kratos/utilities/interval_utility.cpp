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

    if(settings.Has("interval"))
    {
        if(settings["interval"][1].Is<std::string>())
        {
            if(settings["interval"][1].Get<std::string>() == "End")
                settings["interval"][1].Set(std::numeric_limits<TValue>::max());
            else
                KRATOS_ERROR << "the second value of interval can be \"End\" or a number, interval currently: \n"+settings["interval"].PrettyPrintJsonString();
        }
    }
    else
    {
        // No interval provided => set it to [min, max]
        settings.AddEmptyArray("interval");
        Parameters interval = settings["interval"];
        interval.Append(std::numeric_limits<TValue>::min());
        interval.Append(std::numeric_limits<TValue>::max());
    }

    mBegin = settings["interval"][0].Get<TValue>();
    mEnd = settings["interval"][1].Get<TValue>();

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
