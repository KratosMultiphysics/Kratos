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
//                   Riccardo Rossi
//                   Miguel Maso Sotomayor
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
IntervalUtility<TValue>::IntervalUtility()
    : IntervalUtility(IntervalUtility<TValue>::GetDefaultParameters())
{
}


template <class TValue>
IntervalUtility<TValue>::IntervalUtility(Parameters settings)
{
    KRATOS_TRY

    if (!settings.Has("interval")) {
        settings.AddStringArray("interval", {"Begin", "End"});
    }

    auto interval = settings["interval"];
    KRATOS_ERROR_IF_NOT(interval.IsArray() && interval.size() == 2) << "Expecting \"interval\" as an array of size 2, but got:\n" << interval.PrettyPrintJsonString();

    // Replace "Begin" with the minimum representable value
    if(interval[0].Is<std::string>()) {
        if(interval[0].Get<std::string>() == "Begin") {
            interval[0].Set(std::numeric_limits<TValue>::lowest());
        } else {
            KRATOS_ERROR << "the first value of \"interval\" can be \"Begin\" or a number, \"interval\" currently:\n" << interval.PrettyPrintJsonString();
        }
    } else {
        KRATOS_ERROR_IF_NOT(interval[0].Is<TValue>()) << "the first value of \"interval\" can be \"Begin\" or a number, \"interval\" currently:\n" << interval.PrettyPrintJsonString();
    }

    // Replace "End" with the maximum representable value
    /// @todo setting "End" to the maximum representable value will lead to incorrect behaviour
    /// when the test value is also the maximum representable value. This is an unlikely case for
    /// @a double and @a int but still a logical inconsistency.
    if(interval[1].Is<std::string>()) {
        if(interval[1].Get<std::string>() == "End") {
            interval[1].Set(std::numeric_limits<TValue>::max());
        } else {
            KRATOS_ERROR << "the second value of \"interval\" can be \"End\" or a number, \"interval\" currently:\n" << interval.PrettyPrintJsonString();
        }
    } else {
        KRATOS_ERROR_IF_NOT(interval[1].Is<TValue>()) << "the second value of \"interval\" can be \"End\" or a number, \"interval\" currently:\n" << interval.PrettyPrintJsonString();
    }

    this->SetBoundaries(interval[0].Get<TValue>(), interval[1].Get<TValue>());
    KRATOS_ERROR_IF(mEnd < mBegin) << "Invalid \"interval\":\n" << interval.PrettyPrintJsonString();

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

template <class TValue>
Parameters IntervalUtility<TValue>::GetDefaultParameters()
{
    return Parameters(R"({"interval" : ["Begin", "End"]})");
}


template <>
void IntervalUtility<double>::SetBoundaries(double begin, double end) noexcept
{
    const double relative_tolerance = 1e-14;
    const double absolute_tolerance = 1e-30;

    // Apply tolerances and avoid under/overflows
    mBegin = std::min(begin, begin - std::max(relative_tolerance * begin, absolute_tolerance));
    mEnd = std::max(end, end + std::max(relative_tolerance * end, absolute_tolerance));
}


template <>
void IntervalUtility<int>::SetBoundaries(int begin, int end) noexcept
{
    mBegin = begin;
    mEnd = end;
}


template <>
std::string IntervalUtility<double>::Info() const
{
    return "IntervalUtility";
}


template <>
std::string IntervalUtility<int>::Info() const
{
    return "DiscreteIntervalUtility";
}

template <class TValue>
void IntervalUtility<TValue>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}


template <class TValue>
void IntervalUtility<TValue>::PrintData(std::ostream& rOStream) const
{
    rOStream << "[" << GetIntervalBegin() << ", " << GetIntervalEnd() << "]";
}


template <class TValue>
inline std::ostream& operator << (std::ostream& rOStream, const IntervalUtility<TValue>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);
    return rOStream;
}


template class IntervalUtility<double>;

template class IntervalUtility<int>;

} // namespace Detail
} // namespace Kratos
