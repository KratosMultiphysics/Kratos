//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Miguel Maso Sotomayor
//

#include "utilities/interval_utility.h"

namespace Kratos
{

IntervalUtility::IntervalUtility(Parameters Settings)
{
    KRATOS_TRY

    if(Settings.Has("interval")) {
        if(Settings["interval"][1].IsString() ) {
            if(Settings["interval"][1].GetString() == std::string("End"))
                Settings["interval"][1].SetDouble(1e30);
            else
                KRATOS_ERROR << "the second value of interval can be \"End\" or a number, interval currently: \n"+Settings["interval"].PrettyPrintJsonString();
        }
    } else {
        Parameters defaults(R"( {"default_interval": [0.0, 1e30]} )");
        Settings.AddValue("interval", defaults["default_interval"]);
    }

    mIntervalBegin = Settings["interval"][0].GetDouble();
    mIntervalEnd = Settings["interval"][1].GetDouble();

    KRATOS_CATCH("");
}

double IntervalUtility::GetIntervalBegin() const
{
    return mIntervalBegin;
}

double IntervalUtility::GetIntervalEnd() const
{
    return mIntervalEnd;
}

bool IntervalUtility::IsInInterval(double Time)
{
    const double eps = std::max(1e-14*mIntervalBegin, 1e-30);
    return (Time > mIntervalBegin-eps && Time < mIntervalEnd+eps);
}

std::string IntervalUtility::Info() const
{
    return "IntervalUtility";
}

void IntervalUtility::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void IntervalUtility::PrintData(std::ostream& rOStream) const
{
    rOStream << "[" << GetIntervalBegin() << ", " << GetIntervalEnd() << "]";
}

} // namespace Kratos
