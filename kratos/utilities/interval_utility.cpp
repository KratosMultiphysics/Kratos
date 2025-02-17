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

// Project includes
#include "utilities/interval_utility.h"

// System includes
#include <limits> // std::numeric_limits


namespace Kratos {

IntervalUtility::IntervalUtility(Parameters Settings)
{
    KRATOS_TRY
    Settings.ValidateAndAssignDefaults(this->GetSchema());
    if (Settings["interval"][0].Is<std::string>()) Settings["interval"][0].Set<double>(std::numeric_limits<double>::lowest());
    if (Settings["interval"][1].Is<std::string>()) Settings["interval"][1].Set<double>(std::numeric_limits<double>::max());

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

Schema IntervalUtility::GetSchema() const
{
    return Parameters(R"({
"title" : "IntervalUtility",
"description" : "A class for performing inclusive point membership tests in 1D.",
"type" : "object",
"properties" : {
    "interval" : {
        "type" : "array",
        "items" : [
        {
            "type" : "number",
            "pattern" : "Begin"
        },
        {
            "type" : ["number", "string"],
            "pattern" : "End"
        }
        ],
        "additionalItems" : false,
        "default" : ["Begin", "End"]
    }
},
"required" : []
})");
}

std::ostream& operator << (std::ostream& rOStream,
    const IntervalUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos
