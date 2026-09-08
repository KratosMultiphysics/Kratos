//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Carlos Roig
//

// System includes
#include <algorithm>
#include <cstdlib>

// External includes

// Project includes
#include "input_output/logger.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

Logger::Logger(std::string const& TheLabel) : mCurrentMessage(TheLabel)
{
}

Logger::~Logger()
{
    auto outputs = GetOutputsInstance();
    KRATOS_CRITICAL_SECTION
    {
        const bool suppress_message =
            mCurrentMessage.GetCategory() != Logger::Category::CRITICAL
         && std::getenv("KRATOS_QUIET") != nullptr;
        if (!suppress_message) {
            GetDefaultOutputInstance().WriteMessage(mCurrentMessage);
        }

        for (auto it_output = outputs.begin(); it_output != outputs.end(); ++it_output) {
            (*it_output)->WriteMessage(mCurrentMessage);
        }
    }
}

void Logger::AddOutput(LoggerOutput::Pointer pTheOutput)
{
    KRATOS_CRITICAL_SECTION
    {
        GetOutputsInstance().push_back(pTheOutput);
    }
}

void Logger::RemoveOutput(LoggerOutput::Pointer pTheOutput)
{
    KRATOS_TRY

    KRATOS_CRITICAL_SECTION
    {
        auto& r_outputs = GetOutputsInstance();
        auto it_find = std::find(r_outputs.begin(), r_outputs.end(), pTheOutput);
        if (it_find != r_outputs.end()) {
            r_outputs.erase(it_find);
        }
    }

    KRATOS_CATCH("");
}

void Logger::Flush()
{
    auto outputs = GetOutputsInstance();
    GetDefaultOutputInstance().Flush();
    for (auto it_output = outputs.begin(); it_output != outputs.end(); ++it_output) {
        (*it_output)->Flush();
    }
}

std::string Logger::Info() const
{
    return "Logger";
}

/// Print information about this object.
void Logger::PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.
void Logger::PrintData(std::ostream& rOStream) const
{
}

/// Manipulator stream function
Logger& Logger::operator << (std::ostream& (*pf)(std::ostream&))
{
    mCurrentMessage << pf;

    return *this;
}

/// char stream function
Logger& Logger::operator << (const char * rString)
{
    mCurrentMessage << rString;

    return *this;
}

// Location stream function
Logger& Logger::operator << (CodeLocation const& TheLocation)
{
    mCurrentMessage << TheLocation;

    return *this;
}

/// Severity stream function
Logger& Logger::operator << (Severity const& TheSeverity)
{
    mCurrentMessage << TheSeverity;

    return *this;
}

/// Category stream function
Logger& Logger::operator << (Category const& TheCategory)
{
    mCurrentMessage << TheCategory;

    return *this;
}

}  // namespace Kratos.
