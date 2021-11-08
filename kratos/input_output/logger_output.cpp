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
//
//

// System includes
#include <sstream>

// External includes


// Project includes
#include "input_output/logger_output.h"
#include "includes/kratos_version.h"

#if defined(KRATOS_COLORED_OUTPUT)
#include "utilities/color_utilities.h"
#endif

namespace Kratos
{

    KRATOS_CREATE_LOCAL_FLAG( LoggerOutput, WARNING_PREFIX,  0 );
    KRATOS_CREATE_LOCAL_FLAG( LoggerOutput, INFO_PREFIX,     1 );
    KRATOS_CREATE_LOCAL_FLAG( LoggerOutput, DETAIL_PREFIX,   2 );
    KRATOS_CREATE_LOCAL_FLAG( LoggerOutput, DEBUG_PREFIX,    3 );
    KRATOS_CREATE_LOCAL_FLAG( LoggerOutput, TRACE_PREFIX,    4 );

    std::string LoggerOutput::Info() const
    {
        return "LoggerOutput";
    }

    void LoggerOutput::WriteHeader()
    {
        auto& r_stream = GetStream();
        r_stream << " |  /           |             " << std::endl;
        r_stream << " ' /   __| _` | __|  _ \\   __|" << std::endl;
        r_stream << " . \\  |   (   | |   (   |\\__ \\ " << std::endl;
        r_stream << "_|\\_\\_|  \\__,_|\\__|\\___/ ____/" << std::endl;
        r_stream << "           Multi-Physics "<< GetVersionString() << std::endl;
    }

    void LoggerOutput::WriteMessage(LoggerMessage const& TheMessage)
    {
        auto& r_stream = GetStream();
        auto message_severity = TheMessage.GetSeverity();
        if (TheMessage.WriteInThisRank() && message_severity <= mSeverity)
        {
            SetMessageColor(message_severity);

            switch (message_severity)
            {
            case LoggerMessage::Severity::WARNING:
                if (mOptions.Is(WARNING_PREFIX)) r_stream << "[WARNING] ";
                break;
            case LoggerMessage::Severity::INFO:
                if (mOptions.Is(INFO_PREFIX)) r_stream << "[INFO] ";
                break;
            case LoggerMessage::Severity::DETAIL:
                if (mOptions.Is(DETAIL_PREFIX)) r_stream << "[DETAIL] ";
                break;
            case LoggerMessage::Severity::DEBUG:
                if (mOptions.Is(DEBUG_PREFIX)) r_stream << "[DEBUG] ";
                break;
            case LoggerMessage::Severity::TRACE:
                if (mOptions.Is(TRACE_PREFIX)) r_stream << "[TRACE] ";
                break;
            default:
                break;
            }

            if(TheMessage.IsDistributed())
                r_stream << "Rank " << TheMessage.GetSourceRank() << ": ";

            if(TheMessage.GetLabel().size())
                r_stream << TheMessage.GetLabel() << ": " << TheMessage.GetMessage();
            else
                r_stream << TheMessage.GetMessage();

            ResetMessageColor(message_severity);
        }
    }

    void LoggerOutput::Flush()
    {
        GetStream() << std::flush;
    }

    /// Print information about this object.
    void LoggerOutput::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void LoggerOutput::PrintData(std::ostream& rOStream) const
    {
        rOStream << "Max Level : " << mMaxLevel;
    }

    /// char stream function
    LoggerOutput& LoggerOutput::operator << (const char * rString)
    {
        GetStream() << rString;

        return *this;
    }

    LoggerOutput& LoggerOutput::operator << (std::ostream& (*pf)(std::ostream&))
    {
        std::stringstream buffer;
        pf(buffer);

        GetStream() << buffer.str();

        return *this;
    }

    void LoggerOutput::SetMessageColor(LoggerMessage::Severity MessageSeverity)
    {
        #if defined(KRATOS_COLORED_OUTPUT)
        if (MessageSeverity == LoggerMessage::Severity::WARNING) GetStream() << KYEL;
        #endif
    }

    void LoggerOutput::ResetMessageColor(LoggerMessage::Severity MessageSeverity)
    {
        #if defined(KRATOS_COLORED_OUTPUT)
        GetStream() << RST;
        #endif
    }

    /// output stream function
    std::ostream& operator << (std::ostream& rOStream,
        const LoggerOutput& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }


}  // namespace Kratos.
