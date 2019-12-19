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
        mrStream << " |  /           |             " << std::endl;
        mrStream << " ' /   __| _` | __|  _ \\   __|" << std::endl;
        mrStream << " . \\  |   (   | |   (   |\\__ \\ " << std::endl;
        mrStream << "_|\\_\\_|  \\__,_|\\__|\\___/ ____/" << std::endl;
        mrStream << "           Multi-Physics "<< GetVersionString() << std::endl;
    }

    void LoggerOutput::WriteMessage(LoggerMessage const& TheMessage)
    {
        auto message_severity = TheMessage.GetSeverity();
        if (TheMessage.WriteInThisRank() && message_severity <= mSeverity)
        {
            SetMessageColor(message_severity);

            switch (message_severity)
            {
            case LoggerMessage::Severity::WARNING:
                if (mOptions.Is(WARNING_PREFIX)) mrStream << "[WARNING] ";
                break;
            case LoggerMessage::Severity::INFO:
                if (mOptions.Is(INFO_PREFIX)) mrStream << "[INFO] ";
                break;
            case LoggerMessage::Severity::DETAIL:
                if (mOptions.Is(DETAIL_PREFIX)) mrStream << "[DETAIL] ";
                break;
            case LoggerMessage::Severity::DEBUG:
                if (mOptions.Is(DEBUG_PREFIX)) mrStream << "[DEBUG] ";
                break;
            case LoggerMessage::Severity::TRACE:
                if (mOptions.Is(TRACE_PREFIX)) mrStream << "[TRACE] ";
                break;
            default:
                break;
            }

            if(TheMessage.IsDistributed())
                mrStream << "Rank " << TheMessage.GetSourceRank() << ": ";

            if(TheMessage.GetLabel().size())
                mrStream << TheMessage.GetLabel() << ": " << TheMessage.GetMessage();
            else
                mrStream << TheMessage.GetMessage();

            ResetMessageColor(message_severity);
        }
    }

    void LoggerOutput::Flush()
    {
        mrStream << std::flush;
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
        mrStream << rString;

        return *this;
    }

    LoggerOutput& LoggerOutput::operator << (std::ostream& (*pf)(std::ostream&))
    {
        std::stringstream buffer;
        pf(buffer);

        mrStream << buffer.str();

        return *this;
    }

    void LoggerOutput::SetMessageColor(LoggerMessage::Severity MessageSeverity)
    {
        #if defined(KRATOS_COLORED_OUTPUT)
        if (MessageSeverity == LoggerMessage::Severity::WARNING) mrStream << KYEL;
        #endif
    }

    void LoggerOutput::ResetMessageColor(LoggerMessage::Severity MessageSeverity)
    {
        #if defined(KRATOS_COLORED_OUTPUT)
        mrStream << RST;
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
