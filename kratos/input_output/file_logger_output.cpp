//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Armin Geiser
//
//

// System includes


// External includes


// Project includes
#include "input_output/file_logger_output.h"


namespace Kratos
{
    namespace {
    std::ofstream* CreateNewLogFile(const std::string& rFileName)
    {
        // The file stream created here on the heap is deleted in the destructor
        // of the FileLoggerOutput
        return new std::ofstream(rFileName);
    }
    } // namespace

    FileLoggerOutput::FileLoggerOutput(const std::string& rFileName)
        : LoggerOutput(*CreateNewLogFile(rFileName))
    { }

    FileLoggerOutput::~FileLoggerOutput() {
        delete &GetStream();
    }

    void FileLoggerOutput::WriteMessage(LoggerMessage const& TheMessage)
    {
        auto message_severity = TheMessage.GetSeverity();
        if (TheMessage.WriteInThisRank() && message_severity <= GetSeverity())
        {
            WriteSeverityPrefix(message_severity);

            if(TheMessage.IsDistributed())
                GetStream() << "Rank " << TheMessage.GetSourceRank() << ": ";

            if(TheMessage.GetLabel().size())
                GetStream() << TheMessage.GetLabel() << ": " << TheMessage.GetMessage();
            else
                GetStream() << TheMessage.GetMessage();
        }
    }

    std::string FileLoggerOutput::Info() const
    {
        return "FileLoggerOutput";
    }

}  // namespace Kratos.
