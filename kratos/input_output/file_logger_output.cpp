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
    FileLoggerOutput::FileLoggerOutput(const std::string& rFileName)
        : LoggerOutput(),
        mFileStream(rFileName)
    {
        SetStream(&mFileStream);
    }

    std::string FileLoggerOutput::Info() const
    {
        return "FileLoggerOutput";
    }

}  // namespace Kratos.
