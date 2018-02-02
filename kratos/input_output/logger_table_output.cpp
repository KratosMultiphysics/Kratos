//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <sstream>

// External includes


// Project includes
#include "input_output/logger_table_output.h"


namespace Kratos
{
	LoggerTableOutput::LoggerTableOutput(std::ostream& rOutputStream, std::vector<std::string> const& ColumnsNames) : LoggerOutput(rOutputStream), mColumnsNames(ColumnsNames) {

    }

	LoggerTableOutput::LoggerTableOutput(LoggerTableOutput const& Other) : LoggerOutput(Other), mColumnsNames(Other.mColumnsNames), mColumnsPositions(Other.mColumnsPositions) {

    }

	std::string LoggerTableOutput::Info() const
	{
		return "LoggerTableOutput";
    }
    
    void LoggerTableOutput::WriteHeader()
    {
        for(auto& column_name : mColumnsNames)
            *this << column_name << " ";
        
        *this << std::endl;
    }

    void LoggerTableOutput::WriteMessage(LoggerMessage const& TheMessage)
    {
		auto message_severity = TheMessage.GetSeverity();
        if (message_severity <= this->GetSeverity())
        {
            *this << TheMessage.GetMessage();
        }
    }

	/// Print information about this object.
	void LoggerTableOutput::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void LoggerTableOutput::PrintData(std::ostream& rOStream) const
	{
	}

}  // namespace Kratos.
