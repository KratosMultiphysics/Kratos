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
#include <iomanip>

// External includes


// Project includes
#include "input_output/logger_table_output.h"


namespace Kratos
{
	LoggerTableOutput::LoggerTableOutput(std::ostream& rOutputStream, std::vector<std::string> const& ColumnsNames) : LoggerOutput(rOutputStream), mColumnsNames(ColumnsNames) {
        for(auto& column_name : mColumnsNames){
            mColumnsWidth.push_back(column_name.size());
            column_name.erase(column_name.find_last_not_of(" ") + 1);
        }
    }

	LoggerTableOutput::LoggerTableOutput(LoggerTableOutput const& Other) : LoggerOutput(Other), mColumnsNames(Other.mColumnsNames), mColumnsWidth(Other.mColumnsWidth) {

    }

	std::string LoggerTableOutput::Info() const
	{
		return "LoggerTableOutput";
    }
    
    void LoggerTableOutput::WriteHeader()
    {
        for(int i = 0 ; i <  mColumnsNames.size() ; i++)
            std::cout << std::left << std::setw(mColumnsWidth[i]+1) << mColumnsNames[i];
        
        std::cout << std::endl;

        for(int i = 0 ; i <  mColumnsNames.size() ; i++)
            this->GetStream() << std::left << std::setw(mColumnsWidth[i]+1) << mColumnsNames[i];
        
        this->GetStream() << std::endl;
    }

    void LoggerTableOutput::WriteMessage(LoggerMessage const& TheMessage)
    {
      int column_index = -1;

      auto message_severity = TheMessage.GetSeverity();
      if (message_severity <= this->GetSeverity()) {
        for (int i = 0; i < mColumnsNames.size(); i++) {
          if (mColumnsNames[i] == TheMessage.GetLabel()) {
            column_index = i;
            break;
          }
        }
      }
      KRATOS_WATCH(column_index)
      if (column_index >= 0) {
        this->GetStream() << TheMessage.GetMessage();
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
