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
#include <iomanip>

// External includes


// Project includes
#include "includes/parallel_environment.h"
#include "input_output/logger_table_output.h"

namespace Kratos
{
  LoggerTableOutput::LoggerTableOutput(std::ostream& rOutputStream, std::vector<std::string> const& ColumnsNames) : LoggerOutput(rOutputStream), mCurrentColumnIndex(0), mColumnsNames(ColumnsNames) {
    for(auto& column_name : mColumnsNames){
      mColumnsWidth.push_back(column_name.size() + 1);
      column_name.erase(column_name.find_last_not_of(" ") + 1);
    }

    const DataCommunicator& r_communicator = ParallelEnvironment::GetDefaultDataCommunicator();
    if (r_communicator.Rank() == 0) {
      WriteHeader();
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
    for(std::size_t i = 0 ; i <  mColumnsNames.size() ; i++)
      std::cout << std::left << std::setw(mColumnsWidth[i]) << mColumnsNames[i];

    std::cout << std::endl;

    for(std::size_t i = 0 ; i <  mColumnsNames.size() ; i++)
      this->GetStream() << std::left << std::setw(mColumnsWidth[i]) << mColumnsNames[i];

    this->GetStream() << std::endl;
    mHeaderIsWritten = true;
  }

  void LoggerTableOutput::WriteMessage(LoggerMessage const& TheMessage)
  {
    int column_index = -1;

    auto message_severity = TheMessage.GetSeverity();
    if (message_severity <= this->GetSeverity() && TheMessage.WriteInThisRank()) {
      if (!mHeaderIsWritten) {
        WriteHeader();
      }
      for (std::size_t i = 0; i < mColumnsNames.size(); i++) {
        if (mColumnsNames[i] == TheMessage.GetLabel()) {
          column_index = static_cast<int>(i);
          break;
        }
      }
    }

    if (column_index >= 0) { // The label found in columns
      MoveCursorToColumn(column_index);
      auto message = TheMessage.GetMessage();
      message.erase(message.find_last_not_of(" \n\t") + 1);
      this->GetStream()  << std::left << std::setw(mColumnsWidth[column_index]) << message;
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

  /// Print object's data.
  void LoggerTableOutput::MoveCursorToColumn(std::size_t ColumnIndex)
  {
    if(ColumnIndex < mCurrentColumnIndex){
      this->GetStream() << std::endl;
      mCurrentColumnIndex = 0;
    }
    for(std::size_t i = mCurrentColumnIndex ; i < ColumnIndex ; i++){
      this->GetStream() << std::setw(mColumnsWidth[i]) << " ";
    }
    mCurrentColumnIndex = ColumnIndex + 1;
  }

}  // namespace Kratos.
