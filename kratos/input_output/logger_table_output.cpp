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
  LoggerTableOutput::LoggerTableOutput(Parameters rSettings): LoggerOutput(mMyFileStream),
    mCurrentColumnIndex(0),
    mColumnsHeaders({}),
    mColumnsWidth({}),
    mColumnsTexts({}),
    mColumnsLabels({}),
    mFileHeader("")
    {
      Parameters default_settings(R"({
          "file_header" : "",
          "file_name"   : "",
          "label"       : "",
          "columns" : [
              {
                  "column_label" : "",
                  "column_header": ""
              }
          ]
      })");

      rSettings.ValidateAndAssignDefaults(default_settings);
      std::string filename = rSettings["file_name"].GetString();
      KRATOS_ERROR_IF(filename=="") << "File not defined for info logger. Please check it is correctly defined" << std::endl;
      mMyFileStream.open(filename, std::ofstream::trunc); // trunc specifies to out only
      this->ApplySettings(rSettings);
    }

  LoggerTableOutput::LoggerTableOutput(std::ostream& rMyStream, Parameters rSettings) : LoggerOutput(rMyStream),
    mCurrentColumnIndex(0),
    mColumnsHeaders({}),
    mColumnsWidth({}),
    mColumnsTexts({}),
    mColumnsLabels({}),
    mFileHeader("")

    {
      Parameters default_settings(R"({
          "file_header" : "",
          "label"       : "",
          "columns" : [
              {
                  "column_label" : "",
                  "column_header": ""
              }
          ]
      })");

      rSettings.ValidateAndAssignDefaults(default_settings);
      this->ApplySettings(rSettings);

    }

 void LoggerTableOutput::ApplySettings(Parameters rSettings)
 {
      mFileHeader = rSettings["file_header"].GetString();
      mIdLabel = rSettings["label"].GetString();
      int ncolumns = rSettings["columns"].size();
      for(int i=0; i<ncolumns; i++)
      {
          std::string coltext = rSettings["columns"][i]["column_header"].GetString();
          mColumnsHeaders.push_back(coltext);
          mColumnsLabels.push_back(rSettings["columns"][i]["column_label"].GetString());
          int tmp = std::max<int>(4,coltext.size());
          mColumnsWidth.push_back(tmp + 1);
          mColumnsTexts.push_back("null");
      }
      const DataCommunicator& r_communicator = ParallelEnvironment::GetDefaultDataCommunicator();
      if(r_communicator.Rank() == 0)
      {
        this->WriteHeader();
      }

 }

  LoggerTableOutput::LoggerTableOutput(LoggerTableOutput const& Other) : LoggerOutput(Other), 
    mCurrentColumnIndex(Other.mCurrentColumnIndex),
    mColumnsHeaders(Other.mColumnsHeaders),
    mColumnsWidth({Other.mColumnsWidth}),
    mColumnsTexts(Other.mColumnsTexts),
    mColumnsLabels(Other.mColumnsLabels),
    mFileHeader(Other.mFileHeader)

  {}

  std::string LoggerTableOutput::Info() const
  {
    return "LoggerTableOutput";
  }

  void LoggerTableOutput::WriteHeader()
  {
    if(!mHeaderIsWritten )
    {
        this->GetStream() << mFileHeader << std::endl << std::endl;
        // WriteHashLine();
        this->GetStream() << " ";

        for(std::size_t i = 0 ; i <  mColumnsHeaders.size() ; i++) {
            this->GetStream() << std::left << Centered(mColumnsWidth[i],mColumnsHeaders[i]) <<  (i < mColumnsHeaders.size()-1 ? " " : "");
        }
        this->GetStream() << std::endl;
        WriteHashLine();
        mHeaderIsWritten = true;
    }
  }

  void LoggerTableOutput::WriteMessage(LoggerMessage const& TheMessage)
  {
    int column_index = -1;
    auto message_severity = TheMessage.GetSeverity();
    std::string label = TheMessage.GetLabel();
    std::string global_label;
    std::string column_label;
    StripLabels(label, global_label, column_label);
    if (global_label == mIdLabel && message_severity <= this->GetSeverity() && TheMessage.WriteInThisRank()) {
      if (!mHeaderIsWritten) {
        WriteHeader();
      }
      for (std::size_t i = 0; i < mColumnsLabels.size(); i++) {
        if (mColumnsLabels[i] == column_label) {
          column_index = static_cast<int>(i);
          break;
        }
      }
    }

    if (column_index >= 0) { // The label found in columns
        auto message = TheMessage.GetMessage();
        message.erase(message.find_last_not_of(" \n\t") + 1);
        mColumnsTexts[column_index] = message;
        WriteTableLine();
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

  void LoggerTableOutput::WriteHashLine()
  {
      this->GetStream() << " ";
      for(std::size_t i = 0 ; i <  mColumnsHeaders.size() ; i++) {
          this->GetStream() << std::string(mColumnsWidth[i]-1, '-') << (i < mColumnsHeaders.size()-1 ? "  " : " ");
      }
      this->GetStream() << std::endl;
  }
  void LoggerTableOutput::WriteTableLine()
    {
      bool has_complete_line = true;
      for(std::size_t i = 0 ; i <  mColumnsTexts.size() ; i++)
      {
          if(mColumnsTexts[i]=="null")
              has_complete_line = false;
      }
      if(has_complete_line)
      {
          for(std::size_t i = 0 ; i <  mColumnsTexts.size() ; i++)
          {
              this->GetStream() << " ";
              this->GetStream() << std::left << Centered(mColumnsWidth[i],mColumnsTexts[i]);
              mColumnsTexts[i]="null";
          }
          this->GetStream() << "" << std::endl;
      }
  }
  std::string LoggerTableOutput::Centered(int width, const std::string& str)
    {
      int len = str.length();
      if(width < len) { return str; }
      return std::string((width-len)/2,' ') + str + std::string(width-len- ((width-len)/2),' ');
  }
  void LoggerTableOutput::StripLabels(std::string& inLabel, std::string& outGlobalLabel, std::string& outColumnLabel)
    {
      outColumnLabel = "";
      std::string delimiter = ".";
      size_t pos = inLabel.find(delimiter);
      outGlobalLabel = (pos != std::string::npos) ? inLabel.substr(0, pos) : ""; // this is to ensure the format is global.column
      outColumnLabel = (pos != std::string::npos) ? inLabel.substr(pos + delimiter.length(), std::string::npos) : "";
  }
}  // namespace Kratos.
