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


// External includes


// Project includes
#include "includes/define.h"
#include "input_output/logger.h"
#include "utilities/openmp_utils.h"
#include "utilities/timer.h"

namespace Kratos
{

  Logger::Logger(std::string const& TheLabel) : mCurrentMessage(TheLabel)
  {
    mCurrentMessage.SetLevel(GetCurrentLevelInstance());
  }

  Logger::~Logger()
  {
    auto outputs = GetOutputsInstance();
    #pragma omp critical
    {
            GetDefaultOutputInstance().WriteMessage(mCurrentMessage);
      for (auto i_output = outputs.begin(); i_output != outputs.end(); ++i_output)
        (*i_output)->WriteMessage(mCurrentMessage);
    }
  }


  Logger& Logger::Start(std::string const& TheSectionLabel){
    KRATOS_ERROR_IF(OpenMPUtils::IsInParallel() != 0) << "The Logger::Start cannot be called in a parallel region" << std::endl;
    mCurrentMessage.SetLevel(GetCurrentLevelInstance());
    mCurrentMessage << LoggerMessage::START << LoggerMessage::PROFILING;
    GetLabelsStackInstance().push_back(TheSectionLabel);
    auto full_label = CreateFullLabel();
    mCurrentMessage.SetFullLabel(full_label);
    GetCurrentLevelInstance()++;
    Timer::Start(full_label);
    return *this;
  }

  Logger& Logger::Stop(std::string const& TheSectionLabel){
    KRATOS_ERROR_IF(OpenMPUtils::IsInParallel() != 0) << "The Logger::Stop cannot be called in a parallel region" << std::endl;

    if(GetCurrentLevelInstance() > 0){
      GetCurrentLevelInstance()--;
    }

    mCurrentMessage.SetLevel(GetCurrentLevelInstance());
    mCurrentMessage << LoggerMessage::STOP << LoggerMessage::PROFILING;
    Timer::Stop(CreateFullLabel());
    GetLabelsStackInstance().pop_back();
    return *this;
  }

  void Logger::AddOutput(LoggerOutput::Pointer pTheOutput)
  {
    #pragma omp critical
    {
      GetOutputsInstance().push_back(pTheOutput);
    }
  }

  void Logger::Flush() {
    auto outputs = GetOutputsInstance();
    GetDefaultOutputInstance().Flush();
    for (auto i_output = outputs.begin(); i_output != outputs.end(); ++i_output) {
      (*i_output)->Flush();
    }
  }

  std::string Logger::CreateFullLabel(){
    auto& labels_stack = GetLabelsStackInstance();
    std::string result;
    for(auto& label : labels_stack){
      result += "/" + label;
    }
    return result;
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

  /// Flags stream function
  Logger& Logger::operator << (Flags const& TheFlags)
  {
    mCurrentMessage << TheFlags;

    return *this;
  }



}  // namespace Kratos.
