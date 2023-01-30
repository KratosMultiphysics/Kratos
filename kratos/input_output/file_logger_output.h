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


#if !defined(KRATOS_FILE_LOGGER_OUTPUT_H_INCLUDED )
#define  KRATOS_FILE_LOGGER_OUTPUT_H_INCLUDED



// System includes
#include <fstream>

// External includes

// Project includes
#include "input_output/logger_output.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// FileLoggerOutput is a class for all file logger outputs.
/** FileLoggerOutput creates a log file and writes to it.
*/
class KRATOS_API(KRATOS_CORE) FileLoggerOutput : public LoggerOutput
{
public:
  ///@name Type Definitions
  ///@{


  /// Pointer definition of FileLoggerOutput
  KRATOS_CLASS_POINTER_DEFINITION(FileLoggerOutput);

  ///@}
  ///@name Enums
  ///@{

  ///@}
  ///@name Life Cycle
  ///@{

  explicit FileLoggerOutput(const std::string& rName);

  /// Destructor.
  ~FileLoggerOutput() {};

  ///@}
  ///@name Operators
  ///@{

  FileLoggerOutput& operator=(FileLoggerOutput const& Other) = delete;

  ///@}
  ///@name Operations
  ///@{

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const override;

  ///@}
protected:

  void SetMessageColor(LoggerMessage::Severity MessageSeverity) override {};
  void ResetMessageColor(LoggerMessage::Severity MessageSeverity) override {};

private:
  ///@name Life Cycle
  ///@{

  ///@}
  ///@name Member Variables
  ///@{
  std::ofstream mFileStream;

  ///@}
}; // Class FileLoggerOutput

///@}

///@name Input and output
///@{

///@}
///@name macros
///@{


///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_FILE_LOGGER_OUTPUT_H_INCLUDED  defined
