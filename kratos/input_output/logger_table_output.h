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


#if !defined(KRATOS_LOGGER_TABLE_OUTPUT_H_INCLUDED )
#define  KRATOS_LOGGER_TABLE_OUTPUT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>


// External includes


// Project includes
#include "input_output/logger_output.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// LoggerTableOutput takes columns names and only prints the messages with lable given with column name.
/** The columns width will be at least equal to the lable size and can be extended by adding additional
 *  spaces to the end of the column name: "Time Step" width would be 9 and "Time Step   " width would be 12
 */
class KRATOS_API(KRATOS_CORE) LoggerTableOutput : public LoggerOutput
{
public:
  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Enums
  ///@{

  ///@}
  ///@name Life Cycle
  ///@{

  LoggerTableOutput(std::ostream& rOutputStream, std::vector<std::string> const& ColumnsNames);

  LoggerTableOutput(LoggerTableOutput const& Other);

  /// Destructor.
  virtual ~LoggerTableOutput() {}


  ///@}
  ///@name Operators
  ///@{

  LoggerTableOutput& operator=(LoggerTableOutput const& Other) = delete;

  ///@}
  ///@name Operations
  ///@{


  ///@}
  ///@name Access
  ///@{

  void WriteHeader() override;

  void WriteMessage(LoggerMessage const& TheMessage) override;

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const override;

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override;

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override;


  ///@}

private:
  ///@name Life Cycle
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  std::size_t mCurrentColumnIndex;
  std::vector<std::string> mColumnsNames;
  std::vector<std::size_t> mColumnsWidth;

  bool mHeaderIsWritten = false;

  ///@}
  ///@name Private Operations
  ///@{

  void MoveCursorToColumn(std::size_t ColumnIndex);

  ///@}
}; // Class LoggerTableOutput

///@}

///@name Input and output
///@{

/// output stream function
std::ostream& operator << (std::ostream& rOStream,
  const LoggerTableOutput& rThis);

///@}
///@name macros
///@{


///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_LOGGER_TABLE_OUTPUT_H_INCLUDED  defined
