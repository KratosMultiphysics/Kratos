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


#if !defined(KRATOS_LOGGER_OUTPUT_H_INCLUDED )
#define  KRATOS_LOGGER_OUTPUT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <map>


// External includes


// Project includes
#include "includes/define.h"
#include "input_output/logger_message.h"
#include "containers/flags.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// LoggerOutput is the base class for all logger outputs.
/** LoggerOutput defines the interface for all logger outputs
  and also provides the basic (and default) functionalities
  to be extended in other outputs
*/
class KRATOS_API(KRATOS_CORE) LoggerOutput
{
public:
  ///@name Type Definitions
  ///@{

  KRATOS_DEFINE_LOCAL_FLAG(WARNING_PREFIX);
  KRATOS_DEFINE_LOCAL_FLAG(INFO_PREFIX);
  KRATOS_DEFINE_LOCAL_FLAG(DETAIL_PREFIX);
  KRATOS_DEFINE_LOCAL_FLAG(DEBUG_PREFIX);
  KRATOS_DEFINE_LOCAL_FLAG(TRACE_PREFIX);

  /// Pointer definition of LoggerOutput
  KRATOS_CLASS_POINTER_DEFINITION(LoggerOutput);

  ///@}
  ///@name Enums
  ///@{

  ///@}
  ///@name Life Cycle
  ///@{

  explicit LoggerOutput(std::ostream& rOutputStream)
    : mpStream(&rOutputStream),
      mMaxLevel(1),
      mSeverity(LoggerMessage::Severity::INFO),
      mCategory(LoggerMessage::Category::STATUS)
  {
    mOptions.Set(WARNING_PREFIX, true);
    mOptions.Set(INFO_PREFIX, false);
    mOptions.Set(DETAIL_PREFIX, false);
    mOptions.Set(DEBUG_PREFIX, false);
    mOptions.Set(TRACE_PREFIX, false);
  }

  LoggerOutput(LoggerOutput const& Other)
    : mpStream(Other.mpStream), mMaxLevel(Other.mMaxLevel), mSeverity(Other.mSeverity), mCategory(Other.mCategory) {}

  /// Destructor.
  virtual ~LoggerOutput() {}


  ///@}
  ///@name Operators
  ///@{

  LoggerOutput& operator=(LoggerOutput const& Other) = delete;

  ///@}
  ///@name Operations
  ///@{


  ///@}
  ///@name Access
  ///@{

  virtual void WriteHeader();

  virtual void WriteMessage(LoggerMessage const& TheMessage);

  virtual void Flush();

  void SetMaxLevel(std::size_t TheLevel) {
    mMaxLevel = TheLevel;
  }

  std::size_t GetMaxLevel() const {
    return mMaxLevel;
  }

  void SetSeverity(LoggerMessage::Severity const& TheSeverity) {
    mSeverity = TheSeverity;
  }

  LoggerMessage::Severity GetSeverity() const {
    return mSeverity;
  }

  void SetCategory(LoggerMessage::Category const& TheCategory) {
    mCategory = TheCategory;
  }

  LoggerMessage::Category GetCategory() const {
    return mCategory;
  }

  void SetOption(Kratos::Flags ThisFlag, bool Value) {
    mOptions.Set(ThisFlag, Value);
  }

  bool GetOption(Kratos::Flags ThisFlag) {
    return mOptions.Is(ThisFlag);
  }

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const;

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const;

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const;


  /// string stream function
  template<class StreamValueType>
  LoggerOutput& operator << (StreamValueType const& rValue)
  {
    std::stringstream buffer;
    buffer << rValue;

    GetStream() << buffer.str();

    return *this;
  }

  /// Manipulator stream function
  LoggerOutput& operator << (std::ostream& (*pf)(std::ostream&));

  /// char stream function
  LoggerOutput& operator << (const char * rString);



  ///@}
protected:

  LoggerOutput()
    : mpStream(nullptr),
      mMaxLevel(1),
      mSeverity(LoggerMessage::Severity::INFO),
      mCategory(LoggerMessage::Category::STATUS)
  {
    mOptions.Set(WARNING_PREFIX, true);
    mOptions.Set(INFO_PREFIX, false);
    mOptions.Set(DETAIL_PREFIX, false);
    mOptions.Set(DEBUG_PREFIX, false);
    mOptions.Set(TRACE_PREFIX, false);
  }

  std::ostream& GetStream() {return *mpStream;}
  std::ostream* pGetStream() {return mpStream;}
  void SetStream(std::ostream* pStream) {mpStream = pStream;}

  virtual void SetMessageColor(LoggerMessage::Severity MessageSeverity);
  virtual void ResetMessageColor(LoggerMessage::Severity MessageSeverity);

private:
  ///@name Life Cycle
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  std::ostream* mpStream;
  std::size_t mMaxLevel;
  LoggerMessage::Severity mSeverity;
  LoggerMessage::Category mCategory;
  Kratos::Flags mOptions;

  ///@}
}; // Class LoggerOutput

///@}

///@name Input and output
///@{

/// output stream function
std::ostream& operator << (std::ostream& rOStream,
  const LoggerOutput& rThis);

///@}
///@name macros
///@{


///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_LOGGER_OUTPUT_H_INCLUDED  defined
