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


#if !defined(KRATOS_LOGGER_H_INCLUDED )
#define  KRATOS_LOGGER_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "input_output/logger_message.h"
#include "input_output/logger_output.h"
#include "includes/exception.h"



namespace Kratos
{
  ///@addtogroup Kratos
  ///@{

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Logger is in charge of writing the messages to output streams.
  /** Logger is the main class in message writing pipeline which holds an
    array of logger outputs and dispach the arriving logger messages
    to them. Implements a singletone for the list of the outputs and
    also has public constructors and destructors to perform the
    streaming.
  */
  class KRATOS_API(KRATOS_CORE) Logger
    {
    public:
      ///@name Type Definitions
      ///@{

    using LoggerOutputContainerType = std::vector<LoggerOutput::Pointer>;
    ///@}
    ///@name Enums
    ///@{

    using Severity = LoggerMessage::Severity;

    using Category = LoggerMessage::Category;

    using DistributedFilter = LoggerMessage::DistributedFilter;

    ///@}
      ///@name Life Cycle
      ///@{

      explicit Logger(std::string const& TheLabel);

      Logger();

    /// Avoiding Logger to be copied
    Logger(Logger const& rOther) = delete;


      /// Destructor is in charge of passing the message into outputs
      virtual ~Logger();


      ///@}
      ///@name Operators
      ///@{

    /// Loggers can not be assigned.
    Logger& operator=(Logger const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{


      ///@}
      ///@name Static Methods
      ///@{

    static LoggerOutputContainerType& GetOutputsInstance()
    {
      static LoggerOutputContainerType instance;
      return instance;
      }

      static LoggerOutput& GetDefaultOutputInstance()
      {
          static LoggerOutput defaultOutputInstance(std::cout);
          return defaultOutputInstance;
      }

    static void AddOutput(LoggerOutput::Pointer pTheOutput);

    static void Flush();


      ///@}
      ///@name Access
      ///@{

    std::string const& GetCurrentMessage() {
      return mCurrentMessage.GetMessage();
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
    Logger& operator << (StreamValueType const& rValue)
    {
      mCurrentMessage << rValue;

      return *this;
    }

    /// Manipulator stream function
    Logger& operator << (std::ostream& (*pf)(std::ostream&));

    /// char stream function
    Logger& operator << (const char * rString);

    // Location stream function
    Logger& operator << (CodeLocation const& TheLocation);

    /// Severity stream function
    Logger& operator << (Severity const& TheSeverity);

    /// Category stream function
    Logger& operator << (Category const& TheCategory);


      ///@}
     private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

    LoggerMessage mCurrentMessage;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{


      ///@}

    }; // Class Logger

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
            Logger& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
            const Logger& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@name Kratos Macros
  ///@{
// Each-N block
#define KRATOS_LOG_OCCURRENCES_LINE(line) kratos_log_loop_counter##line
#define KRATOS_LOG_OCCURRENCES KRATOS_LOG_OCCURRENCES_LINE(__LINE__)

#define KRATOS_INFO(label) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::INFO
#define KRATOS_INFO_IF(label, conditional) if(conditional) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::INFO
#ifdef KRATOS_DEBUG
  #define KRATOS_INFO_ONCE(label) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES == 0) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::INFO
  #define KRATOS_INFO_FIRST_N(label, logger_count) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES < logger_count) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::INFO
#else
  #define KRATOS_INFO_ONCE(label) if(false) KRATOS_INFO(label)
  #define KRATOS_INFO_FIRST_N(label, logger_count) if(false) KRATOS_INFO(label)
#endif

#define KRATOS_INFO_ALL_RANKS(label) KRATOS_INFO(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_INFO_IF_ALL_RANKS(label, conditional) KRATOS_INFO_IF(label, conditional) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_INFO_ONCE_ALL_RANKS(label) KRATOS_INFO_ONCE(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_INFO_FIRST_N_ALL_RANKS(label, logger_count) KRATOS_INFO_FIRST_N(label, logger_count) << Kratos::Logger::DistributedFilter::FromAllRanks()

#define KRATOS_WARNING(label) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING
#define KRATOS_WARNING_IF(label, conditional) if(conditional) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING
#ifdef KRATOS_DEBUG
  #define KRATOS_WARNING_ONCE(label) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES == 0) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING
  #define KRATOS_WARNING_FIRST_N(label, logger_count) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES < logger_count) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING
#else
  #define KRATOS_WARNING_ONCE(label) if(false) KRATOS_WARNING(label)
  #define KRATOS_WARNING_FIRST_N(label, logger_count) if(false) KRATOS_WARNING(label)
#endif

#define KRATOS_WARNING_ALL_RANKS(label) KRATOS_WARNING(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_WARNING_IF_ALL_RANKS(label, conditional) KRATOS_WARNING_IF(label, conditional) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_WARNING_ONCE_ALL_RANKS(label) KRATOS_WARNING_ONCE(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_WARNING_FIRST_N_ALL_RANKS(label, logger_count) KRATOS_WARNING_FIRST_N(label, logger_count) << Kratos::Logger::DistributedFilter::FromAllRanks()

#define KRATOS_DETAIL(label) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::DETAIL
#define KRATOS_DETAIL_IF(label, conditional) if(conditional) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::DETAIL
#ifdef KRATOS_DEBUG
  #define KRATOS_DETAIL_ONCE(label) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES == 0) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::DETAIL
  #define KRATOS_DETAIL_FIRST_N(label, logger_count) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES < logger_count) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::DETAIL
#else
  #define KRATOS_DETAIL_ONCE(label) if(false) KRATOS_DETAIL(label)
  #define KRATOS_DETAIL_FIRST_N(label, logger_count) if(false) KRATOS_DETAIL(label)
#endif

#define KRATOS_DETAIL_ALL_RANKS(label) KRATOS_DETAIL(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_DETAIL_IF_ALL_RANKS(label, conditional) KRATOS_DETAIL_IF(label, conditional) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_DETAIL_ONCE_ALL_RANKS(label) KRATOS_DETAIL_ONCE(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_DETAIL_FIRST_N_ALL_RANKS(label, logger_count) KRATOS_DETAIL_FIRST_N(label, logger_count) << Kratos::Logger::DistributedFilter::FromAllRanks()

#ifdef KRATOS_DEBUG
#define KRATOS_TRACE(label) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::TRACE
#define KRATOS_TRACE_IF(label, conditional) if(conditional) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::TRACE
#define KRATOS_TRACE_ONCE(label) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES == 0) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::TRACE
#define KRATOS_TRACE_FIRST_N(label, logger_count) static int KRATOS_LOG_OCCURRENCES = -1; if (++KRATOS_LOG_OCCURRENCES < logger_count) Kratos::Logger(label) << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::TRACE
#else
#define KRATOS_TRACE(label) if(false) KRATOS_WARNING(label)
#define KRATOS_TRACE_IF(label, conditional) if(false) KRATOS_WARNING(label)
#define KRATOS_TRACE_ONCE(label) if(false) KRATOS_WARNING(label)
#define KRATOS_TRACE_FIRST_N(label, logger_count) if(false) KRATOS_WARNING(label)
#endif

#define KRATOS_TRACE_ALL_RANKS(label) KRATOS_TRACE(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_TRACE_IF_ALL_RANKS(label, conditional) KRATOS_TRACE_IF(label, conditional) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_TRACE_ONCE_ALL_RANKS(label) KRATOS_TRACE_ONCE(label) << Kratos::Logger::DistributedFilter::FromAllRanks()
#define KRATOS_TRACE_FIRST_N_ALL_RANKS(label, logger_count) KRATOS_TRACE_FIRST_N(label, logger_count) << Kratos::Logger::DistributedFilter::FromAllRanks()

#if defined(KRATOS_ENABLE_CHECK_POINT)
#define KRATOS_CHECK_POINT(label) Kratos::Logger(label) << Kratos::Logger::Category::CHECKING
#else
#define KRATOS_CHECK_POINT(label) \
  if (false)                      \
    Kratos::Logger(label) << Kratos::Logger::Category::CHECKING
#endif
    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LOGGER_H_INCLUDED  defined


