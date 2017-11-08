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


#if !defined(KRATOS_KRATOS_LOG_H_INCLUDED )
#define  KRATOS_KRATOS_LOG_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

// External includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "utilities/stl_io.h"

// Filters
#define KRATOS_LOG_ADD_CUSTOM_FILTER(STRING) KratosLogUtils::GetInstance().AddCustomFilter(STRING);
#define KRATOS_LOG_FILTER(STRING)            KratosLogUtils::GetInstance().Filter(STRING)  

// Stamps
#define KRATOS_TIME_STAMP "[" << KratosLogUtils::GetInstance().CurrentDateTime() << "]"
#define KRATOS_PROCESS_ID "[PID=" << getpid() << "]"

#define KRATOS_LOG_ERROR_STAMP_DETAIL(file,line,function) "ERROR:" << (file) << ":" << (line) << ":" << "\n" << "\t in function: " << function << ":\n"
#define KRATOS_LOG_ERROR_STAMP KRATOS_LOG_ERROR_STAMP_DETAIL(__FILE__,__LINE__,KRATOS_LOG_FILTER(BOOST_CURRENT_FUNCTION))

#define KRATOS_LOG_WARNING_STAMP_DETAIL(file,line) "WARNING:" << (file) << ":" << (line) << ":"
#define KRATOS_LOG_WARNING_STAMP KRATOS_LOG_WARNING_STAMP_DETAIL(__FILE__,__LINE__)

#define KRATOS_LOG_INFO_STAMP_DETAIL ""
#define KRATOS_LOG_INFO_STAMP KRATOS_LOG_INFO_STAMP_DETAIL

#define KRATOS_LOG_DETAI_STAMP_DETAIL KRATOS_TIME_STAMP << ":"
#define KRATOS_LOG_DETAI_STAMP KRATOS_LOG_DETAI_STAMP_DETAIL

#define KRATOS_LOG_DEBUG_STAMP_DETAIL(file,line) KRATOS_TIME_STAMP << ":" << (file) << ":" << (line) << ":"
#define KRATOS_LOG_DEBUG_STAMP KRATOS_LOG_DEBUG_STAMP_DETAIL(__FILE__,__LINE__)

#define KRATOS_LOG_TRACE_STAMP_DETAIL(file,line,function) KRATOS_TIME_STAMP << KRATOS_PROCESS_ID << ":" << (file) << ":" << (line) << ":" << "\n" << "\t in function: " << function << ":\n"
#define KRATOS_LOG_TRACE_STAMP KRATOS_LOG_TRACE_STAMP_DETAIL(__FILE__,__LINE__,KRATOS_LOG_FILTER(BOOST_CURRENT_FUNCTION))

// Logs
#define KRATOS_ASSERT(C)              KRATOS_LOG_IF(KRATOS_SEVERITY_ERROR,C)        << KRATOS_LOG_ERROR_STAMP

#define KRATOS_LOG_ERROR              KRATOS_LOG(KRATOS_SEVERITY_ERROR)             << KRATOS_LOG_ERROR_STAMP
#define KRATOS_LOG_ERROR_N(N)         KRATOS_LOG_N(KRATOS_SEVERITY_ERROR,N)         << KRATOS_LOG_ERROR_STAMP
#define KRATOS_LOG_ERROR_IF(C)        KRATOS_LOG_IF(KRATOS_SEVERITY_ERROR,C)        << KRATOS_LOG_ERROR_STAMP
#define KRATOS_LOG_ERROR_IF_N(C,N)    KRATOS_LOG_IF_N(KRATOS_SEVERITY_ERROR,C,N)    << KRATOS_LOG_ERROR_STAMP
#define KRATOS_LOG_ERROR_FIRST_N(N)   KRATOS_LOG_FIRST_N(KRATOS_SEVERITY_ERROR,N)   << KRATOS_LOG_ERROR_STAMP
#define KRATOS_LOG_ERROR_CHECK(C)     KRATOS_LOG_IF(KRATOS_SEVERITY_ERROR,C)        << KRATOS_LOG_ERROR_STAMP

#define KRATOS_LOG_WARNING            KRATOS_LOG(KRATOS_SEVERITY_WARNING)           << KRATOS_LOG_WARNING_STAMP
#define KRATOS_LOG_WARNING_N(N)       KRATOS_LOG_N(KRATOS_SEVERITY_WARNING,N)       << KRATOS_LOG_WARNING_STAMP
#define KRATOS_LOG_WARNING_IF(C)      KRATOS_LOG_IF(KRATOS_SEVERITY_WARNING,C)      << KRATOS_LOG_WARNING_STAMP
#define KRATOS_LOG_WARNING_IF_N(C,N)  KRATOS_LOG_IF_N(KRATOS_SEVERITY_WARNING,C,N)  << KRATOS_LOG_WARNING_STAMP
#define KRATOS_LOG_WARNING_FIRST_N(N) KRATOS_LOG_FIRST_N(KRATOS_SEVERITY_WARNING,N) << KRATOS_LOG_WARNING_STAMP
#define KRATOS_LOG_WARNING_CHECK(C)   KRATOS_LOG_IF(KRATOS_SEVERITY_WARNING,C)      << KRATOS_LOG_WARNING_STAMP

#define KRATOS_LOG_INFO               KRATOS_LOG(KRATOS_SEVERITY_INFO)              << KRATOS_LOG_INFO_STAMP
#define KRATOS_LOG_INFO_N(N)          KRATOS_LOG_N(KRATOS_SEVERITY_INFO,N)          << KRATOS_LOG_INFO_STAMP
#define KRATOS_LOG_INFO_IF(C)         KRATOS_LOG_IF(KRATOS_SEVERITY_INFO,C)         << KRATOS_LOG_INFO_STAMP
#define KRATOS_LOG_INFO_IF_N(C,N)     KRATOS_LOG_IF_N(KRATOS_SEVERITY_INFO,C)       << KRATOS_LOG_INFO_STAMP
#define KRATOS_LOG_INFO_FIRST_N(N)    KRATOS_LOG_FIRST_N(KRATOS_SEVERITY_INFO,N)    << KRATOS_LOG_INFO_STAMP

#define KRATOS_LOG_DETAIL             KRATOS_LOG(KRATOS_SEVERITY_DETAIL)            << KRATOS_LOG_DETAIL_STAMP
#define KRATOS_LOG_DETAIL_N(N)        KRATOS_LOG_N(KRATOS_SEVERITY_DETAIL,N)        << KRATOS_LOG_DETAIL_STAMP
#define KRATOS_LOG_DETAIL_IF(C)       KRATOS_LOG_IF(KRATOS_SEVERITY_DETAIL,C)       << KRATOS_LOG_DETAIL_STAMP
#define KRATOS_LOG_DETAIL_IF_N(C,N)   KRATOS_LOG_IF_N(KRATOS_SEVERITY_DETAIL,C)     << KRATOS_LOG_DETAIL_STAMP
#define KRATOS_LOG_DETAIL_FIRST_N(N)  KRATOS_LOG_FIRST_N(KRATOS_SEVERITY_DETAIL,N)  << KRATOS_LOG_DETAIL_STAMP

#ifdef  KRATOS_DEBUG
#define KRATOS_LOG_DEBUG              KRATOS_LOG(KRATOS_SEVERITY_DEBUG)             << KRATOS_LOG_DEBUG_STAMP
#define KRATOS_LOG_DEBUG_N(N)         KRATOS_LOG_N(KRATOS_SEVERITY_DEBUG,N)         << KRATOS_LOG_DEBUG_STAMP
#define KRATOS_LOG_DEBUG_IF(C)        KRATOS_LOG_IF(KRATOS_SEVERITY_DEBUG,C)        << KRATOS_LOG_DEBUG_STAMP
#define KRATOS_LOG_DEBUG_IF_N(C,N)    KRATOS_LOG_IF_N(KRATOS_SEVERITY_DEBUG,C)      << KRATOS_LOG_DEBUG_STAMP
#define KRATOS_LOG_DEBUG_FIRST_N(N)   KRATOS_LOG_FIRST_N(KRATOS_SEVERITY_DEBUG,N)   << KRATOS_LOG_DEBUG_STAMP

#define KRATOS_LOG_TRACE              KRATOS_LOG(KRATOS_SEVERITY_TRACE)             << KRATOS_LOG_TRACE_STAMP
#define KRATOS_LOG_TRACE_N(N)         KRATOS_LOG_N(KRATOS_SEVERITY_TRACE,N)         << KRATOS_LOG_TRACE_STAMP
#define KRATOS_LOG_TRACE_IF(C)        KRATOS_LOG_IF(KRATOS_SEVERITY_TRACE,C)        << KRATOS_LOG_TRACE_STAMP
#define KRATOS_LOG_TRACE_IF_N(C,N)    KRATOS_LOG_IF_N(KRATOS_SEVERITY_TRACE,C)      << KRATOS_LOG_TRACE_STAMP
#define KRATOS_LOG_TRACE_FIRST_N(N)   KRATOS_LOG_FIRST_N(KRATOS_SEVERITY_TRACE,N)   << KRATOS_LOG_TRACE_STAMP
#else
#define KRATOS_LOG_DEBUG              
#define KRATOS_LOG_DEBUG_N(N)         
#define KRATOS_LOG_DEBUG_IF(C)        
#define KRATOS_LOG_DEBUG_IF_N(C,N)    
#define KRATOS_LOG_DEBUG_FIRST_N(N)   

#define KRATOS_LOG_TRACE              
#define KRATOS_LOG_TRACE_N(N)         
#define KRATOS_LOG_TRACE_IF(C)        
#define KRATOS_LOG_TRACE_IF_N(C,N)    
#define KRATOS_LOG_TRACE_FIRST_N(N)   
#endif

// Normal log
#define KRATOS_LOG(SEVERITY) KratosLog<SEVERITY>::GetInstance()

// Each-N block
#define KRATOS_LOG_OCCURRENCES_LINE(line) kratos_log_loop_counter##line
#define KRATOS_LOG_OCCURRENCES KRATOS_LOG_OCCURRENCES_LINE(__LINE__)

#define KRATOS_LOG_N(KRATOS_SEVERITY,__N__) \
  static int KRATOS_LOG_OCCURRENCES = 0;\
  if (++KRATOS_LOG_OCCURRENCES > __N__) KRATOS_LOG_OCCURRENCES -= __N__; \
  if (KRATOS_LOG_OCCURRENCES == 1) \
    KRATOS_LOG(KRATOS_SEVERITY)
    
// If block
#define KRATOS_LOG_IF(KRATOS_SEVERITY,__C__) \
  if (__C__) \
    KRATOS_LOG(KRATOS_SEVERITY)
    
// If-Each-N block
#define KRATOS_LOG_IF_N(KRATOS_SEVERITY,__C__,__N__) \
  static int KRATOS_LOG_OCCURRENCES = 0;\
  if (++KRATOS_LOG_OCCURRENCES > __N__) KRATOS_LOG_OCCURRENCES -= __N__; \
  if ((__C__) && KRATOS_LOG_OCCURRENCES == 1) \
    KRATOS_LOG(KRATOS_SEVERITY)
    
// First-N block
#define KRATOS_LOG_FIRST_N(KRATOS_SEVERITY,__N__) \
  static int KRATOS_LOG_OCCURRENCES = 0;\
  if (++KRATOS_LOG_OCCURRENCES < __N__) \
    KRATOS_LOG(KRATOS_SEVERITY)

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
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
  
enum KRATOS_SEVERITY {
  KRATOS_SEVERITY_ERROR,
  KRATOS_SEVERITY_WARNING,
  KRATOS_SEVERITY_INFO,
  KRATOS_SEVERITY_DETAIL,
  KRATOS_SEVERITY_DEBUG,
  KRATOS_SEVERITY_TRACE,
  KRATOS_LOG_SOLVER
};
 
template<KRATOS_SEVERITY S>  
class KratosLog : public std::ostream
{
  private:
    
    KratosLog<S>();
    
    static bool         mInstanceFlag;
    static KratosLog *  mpInstance;
    
  public:
         
    ~KratosLog<S>();
    
    static KratosLog<S>& GetInstance();
    
    void SetSeverityLevel(int severity);
    void SetOutput(std::string const& LogFileName);
    
    int getSeverity();
    std::ostream& getOutput();
    
    bool            mPrintLogOnScreen;
    int             mSeverityLevel;
    std::ofstream * mOutputFile;
    
};  

class KratosLogUtils
{
  private:
  
    KratosLogUtils();
    
    static bool mInstanceFlag;
    static KratosLogUtils * mpInstance;
    
    std::vector<std::string> * mFilterString;
    
  public:
    
    ~KratosLogUtils();
    
    static KratosLogUtils& GetInstance();
    
    void AddCustomFilter(const std::string filter);
    
    const std::string CurrentDateTime();
    const std::string CurrentDateTime(const char * format);
    const char * Filter(const char * inputString);
    const char * FilterNamespace(const char * inputString);
    const char * FilterCustom(const char * inputString, std::vector<std::string> remove);
    
};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{
  
// Generic Specialization
template<KRATOS_SEVERITY S>
KratosLog<S>::KratosLog()
{
    mOutputFile = new std::ofstream();
}

template<KRATOS_SEVERITY S>
KratosLog<S>::~KratosLog()
{
    free(mOutputFile);
}

template<KRATOS_SEVERITY S>
void KratosLog<S>::SetOutput(std::string const& LogFileName)
{
    if(mOutputFile->is_open())

    mOutputFile->close();

    mOutputFile->open(LogFileName.c_str());
}

template<KRATOS_SEVERITY S>
KratosLog<S>& KratosLog<S>::GetInstance() 
{
    if(!mInstanceFlag)
    {   
        mInstanceFlag = 1;
        mpInstance = new KratosLog<S>();
        mpInstance->SetOutput("KratosLog.log");
    }
    
    return (*mpInstance);
}

// Error specialization
inline KratosLog<KRATOS_SEVERITY_ERROR> & operator << (KratosLog<KRATOS_SEVERITY_ERROR> &out, std::ostream& (*fn)(std::ostream&))
{ 
  if(out.mPrintLogOnScreen)
    fn(std::cout);
  
  fn((*out.mOutputFile));
  (*out.mOutputFile).flush();
  
    abort();
  //TODO: Throw exception abort();

  return out;
}
      
template<class T>
inline KratosLog<KRATOS_SEVERITY_ERROR> & operator << (KratosLog<KRATOS_SEVERITY_ERROR> &out, const T& data)
{
  if(out.mPrintLogOnScreen)
  {
      std::cout << data;
  }

  (*out.mOutputFile) << data;
  
  //TODO: Throw exception abort();
  
  return out;  
}
    
// Error specialization
template<KRATOS_SEVERITY S>
inline KratosLog<S> & operator << (KratosLog<S> &out, std::ostream& (*fn)(std::ostream&))
{ 
  if(out.mPrintLogOnScreen)
    fn(std::cout);
  
  fn((*out.mOutputFile));
  (*out.mOutputFile).flush();

  return out;
}
  
template<class T, KRATOS_SEVERITY S>
inline KratosLog<S> & operator << (KratosLog<S> &out, const T& data)
{
  if(out.mPrintLogOnScreen)
  {
      std::cout << data;
  }

  (*out.mOutputFile) << data;
  
  return out;  
}

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_LOG_H_INCLUDED  defined 


