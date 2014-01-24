/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
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
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"

// Macros
#define KRATOS_ERROR_LOG            KRATOS_LOG(KRATOS_ERROR)
#define KRATOS_ERROR_LOG_N          KRATOS_LOG_N(KRATOS_ERROR)
#define KRATOS_ERROR_LOG_IF         KRATOS_LOG_IF(KRATOS_ERROR)
#define KRATOS_ERROR_LOG_CHECK      KRATOS_LOG_CHECK(KRATOS_ERROR)

#ifdef KRATOS_DEBUG
#define KRATOS_DEBUG_LOG    KRATOS_LOG(KRATOS_DEBUG)
#define KRATOS_TRACE_LOG    KRATOS_LOG(KRATOS_TRACE)
#else
#define KRATOS_DEBUG_LOG    
#define KRATOS_TRACE_LOG    
#endif

#define KRATOS_LOG(SEVERITY) KratosLogger<SEVERITY>::KratosLogStream()

# define LOG_OCCURRENCES(line) #KRATOS_LOG_N_ ## line

// Each N
#define KRATOS_LOG_N(KRATOS_SEVERITY) \
  {static int LOG_OCCURRENCES(__FILE__,__LINE__) = 0, LOG_OCCURRENCES_MOD_N = 0; \
  ++LOG_OCCURRENCES; \
  if (++LOG_OCCURRENCES_MOD_N > n) LOG_OCCURRENCES_MOD_N -= n; \
  if (LOG_OCCURRENCES_MOD_N == 1) \
    KRATOS_LOG(KRATOS_SEVERITY)}

#define KRATOS_TIME_STAMP   "[" << CurrentDateTime() << "]"

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
 
template<int S>
class KratosLog
{
  public:
    
    static SetOutput(int severity);
    static SetSeverityLevel(std::string const& LogFileName);

    Log& operator<<(std::ostream& (*fn)(std::ostream&));
    
    template<typename T>
    Log& operator << (const T& data);
    
  private:
    
    static int mSeverity;
    static std::ostream mOutputFile;
 
  
};  
  
/// Logger clas for kratos
/** Detail class definition.
*/
class KratosLog : public std::ostream
{
  private:
    
    /**
     * Default constructor
     */
    Log();
    
  public:
  
    ///@name Type Definitions
    ///@{
      
    enum{
      KRATOS_SEVERITY_ERR,
      KRATOS_SEVERITY_WAR,
      KRATOS_SEVERITY_DET,
      KRATOS_SEVERITY_INF,
      KRATOS_SEVERITY_DEB,
      KRATOS_SEVERITY_TRC
    }

    /// Pointer definition of Timer
    KRATOS_CLASS_POINTER_DEFINITION(KratosLog);

    ///@}
    ///@name Life Cycle
    ///@{
  
    /**
     * Returns the instance to the log class and creates it in case it not been yet created 
     */
    static KratosLogger& GetInstance();
    static Kratos

    /// Destructor.
    /*virtual*/ ~KratosLogger()
    {
    }
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
      
    int SetLogFile(std::string const& LogFileName);
   
    void SetPrintLogOnScreen(bool print) ;

    ///@}
    ///@name Access
    ///@{
      
    KratosLog& GetKratosErrorLog();
    KratosLog& GetKratosWarningLog();
    KratosLog& GetKratosDetailLog();
    KratosLog& GetKratosInfoLog();
    KratosLog& GetKratosDebugLog();
    KratosLog& GetKratosTraceLog();
    
    static const std::string CurrentDateTime();
    static const std::string CurrentDateTime(const char * format);

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Kratos Logger";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
//         PrintTimingInformation(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
      
    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

  private:

    ///@name Static Member Variables
    ///@{

    static bool           mInstanceFalg;
    static KratosLogger * mpInstance;
    
    static KratosLog<> mError

    std::ofstream msNormalLogFile;
    std::ofstream msTimingLogFile;
    std::ofstream msSolverLogFile;
    std::ofstream msExtrasLogFile;

    bool mPrintLogOnScreen;
    int  mSeverityLevel;

    ///@}
    ///@name Member Variables
    ///@{


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
  
// template<typename T>
// Log& operator << (Log& log, const T& data)
// {
//     log.GetOutputStream() << data;
//     
//     return log;
// }

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   Log& rThis)
// {
//     return rIStream;
// }

/// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const Log& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }
///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_LOG_H_INCLUDED  defined 


