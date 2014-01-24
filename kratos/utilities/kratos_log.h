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

#define KRATOS_LOG(SEVERITY) KratosLog<SEVERITY>::KratosLogStream()

# define KRATOS_LOG_OCCURRENCES(line) kratos_log_loop_counter ## line

// Each N
#define KRATOS_LOG_N(KRATOS_SEVERITY) \
  {static int KRATOS_LOG_OCCURRENCES(__LINE__) = 0, KRATOS_LOG_OCCURRENCES_MOD_N = 0; \
  ++LOG_OCCURRENCES; \
  if (++KRATOS_LOG_OCCURRENCES > n) KRATOS_LOG_OCCURRENCES -= n; \
  if (KRATOS_LOG_OCCURRENCES == 1) \
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
    
    static void SetOutput(int severity);
    static void SetSeverityLevel(std::string const& LogFileName);

    KratosLog& operator<<(std::ostream& (*fn)(std::ostream&));
    
    template<typename T>
    KratosLog& operator << (const T& data);
    
    static std::ostream mOutputFile;
    
  private:
    
    static int  mSeverity;
    
    static bool mPrintLogOnScreen;
    static int  mSeverityLevel;
   
 
  
};  
  
/// Logger clas for kratos
/** Detail class definition.
*/
class KratosLogger : public std::ostream
{
  private:
    
    /**
     * Default constructor
     */
    KratosLogger();
    
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
    };

    /// Pointer definition of Timer
    KRATOS_CLASS_POINTER_DEFINITION(KratosLogger);

    ///@}
    ///@name Life Cycle
    ///@{
  
    /**
     * Returns the instance to the log class and creates it in case it not been yet created 
     */
    static KratosLogger& GetInstance();

    /// Destructor.
    ~KratosLogger();
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
      
    int SetLogFile(std::ofstream&, std::string const&);
   
    void SetPrintLogOnScreen(bool print) ;

    ///@}
    ///@name Access
    ///@{
      
    KratosLog<KRATOS_SEVERITY_ERR>& GetKratosLog();
    
    template<int S>
    KratosLog& KratosLogStream(const char * file, const char * funct, int line, int severity);
    
    std::ofstream& GetLogNormalFile();
    std::ofstream& GetLogTimingFile();
    std::ofstream& GetLogSolverFile();   
    std::ofstream& GetLogExtrasFile();
    
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

    std::ofstream msNormalLogFile;
    std::ofstream msTimingLogFile;
    std::ofstream msSolverLogFile;
    std::ofstream msExtrasLogFile;

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


