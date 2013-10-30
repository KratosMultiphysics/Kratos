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


#if !defined(KRATOS_LOG_H_INCLUDED )
#define  KRATOS_LOG_H_INCLUDED

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

// Defines
int const KRATOS_ERROR      = 0; 
int const KRATOS_WARNING    = 1; 
int const KRATOS_INFO       = 2;
int const KRATOS_DEBUG      = 3; 
int const KRATOS_TRACE      = 4;

#define KRATOS_LOG(severity) Logger::severity##_stream(__FILE__,BOOST_CURRENT_FUNCTION,__LINE__,severity)

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
  
/// Logger clas for kratos
/** Detail class definition.
*/
class Log : public std::ostream
{
private:
  
    static bool     mInstanceFalg;
    static Log *    mInstance;
    
    /**
     * Default constructor
     */
    Log() 
    {    
        std::stringstream LogFileName;
        
        mPrintLogOnScreen = 0;
        mSeverityLevel    = KRATOS_INFO;
        
        LogFileName << "KratosMultiphysics" << '.' << getpid() << '.' << "log";
        
        SetLogFile(LogFileName.str());
    }
    
public:
  
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Timer
    KRATOS_CLASS_POINTER_DEFINITION(Log);

    ///@}
    ///@name Life Cycle
    ///@{
  
    /**
     * Returns the instance to the log class and creates it in case it not been yet created 
     */
    static Log& GetInstance() 
    {
        if(!mInstanceFalg)
        {   
            mInstanceFalg = 1;
            mInstance = new Log();
        }
        
        return (*mInstance);
    }

    /// Destructor.
    virtual ~Log()
    {
        msLogFile << std::endl;
      
        if(msLogFile.is_open())
            msLogFile.close();
        
        mInstanceFalg = 0;
            
        delete mInstance;
    }
    
    ///@}
    ///@name Operators
    ///@{
      
    // NOT NEEDED ANYMORE!!!
    /// Class operator only writes in the generic log file 
 
    /// This function must be overloaded twice in order to adress std::endl as it is a function
    Log& operator<<(std::ostream& (*fn)(std::ostream&))
    {
        fn(msLogFile);
        
        if(mPrintLogOnScreen)
            fn(std::cout);
        
        return * this;
    }
    
    /// Rest of the inputs
    template<typename T>
    Log& operator << (const T& data)
    {
        msLogFile << data;
        
        if(mPrintLogOnScreen)
            std::cout << msLogFile;
        
        if(mSeverityLevel == KRATOS_ERROR)
            abort();
        
        return * this;
    }

    ///@}
    ///@name Operations
    ///@{
      
    int SetLogFile(std::string const& LogFileName)
    {
        if(msLogFile.is_open())
            msLogFile.close();

        msLogFile.open(LogFileName.c_str());

        return msLogFile.is_open();
    }
    
    void SetPrintLogOnScreen(bool print) 
    {
        mPrintLogOnScreen = print;
    }

    ///@}
    ///@name Access
    ///@{
    
    std::ofstream& GetLogStream()
    {
        return msLogFile;
    }

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

    std::ofstream msLogFile;

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

}; // Class Log

class Logger
{
public:
  
    /**
     * Preamble for log output. This version should never be called.
     * 
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& _stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
         
        stream << "[" << CurrentDateTime() << "]:" << "MISSING_SEVERITY: " << "\':";
        return stream;
    }
  
    /**
     * Preamble for log output. This version prints:
     * 
     * Timestamp
     * Severity
     * file
     * line
     * function
     * 
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KRATOS_ERROR_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
          
        stream << "[" << CurrentDateTime() << "]:" <<  "KRATOS_ERROR:" << file << ":" << line << " In function \'" << funct << "\':";
        return stream;
    }
    
    /**
     * Preamble for log output. This version should never be called.
     * 
     * Timestamp
     * Severity
     * file
     * line
     * 
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/ 
    static Log& KRATOS_WARNING_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
            
        stream << "[" << CurrentDateTime() << "]:" <<  "KRATOS_WARNING:" << file << ":" << line << "\':";
        return stream;
    }
    
    /**
     * Preamble for log output. This version should never be called.
     * 
     * Timestamp
     * 
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KRATOS_INFO_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
         
        stream << "[" << CurrentDateTime() << "]:";
        return stream;
    }
    
    /**
     * Preamble for log output. This version should never be called.
     * 
     * Timestamp
     * 
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KRATOS_DEBUG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
        
        stream << "[" << CurrentDateTime() << "]:";
        return Log::GetInstance();
    }
    
    /**
     * Preamble for log output. This version print trace-level information
     * 
     * Timestamp
     * Severity
     * file
     * line
     * function
     * 
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KRATOS_TRACE_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
        
        stream << "[" << CurrentDateTime() << "]:" <<  "KRATOS_TRACE:" << file << ":" << line << " In function \'" << funct << "\':";
        return stream;
    }
    
private:
  
    /**
     * Returns the current time in HH:MM:SS format
     */
    static const std::string CurrentDateTime() 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        
        char       buf[80];
        tstruct = *localtime(&now);
        
        strftime(buf, sizeof(buf), "%X", &tstruct);

        return buf;
    }
    
    /**
    * Returns the current time in the specified format
    * @param format:    valid format for time. Please check: "http://www.cplusplus.com/reference/ctime/strftime/"
    *                   for more details.
    */
    static const std::string CurrentDateTime(const char * format) 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        
        char       buf[80];
        tstruct = *localtime(&now);
        
        strftime(buf, sizeof(buf), format, &tstruct);

        return buf;
    }
  
};

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

#endif // KRATOS_LOG_H_INCLUDED  defined 


