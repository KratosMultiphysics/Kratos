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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 11:35:13 $
//   Revision:            $Revision: 1.6 $
//
//

#include "utilities/kratos_log.h"


namespace Kratos
{
  
    /// KratosLog
  
    KratosLog()
    {    
        std::stringstream NormalLogFileName;
        std::stringstream TimingLogFileName;
        std::stringstream SolverLogFileName;
        std::stringstream ExtrasLogFileName;
        
        mPrintLogOnScreen = 0;
        mSeverityLevel    = KRATOS_INFO;
        
        // LogFileName << "KratosMultiphysics" << '.' << getpid() << '.' << "log";
        NormalLogFileName << "KratosMultiphysics" << '.' << "log";
        TimingLogFileName << "KratosMultiphysics" << '.' << "Timing" << '.' << "log";
        SolverLogFileName << "KratosMultiphysics" << '.' << "Solver" << '.' << "log";
        ExtrasLogFileName << "KratosMultiphysics" << '.' << "Extras" << '.' << "log";
        
        SetLogFile(msNormalLogFile,NormalLogFileName.str());
        SetLogFile(msTimingLogFile,TimingLogFileName.str());
        SetLogFile(msSolverLogFile,SolverLogFileName.str());
        SetLogFile(msExtrasLogFile,ExtrasLogFileName.str());
    }

    ~KratosLog()
    {
        msLogFile << std::endl;
      
        if(msLogFile.is_open())
            msLogFile.close();
        
        mInstanceFalg = 0;
            
        delete mInstance;
    }

    static KratosLog& KratosLog::GetInstance() 
    {
        if(!mInstanceFalg)
        {   
            mInstanceFalg = 1;
            mInstance = new KratosLog();
        }
        
        return (*mInstance);
    }
    
    int KratosLog::SetLogFile(std::ofstream & fileStream, std::string const& LogFileName)
    {
        if(fileStream.is_open())
            fileStream.close();

        msNormalLogFile.open(LogFileName.c_str());

        return fileStream.is_open();
    }
    
    void KratosLog::SetPrintLogOnScreen(bool print) 
    {
        mPrintLogOnScreen = print;
    }
    
    void KratosLog::SetSeverityLevel(int severity)
    {
        mSeverityLevel = severity;
    }
    
    std::ofstream& GetNormalLogFile()
    {
        return msNormalLogFile;
    }
    
    std::ofstream& GetTimingLogFile()
    {
        return msTimingLogFile;
    }
    
    std::ofstream& GetSolverLogFile()
    {
        return msSolverLogFile;
    }
    
    std::ofstream& GetExtrasLogFile()
    {
        return msExtrasLogFile;
    }
    
    /// Class operator only writes in the generic log file 
    /// This function must be overloaded twice in order to adress std::endl as it is a function
    Log& KratosLog::operator<<(std::ostream& (*fn)(std::ostream&))
    {   
        if(mPrintLogOnScreen)
            fn(std::cout);
     
        fn(msLogFile);
        
        return * this;
    }
    
    /// Rest of the inputs
    template<typename T>
    Log& KratosLog::operator << (const T& data)
    {
        if(mPrintLogOnScreen)
        {
            std::cout << msLogFile;
        }
        
        msLogFile << data;
        
        if(mSeverityLevel == KRATOS_ERROR)
        {
            abort();
        }
        
        return * this;
    }
    
    bool Log::mInstanceFalg = 0;
    Log* Log::mInstance     = NULL;
    
    /// KratosLogger
    
    /**
     * This should never be called or instanced.
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KratosLogger::_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
         
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":" << "MISSING_SEVERITY: " << "\':";
        return stream;
    }
  
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KratosLogger::KRATOS_ERROR_LOG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
          
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":" << "KRATOS_ERROR:" << file << ":" << line << " In function \'" << funct << "\':";
        return stream;
    }
    
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/ 
    static Log& KratosLogger::KRATOS_WARNING_LOG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
            
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":" << "KRATOS_WARNING:" << file << ":" << line << "\':";
        return stream;
    }
    
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KratosLogger::KRATOS_DETAIL_LOG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
         
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":";
        return stream;
    }
    
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KratosLogger::KRATOS_INFO_LOG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
         
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":";
        return stream;
    }
    
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KratosLogger::KRATOS_DEBUG_LOG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
        
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":";
        return Log::GetInstance();
    }
    
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    static Log& KratosLogger::KRATOS_TRACE_LOG_stream(const char * file, const char * funct, int line, int severity)
    {
        Log& stream = Log::GetInstance();
        
        stream.GetLogFile() << KRATOS_TIME_STAMP << ":" << "KRATOS_TRACE:" << file << ":" << line << " In function \'" << funct << "\':";
        return stream;
    }
    
        /**
     * Returns the current time in HH:MM:SS format
     */
    static const std::string KratosLogger::CurrentDateTime() 
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
    static const std::string KratosLogger::CurrentDateTime(const char * format) 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        
        char       buf[80];
        tstruct = *localtime(&now);
        
        strftime(buf, sizeof(buf), format, &tstruct);

        return buf;
    }
}
