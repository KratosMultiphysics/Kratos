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
  
    /// KratosLogger
  
    KratosLogger::KratosLogger()
    {    
        std::stringstream NormalLogFileName;
        std::stringstream TimingLogFileName;
        std::stringstream SolverLogFileName;
        std::stringstream ExtrasLogFileName;
        
        //mSeverityLevel    = KRATOS_INFO;
        
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

    KratosLogger::~KratosLogger()
    {
        msNormalLogFile << std::endl;
      
        if(msNormalLogFile.is_open())
            msNormalLogFile.close();
        
        mInstanceFalg = 0;
            
        delete mpInstance;
    }

    static KratosLogger& KratosLogger::GetInstance() 
    {
        if(!mInstanceFalg)
        {   
            mInstanceFalg = 1;
            mpInstance = new KratosLogger();
        }
        
        return (*mpInstance);
    }
    
    int KratosLogger::SetLogFile(std::ofstream & fileStream, std::string const& LogFileName)
    {
        if(fileStream.is_open())
            fileStream.close();

        msNormalLogFile.open(LogFileName.c_str());

        return fileStream.is_open();
    }
    
//     void KratosLogger::SetPrintLogOnScreen(bool print) 
//     {
//         mPrintLogOnScreen = print;
//     }
    
//     void KratosLogger::SetSeverityLevel(int severity)
//     {
//         mSeverityLevel = severity;
//     }

    std::ofstream& KratosLogger::GetLogNormalFile()
    {
        return msNormalLogFile;
    }
    
    std::ofstream& KratosLogger::GetLogTimingFile()
    {
        return msTimingLogFile;
    }
    
    std::ofstream& KratosLogger::GetLogSolverFile()
    {
        return msSolverLogFile;
    }
    
    std::ofstream& KratosLogger::GetLogExtrasFile()
    {
        return msExtrasLogFile;
    }
    
    /////////////// KRATOS LOG ////////////////////////////////////////////////////////////////

//     void KratosLog::SetSeverityLevel(int severity)
//     {
//         //mSeverity = severity;
//     }
//     
//     void KratosLog::SetOutput(std::string const& LogFileName)
//     {
//         //mOutputFile = 
//     }

    
    /// Class operator only writes in the generic log file 
    /// This function must be overloaded twice in order to adress std::endl as it is a function
    KratosLog& KratosLog::operator<<(std::ostream& (*fn)(std::ostream&))
    {   
        if(mPrintLogOnScreen)
            fn(std::cout);
     
        fn(mOutputFile);
        
        return * this;
    }
    
    /// Rest of the inputs
    template<typename T>
    KratosLog& KratosLog::operator << (const T& data)
    {
        if(mPrintLogOnScreen)
        {
            std::cout << data;
        }
        
        mOutputFile << data;
        
        if(mSeverityLevel == KratosLogger::KRATOS_SEVERITY_ERR)
        {
            abort();
        }
        
        return * this;
    }
    
    /// KratosLogger
  
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    template<KratosLogger::KRATOS_SEVERITY_ERR>
    static inline KratosLog& KratosLogger::KratosLogStream(const char * file, const char * funct, int line, int severity)
    {
        KratosLog<KratosLogger::KRATOS_SEVERITY_ERR> << KRATOS_TIME_STAMP << ":" << "KRATOS_ERROR:" << file << ":" << line << " In function \'" << funct << "\':";
        return KratosLog<KratosLogger::KRATOS_SEVERITY_ERR>::mOutputFile;
    }
    
//     /**
//      * @param file:     File where the function was called
//      * @param funct:    Function where the function was called
//      * @param line:     Line where the function was called
//      * @param severity: Severity of the error **TODO: To be removed?
//      **/ 
//     static inline KratosLog& KratosLogger::KratosLogStream<KRATOS_SEVERITY_WAR>(const char * file, const char * funct, int line, int severity)
//     {
//         KratosLog& stream = KratosLog::GetInstance();
//             
//         stream.GetLogFile() << KRATOS_TIME_STAMP << ":" << "KRATOS_WARNING:" << file << ":" << line << "\':";
//         return stream;
//     }
//     
//     /**
//      * @param file:     File where the function was called
//      * @param funct:    Function where the function was called
//      * @param line:     Line where the function was called
//      * @param severity: Severity of the error **TODO: To be removed?
//      **/
//     static inline KratosLog& KratosLogger::KratosLogStream<KRATOS_SEVERITY_DET>(const char * file, const char * funct, int line, int severity)
//     {
//         KratosLog& stream = KratosLog::GetInstance();
//          
//         stream.GetLogFile() << KRATOS_TIME_STAMP << ":";
//         return stream;
//     }
//     
//     /**
//      * @param file:     File where the function was called
//      * @param funct:    Function where the function was called
//      * @param line:     Line where the function was called
//      * @param severity: Severity of the error **TODO: To be removed?
//      **/
//     static inline KratosLog& KratosLogger::KratosLogStream<KRATOS_SEVERITY_INF>(const char * file, const char * funct, int line, int severity)
//     {
//         KratosLog& stream = KratosLog::GetInstance();
//          
//         stream.GetLogFile() << KRATOS_TIME_STAMP << ":";
//         return stream;
//     }
//     
//     /**
//      * @param file:     File where the function was called
//      * @param funct:    Function where the function was called
//      * @param line:     Line where the function was called
//      * @param severity: Severity of the error **TODO: To be removed?
//      **/
//     static inline KratosLog& KratosLogger::KratosLogStream<KRATOS_SEVERITY_DEB>(const char * file, const char * funct, int line, int severity)
//     {
//         KratosLog& stream = KratosLog::GetInstance();
//         
//         stream.GetLogFile() << KRATOS_TIME_STAMP << ":";
//         return KratosLog::GetInstance();
//     }
//     
//     /**
//      * @param file:     File where the function was called
//      * @param funct:    Function where the function was called
//      * @param line:     Line where the function was called
//      * @param severity: Severity of the error **TODO: To be removed?
//      **/
//     static inline KratosLog& KratosLogger::KratosLogStream<KRATOS_SEVERITY_TRC>(const char * file, const char * funct, int line, int severity)
//     {
//         KratosLog& stream = KratosLog::GetInstance();
//         
//         stream.GetLogFile() << KRATOS_TIME_STAMP << ":" << "KRATOS_TRACE:" << file << ":" << line << " In function \'" << funct << "\':";
//         return stream;
//     }
    
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
