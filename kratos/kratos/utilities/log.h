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
#include <map>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"

// Defines
int const KratosError      = 0; 
int const KratosWarning    = 1; 
int const KratosInfo       = 2;
int const KratosDebug      = 3; 
int const KratosTrace      = 4;

#define KRATOS_LOG(severity) Log::severity##_stream(__FILE__,BOOST_CURRENT_FUNCTION,__LINE__,severity)

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
    
    Log() 
    {    
        std::stringstream LogFileName;
        
        std::stringstream LogStandard;
        std::stringstream LogWarning;
        std::stringstream LogError;
        
        mPrintLogOnScreen = 0;
        mPrintWngOnScreen = 0;
        mPrintErrOnScreen = 1;
        
        //Default filename for logs
        LogFileName << "KratosMultiphysics" << '-' << getpid();
        
        LogStandard << LogFileName << '.' << "log";
        LogWarning  << LogFileName << '.' << "wng";
        LogError    << LogFileName << '.' << "err";
        
        SetLogFile(LogFileName.str());
        SetWngFile(LogFileName.str());
        SetErrFile(LogFileName.str());
    }
    
public:
  
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Timer
    KRATOS_CLASS_POINTER_DEFINITION(Log);

    ///@}
    ///@name Life Cycle
    ///@{
  
    static Log& GetInstance() 
    {
        if(!mInstanceFalg)
            mInstance = new Log();
        
        return (*mInstance);
    }
    
    static std::ostream& KRATOS_ERROR_stream(const char * file, const char * funct, int line, int severity)
    {
        msErrFile << "KRATOS_ERROR:" << file << ":" << line << " In function \'" << funct << "\':";
        return msErrFile;
    }

    /// Destructor.
    virtual ~Log()
    {
        msLogFile << std::endl;
      
        if(msLogFile.is_open())
            msLogFile.close();
            
        delete mInstance;
    }


    ///@}
    ///@name Operators
    ///@{
      
    // NOT NEEDED ANYMORE!!!
    /// Class operator only writes in the generic log file 
    
    /// This function must be overloaded twice in order to adress std::endl as it is a function
//     Log& operator<<(std::ostream& (*fn)(std::ostream&))
//     {
//         fn(msLogFile);
//         
//         if(mPrintLogOnScreen)
//             fn(std::cout);
//         
//         return * this;
//     }
//     
//     /// Rest of the inputs
//     template<typename T>
//     Log& operator << (const T& data)
//     {
//         msLogFile << data;
//         
//         if(mPrintLogOnScreen)
//             std::cout << msLogFile;
//         
//         return * this;
//     }

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
    
    int SetWngFile(std::string const& WngFileName)
    {
        if(msWngFile.is_open())
            msWngFile.close();

        msWngFile.open(WngFileName.c_str());

        return msWngFile.is_open();
    }
    
    int SetErrFile(std::string const& ErrFileName)
    {
        if(msErrFile.is_open())
            msErrFile.close();

        msErrFile.open(ErrFileName.c_str());

        return msErrFile.is_open();
    }
    
    void SetPrintLogOnScreen(bool print) 
    {
        mPrintLogOnScreen = print;
    }
    
    void SetPrintWngOnScreen(bool print) 
    {
        mPrintLogOnScreen = print;
    }
    
    void SetPrintErrOnScreen(bool print) 
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
    
    std::ofstream& GetWngStream()
    {
        return msWngFile;
    }
    
    std::ofstream& GetErrStream()
    {
        return msErrFile;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
      
    //TODO: Ask Pooyna if we allow write of multiple lines or just a line per call
      
    template<typename T>
    void WriteLog(const T& data) 
    {
        msLogFile << data;
        
        if(mPrintLogOnScreen)
            std::cout << msLogFile;
    }
    
    template<typename T>
    void WriteWng(const T& data)
    {
        msWngFile << data;
        
        if(mPrintWngOnScreen)
            std::cout << msWngFile;
    }
    
    template<typename T>
    void WriteErr(const T& data)
    {
        msErrFile << data;
        
        if(mPrintErrOnScreen)
            std::cout << msErrFile;      
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Log";
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

    static std::ofstream msLogFile;
    static std::ofstream msWngFile;
    static std::ofstream msErrFile;

    bool mPrintLogOnScreen;
    bool mPrintWngOnScreen;
    bool mPrintErrOnScreen;

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


