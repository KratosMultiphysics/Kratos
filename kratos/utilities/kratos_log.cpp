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
  
    template<> bool                             KratosLog<KRATOS_SEVERITY_ERROR>::mInstanceFalg = 0;
    template<> KratosLog<KRATOS_SEVERITY_ERROR> * KratosLog<KRATOS_SEVERITY_ERROR>::mpInstance    = NULL;
    
    template<> KratosLog<KRATOS_SEVERITY_ERROR>::KratosLog()
    {
        mOutputFile = new std::ofstream();
    }
    
    template<> KratosLog<KRATOS_SEVERITY_ERROR>::~KratosLog()
    {
        free(mOutputFile);
    }

    template<> void KratosLog<KRATOS_SEVERITY_ERROR>::SetOutput(std::string const& LogFileName)
    {
        if(mOutputFile->is_open())

            mOutputFile->close();

        mOutputFile->open(LogFileName.c_str());
    }
    
    template<> KratosLog<KRATOS_SEVERITY_ERROR>& KratosLog<KRATOS_SEVERITY_ERROR>::GetInstance() 
    {
        if(!mInstanceFalg)
        {   
            mInstanceFalg = 1;
            mpInstance = new KratosLog<KRATOS_SEVERITY_ERROR>();
            mpInstance->SetOutput("KRATOS_LOG_TEST.log");
        }
        
        return (*mpInstance);
    }
  
    /**
     * @param file:     File where the function was called
     * @param funct:    Function where the function was called
     * @param line:     Line where the function was called
     * @param severity: Severity of the error **TODO: To be removed?
     **/
    //inline KratosLog<KRATOS_SEVERITY_ERROR>& KratosLogger::KratosLogStream(const char * file, const char * funct, int line, int severity)
    //{
        //KratosLog<KRATOS_SEVERITY_ERROR>& log = KratosLog<KRATOS_SEVERITY_ERROR>::GetInstance();
        
        //log << KRATOS_TIME_STAMP << ":" << "KRATOS_ERROR:" << file << ":" << line << " In function \'" << funct << "\':";
        
        //return log;
    //}
    
    
//         /**
//      * Returns the current time in HH:MM:SS format
//      */
//     const std::string KratosLogger::CurrentDateTime() 
//     {
//         time_t     now = time(0);
//         struct tm  tstruct;
//         
//         char       buf[80];
//         tstruct = *localtime(&now);
//         
//         strftime(buf, sizeof(buf), "%X", &tstruct);
// 
//         return buf;
//     }
//     
//     /**
//     * Returns the current time in the specified format
//     * @param format:    valid format for time. Please check: "http://www.cplusplus.com/reference/ctime/strftime/"
//     *                   for more details.
//     */
//     const std::string KratosLogger::CurrentDateTime(const char * format) 
//     {
//         time_t     now = time(0);
//         struct tm  tstruct;
//         
//         char       buf[80];
//         tstruct = *localtime(&now);
//         
//         strftime(buf, sizeof(buf), format, &tstruct);
// 
//         return buf;
//     }
}
