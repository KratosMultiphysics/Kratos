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
    template<> bool                                  KratosLog<KRATOS_SEVERITY_ERROR>::mInstanceFalg    = 0;
    template<> KratosLog<KRATOS_SEVERITY_ERROR> *    KratosLog<KRATOS_SEVERITY_ERROR>::mpInstance       = NULL;
    
    template<> bool                                  KratosLog<KRATOS_SEVERITY_WARNING>::mInstanceFalg  = 0;
    template<> KratosLog<KRATOS_SEVERITY_WARNING> *  KratosLog<KRATOS_SEVERITY_WARNING>::mpInstance     = NULL;
    
    template<> bool                                  KratosLog<KRATOS_SEVERITY_INFO>::mInstanceFalg     = 0;
    template<> KratosLog<KRATOS_SEVERITY_INFO> *     KratosLog<KRATOS_SEVERITY_INFO>::mpInstance        = NULL;
    
    template<> bool                                  KratosLog<KRATOS_SEVERITY_DETAIL>::mInstanceFalg   = 0;
    template<> KratosLog<KRATOS_SEVERITY_DETAIL> *   KratosLog<KRATOS_SEVERITY_DETAIL>::mpInstance      = NULL;
    
    template<> bool                                  KratosLog<KRATOS_SEVERITY_DEBUG>::mInstanceFalg    = 0;
    template<> KratosLog<KRATOS_SEVERITY_DEBUG> *    KratosLog<KRATOS_SEVERITY_DEBUG>::mpInstance       = NULL;
    
    template<> bool                                  KratosLog<KRATOS_SEVERITY_TRACE>::mInstanceFalg    = 0;
    template<> KratosLog<KRATOS_SEVERITY_TRACE> *    KratosLog<KRATOS_SEVERITY_TRACE>::mpInstance       = NULL;
    
    /**
    * Returns the current time in HH:MM:SS format
    */
    const std::string KratosLogUtils::CurrentDateTime() 
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
    const std::string KratosLogUtils::CurrentDateTime(const char * format) 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        
        char       buf[80];
        tstruct = *localtime(&now);
        
        strftime(buf, sizeof(buf), format, &tstruct);

        return buf;
    }
    
    std::vector<std::string> KratosLogUtils::GetFilterNamespaces()
    {
        std::vector<std::string> filters(2);
        
        filters.push_back("Kratos");
        filters.push_back("Kernel");
        
        return filters;
    }
    
    const char * KratosLogUtils::FilterNamespace(const char * inputString, std::vector<std::string> remove)
    {
        std::string s(inputString);
      
        //for(unsigned int i = 0; i < remove.size(); i++)
        s = s.substr(s.find("Kratos::"));

        return s.c_str();
    }
}
