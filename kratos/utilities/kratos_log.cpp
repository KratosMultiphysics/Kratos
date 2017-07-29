//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    janosch
//                    
//

#include "utilities/kratos_log.h"

namespace Kratos
{
    template<> bool                                 KratosLog<KRATOS_SEVERITY_ERROR>::mInstanceFlag    = 0;
    template<> KratosLog<KRATOS_SEVERITY_ERROR> *   KratosLog<KRATOS_SEVERITY_ERROR>::mpInstance       = NULL;
    
    template<> bool                                 KratosLog<KRATOS_SEVERITY_WARNING>::mInstanceFlag  = 0;
    template<> KratosLog<KRATOS_SEVERITY_WARNING> * KratosLog<KRATOS_SEVERITY_WARNING>::mpInstance     = NULL;
    
    template<> bool                                 KratosLog<KRATOS_SEVERITY_INFO>::mInstanceFlag     = 0;
    template<> KratosLog<KRATOS_SEVERITY_INFO> *    KratosLog<KRATOS_SEVERITY_INFO>::mpInstance        = NULL;
    
    template<> bool                                 KratosLog<KRATOS_SEVERITY_DETAIL>::mInstanceFlag   = 0;
    template<> KratosLog<KRATOS_SEVERITY_DETAIL> *  KratosLog<KRATOS_SEVERITY_DETAIL>::mpInstance      = NULL;
    
    template<> bool                                 KratosLog<KRATOS_SEVERITY_DEBUG>::mInstanceFlag    = 0;
    template<> KratosLog<KRATOS_SEVERITY_DEBUG> *   KratosLog<KRATOS_SEVERITY_DEBUG>::mpInstance       = NULL;
    
    template<> bool                                 KratosLog<KRATOS_SEVERITY_TRACE>::mInstanceFlag    = 0;
    template<> KratosLog<KRATOS_SEVERITY_TRACE> *   KratosLog<KRATOS_SEVERITY_TRACE>::mpInstance       = NULL;
    
    bool             KratosLogUtils::mInstanceFlag    = 0;
    KratosLogUtils * KratosLogUtils::mpInstance       = NULL;
    
    KratosLogUtils& KratosLogUtils::GetInstance() 
    {
        if(!mInstanceFlag)
        {   
            mInstanceFlag = 1;
            mpInstance = new KratosLogUtils();
        }
        
        return (*mpInstance);
    }
    
    KratosLogUtils::KratosLogUtils()
    {
        mFilterString = new std::vector<std::string>(0);
    }
    
    KratosLogUtils::~KratosLogUtils()
    {
        delete mFilterString;
    }
    
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
    * @param format:        valid format for time. Please check: "http://www.cplusplus.com/reference/ctime/strftime/"
    *                       for more details.
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
    
    void KratosLogUtils::AddCustomFilter(const std::string filter)
    {
        std::string s(filter);
        s.append("::");
        mFilterString->push_back(s);
    }
    

    /**
     * Kenrel function for filters
     * @param input:        input string
     * @param filter:       string to be filtered
     */
    void filterKernel(std::string& input, const std::string& filter)
    {
        while(input.find(filter) != std::string::npos)
        {
            input = input.substr(0,input.find(filter)).append(input.substr(input.find(filter)+filter.size()));
        }
    }
    
    /**
     * Filter kratos generic namespaces (Kratos:: boost:: ...)
     * @param inputString:  string containintg the namespace chain for a function call 
     */   
    const char * KratosLogUtils::Filter(const char * inputString)
    {
        std::string s(inputString);
      
        // Apply default filters
        filterKernel(s,"Kratos::");
        
        // apply custom filters
        for(unsigned int i = 0; i < mFilterString->size(); i++)
        {
            filterKernel(s,(*mFilterString)[i]);
        }
        
        s.append(" ");

        return s.c_str();
    }
}
