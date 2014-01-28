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
#define KRATOS_LOG_ERROR_INFO(file,line,function,severity) "[" << #severity << "]:" << (file) << ":" << (line) << ":"

#define KRATOS_ERROR_LOG        KRATOS_LOG(KRATOS_SEVERITY_ERROR) << KRATOS_LOG_ERROR_INFO(__FILE__,__LINE__,BOOST_CURRENT_FUNCTION,KRATOS_SEVERITY_ERROR)
#define KRATOS_ERROR_LOG_N      KRATOS_LOG_N(KRATOS_SEVERITY_ERROR)	
#define KRATOS_ERROR_LOG_IF     KRATOS_LOG_IF(KRATOS_SEVERITY_ERROR)	
#define KRATOS_ERROR_LOG_CHECK	KRATOS_LOG_CHECK(KRATOS_SEVERITY_ERROR)	

#ifdef  KRATOS_DEBUG
#define KRATOS_DEBUG_LOG        KRATOS_LOG(KRATOS_DEBUG)
#define KRATOS_TRACE_LOG        KRATOS_LOG(KRATOS_TRACE)
#else
#define KRATOS_DEBUG_LOG    
#define KRATOS_TRACE_LOG    
#endif

#define KRATOS_LOG(SEVERITY)                                                           \
  KratosLog<SEVERITY>::GetInstance()    \

# define KRATOS_LOG_OCCURRENCES(line) kratos_log_loop_counter ## line

// Each N
#define KRATOS_LOG_N(KRATOS_SEVERITY) \
  {static int KRATOS_LOG_OCCURRENCES(__LINE__) = 0, KRATOS_LOG_OCCURRENCES_MOD_N = 0; \
  ++LOG_OCCURRENCES; \
  if (++KRATOS_LOG_OCCURRENCES > 10) KRATOS_LOG_OCCURRENCES -= n; \
  if (KRATOS_LOG_OCCURRENCES == 1) \
    KRATOS_LOG(KRATOS_SEVERITY) << KRATOS_LOG_ERROR_INFO(__FILE__,__LINE__,BOOST_CURRENT_FUNCTION,KRATOS_SEVERITY_ERROR)}

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
  
enum KRATOS_SEVERITY {
  KRATOS_SEVERITY_ERROR,
  KRATOS_SEVERITY_WAR,
  KRATOS_SEVERITY_INF,
  KRATOS_SEVERITY_DET,
  KRATOS_SEVERITY_DEB,
  KRATOS_SEVERITY_TRC,
  KRATOS_LOG_SOLVER
};
 
template<KRATOS_SEVERITY S>  
class KratosLog : public std::ostream
{
  private:
    
    KratosLog<S>();
    
    static bool         mInstanceFalg;
    static KratosLog *  mpInstance;
    
  public:
    
    static KratosLog<S>& GetInstance();
    
    void SetSeverityLevel(int severity);
    void SetOutput(std::string const& LogFileName);
    
    int getSeverity();
    std::ostream& getOutput();
    
    bool            mPrintLogOnScreen;
    int             mSeverityLevel;
    std::ofstream * mOutputFile;
    
    ~KratosLog<S>();
    
};  


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

inline KratosLog<KRATOS_SEVERITY_ERROR> & operator << (KratosLog<KRATOS_SEVERITY_ERROR> &out, std::ostream& (*fn)(std::ostream&))
{ 
  if(out.mPrintLogOnScreen)
    fn(std::cout);
     
  fn((*out.mOutputFile));
        
  return out;
};
  
template<class T>
inline KratosLog<KRATOS_SEVERITY_ERROR> & operator << (KratosLog<KRATOS_SEVERITY_ERROR> &out, const T& data)
{
  if(out.mPrintLogOnScreen)
  {
      std::cout << data;
  }

  (*out.mOutputFile) << data;
  (*out.mOutputFile).flush();

//   if(out.mSeverityLevel == getKratosSeverityError()?)
//   {
//       (*out.mOutputFile).flush();
//       abort();
//   }

  return out;  
};

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_LOG_H_INCLUDED  defined 


