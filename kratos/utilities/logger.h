//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Pooyan Dadvand 
//
	           

#if !defined(KRATOS_LOGGER_H_INCLUDED )
#define  KRATOS_LOGGER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_exception.h"



namespace Kratos
{
  ///@addtogroup Kratos
  ///@{

  ///@name Kratos Macros
  ///@{ 

#if defined(KRATOS_HERE)
#undef KRATOS_HERE
#endif

#define KRATOS_HERE Logger::CleanFunctionName(BOOST_CURRENT_FUNCTION , __FILE__ , __LINE__)

#define KRATOS_ERROR throw KratosException("Error: ", KRATOS_HERE)   

#define KRATOS_TIME_STAMP "[" << KratosLogUtils::GetInstance().CurrentDateTime() << "]"
#define KRATOS_PROCESS_ID "[PID=" << getpid() << "]"

#define KRATOS_ERROR_STAMP_DETAIL(file,line,function) "ERROR:" << KRATOS_HERE
#define KRATOS_ERROR_STAMP KRATOS_ERROR_STAMP_DETAIL(__FILE__,__LINE__,KRATOS_LOG_FILTER(BOOST_CURRENT_FUNCTION))

//#define KRATOS_ERROR throw KratosException();  
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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) Logger
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Logger
      KRATOS_CLASS_POINTER_DEFINITION(Logger);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Logger();

      /// Destructor.
      virtual ~Logger();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Static Methods
      ///@{
      
	  static std::string CleanFunctionName(const std::string& FunctionName, const std::string& FileName, int LineNumber);

	  static std::string Filter(const std::string& ThisString);
      
	  static std::string ReplaceAll(std::string& ThisString, const std::string& FromString, const std::string& ToString);
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const;
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;
      
            
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
      
      /// Assignment operator.
      Logger& operator=(Logger const& rOther);

      /// Copy constructor.
      Logger(Logger const& rOther);

        
      ///@}    
        
    }; // Class Logger 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Logger& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Logger& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LOGGER_H_INCLUDED  defined 


