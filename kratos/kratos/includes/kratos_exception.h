//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//
	           

#if !defined(KRATOS_KRATOS_EXCEPTION_H_INCLUDED )
#define  KRATOS_KRATOS_EXCEPTION_H_INCLUDED



// System includes
#include <stdexcept>
#include <string>
#include <iostream> 


// External includes 


// Project includes


namespace Kratos
{
  ///@addtogroup Kratos
  ///@{

  ///@name Kratos Macros
  ///@{ 

  ///@} 
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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) KratosException : public std::exception
    {
    public:
      ///@name Type Definitions
      ///@{
      
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  KratosException();

	  KratosException(const std::string& rWhat );

      KratosException(const std::string& rWhat, const std::string& rWhere);
	  
	  /// Copy constructor.
      KratosException(KratosException const& rOther);

      /// Destructor.
      virtual ~KratosException() throw(); //noexcept; // noexcept(true);
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// string stream function
      template<class StreamValueType>
      KratosException& operator << (StreamValueType const& rValue)
      {
          std::stringstream buffer;
          buffer << rValue;

          append_message(buffer.str());

          return *this;
      }

      /// Manipulator stream function
    //  template<class StreamValueType>
      KratosException& operator << (std::ostream& (*pf)(std::ostream&));
      /// char stream function
      KratosException& operator << (const char * rString);

      ///@}
      ///@name Operations
      ///@{

	  void append_message(std::string const& rWhat);

	  void append_where(std::string const& rWhere);
      
      
      ///@}
      ///@name Access
      ///@{ 

	  /// The overide of the base class what method
	  /** This method returns the entire message with where information
	  */

      const char* what() const throw () /*noexcept*/; // Todo: I should change this after switching to c++11. Pooyan.

	  const std::string& message() const;

	  const std::string& where() const;

      
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

		std::string mWhat;
		std::string mWhere;
		std::string mMessage;
        
        
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
      KratosException& operator=(KratosException const& rOther);

        
      ///@}    
        
    }; // Class KratosException 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  std::istream& operator >> (std::istream& rIStream, 
  				    KratosException& rThis);

  /// output stream function
  //std::ostream& operator << (std::ostream& rOStream,
  //	  const KratosException& rThis);

//  /// string stream function
//  template<class StreamValueType>
//  KratosException operator << (KratosException rThis, StreamValueType const& rValue)
//  {
//	  std::stringstream buffer;
//	  buffer << rValue;

//	  rThis.append_message(buffer.str());

//	  return rThis;
//  }

//  /// Manipulator stream function
////  template<class StreamValueType>
//  KRATOS_API(KRATOS_CORE) KratosException operator << (KratosException rThis, std::ostream& (*pf)(std::ostream&));
//  /// char stream function
//  KRATOS_API(KRATOS_CORE) KratosException operator << (KratosException rThis, const char * rString);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_KRATOS_EXCEPTION_H_INCLUDED  defined 


