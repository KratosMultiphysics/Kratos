//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi 
//
	           
#if !defined(KRATOS_GLOBAL_POINTER_UNORDERED_MAP_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_UNORDERED_MAP_H_INCLUDED

// System includes
#include <vector>
#include <iostream> 
#include <utility>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/global_pointer.h"
#include "includes/serializer.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

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
  template< class TDataType, class TValueType >
  class GlobalPointersUnorderedMap : public std::unordered_map< GlobalPointer<TDataType>, 
                                                                TValueType,  
                                                                GlobalPointerHasher<TDataType>, 
                                                                GlobalPointerComparor<TDataType> 
                                                                >
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of GlobalPointersUnorderedMap
      KRATOS_CLASS_POINTER_DEFINITION(GlobalPointersUnorderedMap);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      GlobalPointersUnorderedMap(){}

      /// Destructor.
      virtual ~GlobalPointersUnorderedMap(){}     

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
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
      virtual std::string Info() const
      {
	    std::stringstream buffer;
        buffer << "GlobalPointersUnorderedMap" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GlobalPointersUnorderedMap";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
            
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
        friend class Serializer;

        void save(Serializer& rSerializer) const
        {
            rSerializer.save("Size", this->size());
            for(const auto& item : (*this))
            {
                rSerializer.save("D", item);
            }
        }

        void load(Serializer& rSerializer)
        {
            std::size_t size;
            rSerializer.load("Size", size);
            this->reserve(size);

            for(std::size_t i = 0; i<size; ++i)
            {
                std::pair< GlobalPointer<TDataType>, TValueType>  tmp(nullptr, 0);
                rSerializer.load("D", tmp);
                this->insert(tmp);
            }
        }
        
        
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
      GlobalPointersUnorderedMap& operator=(GlobalPointersUnorderedMap const& rOther){}

      /// Copy constructor.
      GlobalPointersUnorderedMap(GlobalPointersUnorderedMap const& rOther){}

        
      ///@}    
        
    }; // Class GlobalPointersUnorderedMap 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template< class TDataType, class TValueType >
  inline std::istream& operator >> (std::istream& rIStream, 
				    GlobalPointersUnorderedMap<TDataType,TValueType>& rThis)
    {
        return rIStream;
    }

  /// output stream function
  template< class TDataType, class TValueType >
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const GlobalPointersUnorderedMap<TDataType,TValueType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_GLOBAL_POINTER_UNORDERED_MAP_H_INCLUDED  defined 


