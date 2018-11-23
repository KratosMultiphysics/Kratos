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
	           
#if !defined(KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED

// System includes
#include <vector>
#include <iostream> 


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
  template< class TDataType >
  class GlobalPointersVector : public std::vector< GlobalPointer<TDataType> >
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of GlobalPointersVector
      KRATOS_CLASS_POINTER_DEFINITION(GlobalPointersVector);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      GlobalPointersVector(){}

      /// Destructor.
      virtual ~GlobalPointersVector(){}

      template < class TContainerType > 
      void FillFromContainer( TContainerType& container)
      {
          this->reserve(container.size());
          for(auto& item : container)
          {
              this->push_back(GlobalPointer<TDataType>(&item));
          }
      }
      

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
        buffer << "GlobalPointersVector" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GlobalPointersVector";}

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
            std::size_t pointer_size = sizeof(GlobalPointer<TDataType> );

            std::string data;
            data.resize( this->size() * pointer_size);
            for(std::size_t i=0; i<this->size(); ++i)
            {
                (*this)[i].save(&data[0]+i*pointer_size);
            }

            rSerializer.save("Size", this->size());
            rSerializer.save("Data", data);

            // the following version works but it should be less optimized than the version above
            // rSerializer.save("Size", this->size());
            // for(const auto& item : *this)
            //     rSerializer.save("Gp", item);
        }

        void load(Serializer& rSerializer)
        {
            std::size_t pointer_size = sizeof(GlobalPointer<TDataType> );

            std::size_t size;
            rSerializer.load("Size", size);
            this->reserve(size);

            std::string tmp;
            rSerializer.load("Data", tmp);

            for(std::size_t i = 0; i<size; ++i)
            {
                GlobalPointer<TDataType> p(nullptr);
                p.load(&tmp[0]+i*pointer_size);
                this->push_back(p);
           }

            // the following version works but it should be less optimized than the version above
            // std::size_t size;
            // rSerializer.load("Size", size);
            // this->reserve(size);
            // for(std::size_t i = 0; i<size; ++i)
            // {
            //     GlobalPointer<TDataType> p(nullptr);
            //     rSerializer.load("Gp", p);
            //     this->push_back(p);
            // }
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
      GlobalPointersVector& operator=(GlobalPointersVector const& rOther){}

      /// Copy constructor.
      GlobalPointersVector(GlobalPointersVector const& rOther){}

        
      ///@}    
        
    }; // Class GlobalPointersVector 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template< class TDataType >
  inline std::istream& operator >> (std::istream& rIStream, 
				    GlobalPointersVector<TDataType>& rThis){}

  /// output stream function
  template< class TDataType >
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const GlobalPointersVector<TDataType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED  defined 


