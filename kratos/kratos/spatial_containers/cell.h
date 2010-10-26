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
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2010-10-08 16:07:33 $
//   Revision:            $Revision: 1.0$
//
//


#if !defined(KRATOS_CELL_H_INCLUDED)
#define  KRATOS_CELL_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>

// External includes 
//#include "boost/smart_ptr.hpp"


// Project includes
//#include "includes/define.h"
//#include "includes/ublas_interface.h"
//#include "includes/properties.h"


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
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<
  class TPointerType,
  class TContainerType = typename std::vector<TPointerType>, 
  class TIteratorType  = typename TContainerType::iterator
  > 
  class Cell
    {
    public:
      ///@name Type Definitions
      ///@{
      
      //typedef std::vector<TPointerType> ContainerType;
      typedef TContainerType ContainerType;
      
      typedef TIteratorType IteratorType; 
      
      typedef std::size_t  SizeType;
      
      //typedef std::vector< IndexType > CellIndex;
            
      //typedef Cell<TDimension, TPointerType>  CellType; 
      
      /// Pointer definition of Cell
      KRATOS_CLASS_POINTER_DEFINITION(Cell);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Cell()
      {
      }

      /// Destructor.
      virtual ~Cell(){}
      
      
      void Add(TPointerType& ThisObject)
      {
	mObjects.push_back(ThisObject);
      }
      
      void Clear()
      {
	mObjects.clear();
      }
      
      /// WARNING = Solo para ser allocatado
      void AllocateCell(const std::size_t size)
      {
	mObjects.reserve(size);
      }
      
      SizeType GetContainerSize()
      {
	return mObjects.size();  
      }
            
       ContainerType GetContainer() 
       {
	 return mObjects;
       }
       
       IteratorType Begin() 
       {
	 return mObjects.begin();
       }
       
       IteratorType End()  
       {
	 return mObjects.end();
       }
       
       IteratorType Begin() const   
        { 	 
	  return mObjects.begin();
	}
        
       IteratorType End() const  
       {
	 return mObjects.end();
       }
       
       
       
       
      ///@}
      ///@name Operators 
      ///@{
      
      //Cell& operator=(Cell const& rOther){}
      
            

      /// Copy constructor.
      //Cell(Cell const& rOther){}
      
      
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
	return "Cell Class "; 
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
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
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        

      ContainerType mObjects;
        
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
        
    }; // Class Cell 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template< class TPointerType,class TContainerType, class TIteratorType> 
  inline std::istream& operator >> (std::istream& rIStream, 
				    Cell<TPointerType>& rThis){ return rIStream;}

  /// output stream function
  template< class TPointerType, class TContainerType, class TIteratorType> 
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Cell<TPointerType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


