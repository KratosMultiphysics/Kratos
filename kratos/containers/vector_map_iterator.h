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
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_VECTOR_MAP_ITERATOR_H_INCLUDED )
#define  KRATOS_VECTOR_MAP_ITERATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 
#include <boost/iterator/iterator_adaptor.hpp>


// Project includes
#include "includes/define.h"


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
	template<class TIteratorType, class TDataType>
  class VectorMapIterator
	  : public boost::iterator_adaptor<VectorMapIterator<TIteratorType, TDataType>,
	  TIteratorType, TDataType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of VectorMapIterator
      KRATOS_CLASS_POINTER_DEFINITION(VectorMapIterator);

	  typedef boost::iterator_adaptor<VectorMapIterator,
	  TIteratorType,  TDataType> BaseType;

	  typedef typename TIteratorType::value_type::first_type key_type;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  VectorMapIterator(){}

	  VectorMapIterator(TIteratorType NewIterator) :BaseType(NewIterator) {}

	  VectorMapIterator(VectorMapIterator const & NewIterator) :BaseType(NewIterator.base()) {}

	  //template<class TOtherIteratorType>
	  //VectorMapIterator(VectorMapIterator<TIteratorType> const & NewIterator) :BaseType(NewIterator.base()) {}


      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

	key_type key()
	{
		return this->base()->first;
	}
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    private:
   friend class boost::iterator_core_access;
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
        
      typename BaseType::reference dereference() const
      {
# if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x551))
          return const_cast<BaseType::reference>(*this->base()->second);
# else
          return this->base()->second;
# endif 
      }
        
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
        
    }; // Class VectorMapIterator 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
   ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VECTOR_MAP_ITERATOR_H_INCLUDED  defined 


