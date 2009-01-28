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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2007-06-27 11:04:33 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_VARIABLES_LIST_H_INCLUDED )
#define  KRATOS_VARIABLES_LIST_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <vector> 


// External includes 
#include <boost/iterator/indirect_iterator.hpp>


// Project includes
#include "includes/define.h"
//#include "includes/kratos_components.h"
#include "containers/variable.h"


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
  class VariablesList
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of VariablesList
      KRATOS_CLASS_POINTER_DEFINITION(VariablesList);



      typedef std::size_t SizeType;

      typedef std::size_t IndexType;

	  typedef double BlockType;

      typedef std::vector<IndexType> PositionsContainerType;

      typedef std::vector<const VariableData*> VariablesContainerType;

	  typedef VariableData data_type;
	  typedef const VariableData* value_type;
      typedef const VariableData* const_pointer;
	  typedef VariableData const& const_reference;
  
	  typedef boost::indirect_iterator<VariablesContainerType::const_iterator>          const_iterator;                
	  typedef boost::indirect_iterator<VariablesContainerType::const_reverse_iterator>  const_reverse_iterator;              

      typedef VariablesContainerType::size_type size_type;
	  typedef VariablesContainerType::const_iterator ptr_const_iterator;
	  typedef VariablesContainerType::const_reverse_iterator ptr_const_reverse_iterator;
	  typedef VariablesContainerType::difference_type difference_type;
	  //typedef typename VariablesContainerType::iterator iterator;
	  //typedef typename VariablesContainerType::const_iterator const_iterator;
	  //typedef typename VariablesContainerType::reverse_iterator reverse_iterator;
	  //typedef typename VariablesContainerType::const_reverse_iterator const_reverse_iterator;


      ///@}
      ///@name Life Cycle 
	  ///@{ 

	  /// Default constructor.
	  VariablesList() : mDataSize(0), mPositions(), mVariables()
	  {
	  }

      template <class TInputIteratorType>
      VariablesList(TInputIteratorType First, TInputIteratorType Last) 
      {
		for(; First != Last ; First++)
			push_back(*First);
      }

	  /// Copy constructor.
	  VariablesList(VariablesList const& rOther) : mDataSize(rOther.mDataSize)
		  , mPositions(rOther.mPositions)
		  , mVariables(rOther.mVariables) {}      

	  /// Destructor.
      virtual ~VariablesList()
      {
      }
      

      ///@}
      ///@name Operators 
      ///@{

      /// Assignment operator.
      VariablesList& operator=(VariablesList const& rOther)
	  {
		  mDataSize = rOther.mDataSize;
		  mPositions = rOther.mPositions;
		  mVariables = rOther.mVariables;

		return *this;
	  }

      IndexType operator()(IndexType VariableKey) const
      {
		return mPositions[VariableKey];
      }

      template<class TDataType>
      IndexType operator()(Variable<TDataType> const& ThisVariable) const
      {
		return mPositions[ThisVariable.Key()];
      }

      const VariableData* operator[](IndexType Index) const
      {
		return mVariables[Index];
      }

    
    bool operator==( const VariablesList& r ) const // nothrow
    { 
      if( size() != r.size() )
      return false;
      else
      return std::equal(mPositions.begin(), mPositions.end(), r.mPositions.begin()) &&
			 std::equal(mVariables.begin(), mVariables.end(), r.mVariables.begin());
    }
            
      
      ///@}
      ///@name Operations
      ///@{
      
    const_iterator             begin() const      { return const_iterator( mVariables.begin() ); }
    const_iterator             end() const        { return const_iterator( mVariables.end() ); }
    const_reverse_iterator     rbegin() const     { return const_reverse_iterator( mVariables.rbegin() ); } 
    const_reverse_iterator     rend() const       { return const_reverse_iterator( mVariables.rend() ); } 
    ptr_const_iterator         ptr_begin() const  { return mVariables.begin(); }
    ptr_const_iterator         ptr_end() const    { return mVariables.end(); }
    ptr_const_reverse_iterator ptr_rbegin() const { return mVariables.rbegin(); }
    ptr_const_reverse_iterator ptr_rend() const   { return mVariables.rend(); }

    const_reference  front() const /* nothrow */ { assert( !IsEmpty() ); return *mVariables.front(); }
    const_reference  back() const  /* nothrow */ { assert( !IsEmpty() ); return *mVariables.back(); }
    
    size_type size() const {return mVariables.size();}
    
    size_type max_size() const {return mVariables.max_size();}
    
     void swap(VariablesList& rOther) 
	 {
		 SizeType temp = mDataSize;
		 mDataSize = rOther.mDataSize;
		 rOther.mDataSize = temp;

		 mVariables.swap(rOther.mVariables);
		 mPositions.swap(rOther.mPositions);
	 }
    
  template<class TOtherDataType>
    void push_back(TOtherDataType const& x)
    {
      Add(x); 
    }
        
  //template<class TOtherDataType>
  //  iterator insert(iterator Position, const TOtherDataType& rData)
  //  {
  //    return iterator(mVariables.insert(Position, TPointerType(new TOtherDataType(rData))));
  //  }
  //  
  //  iterator insert(iterator Position, const TPointerType pData)
  //  {
  //    return iterator(mVariables.insert(Position, pData));
  //  }

  //  template <class InputIterator>
  //  void insert(InputIterator First, InputIterator Last)
  //  {
  //    for(;First != Last; ++First)
  //      insert(*First);
  //  }
    
                       
    //iterator erase(iterator pos) {return iterator(mVariables.erase(pos.base()));}
    //
    //iterator erase( iterator first, iterator last )
    //{
    //  return iterator( mVariables.erase( first.base(), last.base() ) );
    //}

    void clear() 
	{
		mDataSize = 0;
		mVariables.clear();
		mPositions.clear();
	}

	//void reserve(int dim){mVariables.reserve(dim);}
	
	//int capacity()
	//{
	//	return mVariables.capacity();
	//}    
    
     template<class TDataType>
	void Add(Variable<TDataType> const& ThisVariable)
	{
	  if(ThisVariable.Key()== 0)
		  KRATOS_ERROR(std::logic_error, 
		  "Adding uninitialize variable to this variable list. Check if all variables are registered before kernel initialization","");

	  if(Has(ThisVariable))
		  return;

	  if(mPositions.size() <= ThisVariable.Key())
	    mPositions.resize(ThisVariable.Key()+1/*KratosComponents<VariableData>::Size()*/, static_cast<IndexType>(-1));

	  mPositions[ThisVariable.Key()] = mDataSize;
	  mVariables.push_back(&ThisVariable);
	  const SizeType block_size = sizeof(BlockType);
	  mDataSize += static_cast<SizeType>(((block_size - 1) + sizeof(TDataType)) / block_size);
	}

      IndexType Index(IndexType VariableKey) const
      {
	return mPositions[VariableKey];
      }

      template<class TDataType>
      IndexType Index(Variable<TDataType> const& ThisVariable) const
      {
	return mPositions[ThisVariable.Key()];
      }

      IndexType Index(const VariableData* pThisVariable) const
      {
	return mPositions[pThisVariable->Key()];
      }

      
      ///@}
      ///@name Access
      ///@{

      SizeType DataSize() const
	{
	  return mDataSize;
	}

      
      VariablesContainerType const& Variables()
	{
	  return mVariables;
	}
      
      ///@}
      ///@name Inquiry
      ///@{
      
      template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
	{
	  if(mPositions.empty())
	    return false;

	  if(rThisVariable.Key()== 0)
	    return false;

	  if(mPositions.size() <= rThisVariable.Key())
		  return false;

	  return (Index(rThisVariable) < mDataSize);
	}

      bool IsEmpty() const
	{
	  return mVariables.empty();
	}
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	return "variables list";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
		  rOStream << " with " << size() << " variables";
		  rOStream << " (size : " << mDataSize << " blocks of " << sizeof(BlockType) << " bytes) "<< std::endl;
	for(IndexType i = 0 ; i < mVariables.size() ; ++i)
	  rOStream << "    " << mVariables[i]->Name() << " \t-> " << mPositions[mVariables[i]->Key()] << std::endl;
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

      SizeType mDataSize;
      
      PositionsContainerType mPositions;

      VariablesContainerType mVariables;        
        
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
        
    }; // Class VariablesList 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    VariablesList& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const VariablesList& rThis)
    {
      rThis.PrintInfo(rOStream);
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VARIABLES_LIST_H_INCLUDED  defined 


