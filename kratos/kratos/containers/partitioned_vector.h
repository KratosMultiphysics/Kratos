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
//	 Project Name:		  Kratos	   
//	 Last Modified by:	  $Author: pooyan $
//	 Date:				  $Date: 2008-12-01 10:09:01 $
//	 Revision:			  $Revision: 1.1 $
//
//


#if	!defined(KRATOS_PARTITIONED_VECTOR_H_INCLUDED	)
#define	 KRATOS_PARTITIONED_VECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>	


// External	includes 

// Project includes
#include "includes/define.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{	
  
  ///@}	
  ///@name Type	Definitions
  ///@{	
  
  ///@}	
  ///@name	Enum's
  ///@{
	  
  ///@}
  ///@name	Functions 
  ///@{
	  
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short	class definition.
  /** Detail class definition.
  */
  template<class TVectorType>
  class	partitioned_vector	
	{
	public:
	  ///@name Type	Definitions
	  ///@{
	  
	  /// Pointer definition of	partitioned_vector
	  KRATOS_CLASS_POINTER_DEFINITION(partitioned_vector);
	  
	  typedef std::size_t size_type;
	  typedef TVectorType VectorType;
	  typedef typename VectorType::value_type value_type;
	  typedef value_type& reference;
	  typedef value_type const& const_reference;
	  typedef std::vector<VectorType> VectorsContainerType;
	  typedef std::vector<size_type> PartitionsIndicesVectorType;

	  ///@}
	  ///@name Life	Cycle 
	  ///@{	
	  
	  /// Default constructor.
	  partitioned_vector () : mVectors(),  mPartitionsIndices(1,0) // one partition with size zero
	  {
	  }

	  template<class TSizesVectorType>
	    partitioned_vector (TSizesVectorType const& PartitionSizes) 
	    : mVectors(PartitionSizes.size()),  mPartitionsIndices(PartitionSizes.size())
	    {
	      resize(PartitionSizes);
	    }



	  ///@}
	  ///@name Operators 
	  ///@{

	  // Element access
	  const_reference	operator ()	(size_type i) const	
	  {
	    return GetValue(i);
	  }
	  
	  reference operator () (size_type i)	
	  {
	    return GetValue(i);
	  }

	  const_reference	operator []	(size_type i) const	
	  {
	    return GetValue(i);
	  }

	  reference operator [] (size_type i)	
	  {
	    return GetValue(i);
	  }

	  // Assignment
	  partitioned_vector &operator = (const partitioned_vector &v) 
	  {
	    mVectors= v.mVectors;
	    mPartitionsIndices= v.mPartitionsIndices;
	    return *this;
	  }
	  
	  
	  ///@}
	  ///@name Operations
	  ///@{
	  
	  // Resizing
	  template<class TSizesVectorType>
	  void resize(TSizesVectorType const& PartitionSizes) 
	  {
	    mVectors.resize(PartitionSizes.size());
	    mPartitionsIndices.resize(PartitionSizes.size() + 1);

	    mPartitionsIndices[0] = 0;

	    for(size_type i = 0 ; i < PartitionSizes.size() ; i++)
	      {
		mVectors[i].resize(PartitionSizes[i]);
		mPartitionsIndices[i+1] = PartitionSizes[i] + mPartitionsIndices[i];
	      }
	  }

	  // Swapping
	  void swap (partitioned_vector &v) 
	  {
	    if (this !=	&v)	
	    {
	      mVectors.swap (v.mVectors);
	      mPartitionsIndices.swap (v.mPartitionsIndices);
	    }
	  }



	  ///@}
	  ///@name Access
	  ///@{	
	  
	public:

	  size_type size () const	
	  {
	    return mPartitionsIndices.back();
	  }

	  const VectorsContainerType &data () const	
	  {
	    return mVectors;
	  }

	  VectorsContainerType &data () 
	  {
	    return mVectors;
	  }

	  
	  const VectorType &partition (size_type i) const	
	  {
	    return mVectors[i];
	  }

	  VectorType &partition (size_type i) 
	  {
	    return mVectors[i];
	  }

	  const VectorType &partition_index (size_type i) const	
	  {
	    return mPartitionsIndices[i];
	  }

	  

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

	protected:
	  ///@name Protected static	Member Variables 
	  ///@{	
		
		
	  ///@}	
	  ///@name Protected member	Variables 
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
		
	  VectorsContainerType mVectors;
	  PartitionsIndicesVectorType mPartitionsIndices;
		
	  ///@}	
	  ///@name Private Operators
	  ///@{	
		
	  reference GetValue(size_type i)	
	  {
	    for(size_type partition_index = 0 ; partition_index < mVectors.size() ; partition_index++)
	      if(i < mPartitionsIndices[partition_index+1])
		return mVectors[partition_index][i -  mPartitionsIndices[partition_index]];
	  }

		
	  ///@}	
	  ///@name Private Operations
	  ///@{	
		
		
	  ///@}	
	  ///@name Private	Access 
	  ///@{	
		
		
	  ///@}	   
	  ///@name Private Inquiry 
	  ///@{	
		
		
	  ///@}	   
	  ///@name Un accessible methods 
	  ///@{	
	  
		
	  ///@}	   
		
	}; // Class	partitioned_vector	

  ///@}	
  
  ///@name Type	Definitions		  
  ///@{	
  
  
  ///@}	
  ///@name Input and output	
  ///@{	
		
  ///@}	
  
  
}  // namespace	Kratos.

#endif // KRATOS_PARTITIONED_VECTOR_H_INCLUDED  defined 
