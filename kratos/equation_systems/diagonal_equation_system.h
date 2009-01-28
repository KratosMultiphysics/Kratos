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


#if !defined(KRATOS_DIAGONAL_EQUATION_SYSTEM_H_INCLUDED )
#define  KRATOS_DIAGONAL_EQUATION_SYSTEM_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "containers/pointer_vector_set.h"
#include "utilities/indexed_object.h"

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
  template<class TSystemSpace, class TLocalSpace, class TDofType>
  class DiagonalEquationSystem
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of DiagonalEquationSystem
      KRATOS_CLASS_POINTER_DEFINITION(DiagonalEquationSystem);

      typedef typename TSystemSpace::VectorType SystemMatrixType;

      typedef typename TSystemSpace::VectorType SystemVectorType;
  
      typedef typename TSystemSpace::SizeType IndexType;
  
      typedef typename TSystemSpace::SizeType SizeType;
  
      typedef typename TLocalSpace::MatrixType LocalMatrixType;

      typedef typename TLocalSpace::VectorType LocalVectorType;

      typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType;

      typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
  
      typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      DiagonalEquationSystem() :
	mEquationSystemSize(),
	mDirichletSize(),
	mSystemMatrix(), 
	mResults(),
	mRightHandSide() {}

      DiagonalEquationSystem(SizeType NewSystemSize) : 
	mEquationSystemSize(NewSystemSize),
	mDirichletSize(),
	mSystemMatrix(NewSystemSize), 
	mResults(NewSystemSize),
	mRightHandSide(NewSystemSize) {}

     /// Copy constructor.
      DiagonalEquationSystem(const DiagonalEquationSystem& rOther) : 
	mEquationSystemSize(rOther.mEquationSystemSize),
	mDirichletSize(rOther.mDirichletSize),
	mSystemMatrix(rOther.mSystemMatrix), 
	mResults(rOther.mResults),
	mRightHandSide(rOther.mRightHandSide) {}

      /// Destructor.
      virtual ~DiagonalEquationSystem() {}
      

      ///@}
      ///@name Operators 
      ///@{

      /// Assignment operator.
      DiagonalEquationSystem& operator=(const DiagonalEquationSystem& rOther)
      {
	mEquationSystemSize = rOther.mEquationSystemSize;
	mDirichletSize = rOther.mDirichletSize;
	mSystemMatrix = rOther.mSystemMatrix;
	mResults = rOther.mResults;
	mRightHandSide = rOther.mRightHandSide;

	return *this;
      }

 
      
      
      ///@}
      ///@name Operations
      ///@{

      void Initialize(DofsArrayType const& Dofs)
	{
	  AssignEquationIds(Dofs);
	  Resize(mEquationSystemSize);
// 	  TSystemSpace::SetToZero(mSystemMatrix);
// 	  TSystemSpace::SetToZero(mDirichletMatrix);
	  TSystemSpace::SetToZero(mResults);
	  TSystemSpace::SetToZero(mRightHandSide);

	}

      void Initialize()
	{
	  AssignEquationIds();
	  Resize(mEquationSystemSize);
//   	  TSystemSpace::SetToZero(mSystemMatrix);
// 	  TSystemSpace::SetToZero(mDirichletMatrix);
	  TSystemSpace::SetToZero(mResults);
	  TSystemSpace::SetToZero(mRightHandSide);

	}

      void InitializeAllData()
	{
	}

      void InitializeData()
	{
	}

      IndexType AssignEquationIds(DofsArrayType const& Dofs)
	{
	  SetDofs(Dofs);
	  return AssignEquationIds();
	}

      IndexType AssignEquationIds()
	{
	  mEquationSystemSize = AssignDofsEquationIds(mDofs);
	  return mEquationSystemSize;
	}

      static IndexType AssignDofsEquationIds(DofsArrayType & rDofs)
	{
	  SizeType free_id = 0;
	  SizeType fix_id = rDofs.size();
	  
	  for (typename DofsArrayType::iterator dof_iterator = rDofs.begin(); dof_iterator != rDofs.end(); ++dof_iterator)
	    if (dof_iterator->IsFixed())
	      dof_iterator->SetEquationId(--fix_id);
	    else
	      dof_iterator->SetEquationId(free_id++);

	  return fix_id;	  
	}

      void AssembleSystemMatrix(LocalMatrixType const& LocalMatrix, DofsArrayType const& Dofs)
      {
	AssembleMatrix(mSystemMatrix, LocalMatrix, Dofs);
      }
      
      template<class TIndexArrayType>
      void AssembleSystemMatrix(LocalMatrixType const& LocalMatrix, TIndexArrayType const& EquationIds)
      {
	AssembleMatrix(mSystemMatrix, LocalMatrix, EquationIds);
      }

      void AssembleRightHandSide(LocalVectorType const& LocalVector, DofsArrayType const& Dofs)
      {
	AssembleVector(mRightHandSide, LocalVector, Dofs);
      }
      
      template<class TIndexArrayType>
      void AssembleRightHandSide(LocalVectorType const& LocalVector, TIndexArrayType const& EquationIds)
      {
	AssembleVector(mRightHandSide, LocalVector, EquationIds);
      }

      
      void Assemble(LocalMatrixType const& LocalMatrix, LocalVectorType const& LocalVector, DofsArrayType const& Dofs)
      {
	AssembleEquationSystem(mSystemMatrix, mRightHandSide, LocalMatrix, LocalVector, Dofs); 
      }
      
      template<class TIndexArrayType>
      void Assemble(LocalMatrixType const& LocalMatrix, LocalVectorType const& LocalVector, TIndexArrayType const& EquationIds)
      {
	AssembleEquationSystem(mSystemMatrix, mRightHandSide, LocalMatrix, LocalVector, EquationIds); 
      }
      
      static void AssembleMatrix(SystemMatrixType& rSystemMatrix, LocalMatrixType const& LocalMatrix, DofsArrayType const& Dofs)
      {
	KRATOS_TRY_LEVEL_4
	
	SizeType global_size = TSystemSpace::Size(rSystemMatrix);
			
	SizeType local_size = TLocalSpace::Size1(LocalMatrix);
			
	for(SizeType i_local=0; i_local<local_size; i_local++)
	  {
	    SizeType i_global= Dofs[i_local].EquationId();
	    if ( i_global <  global_size)
	      {
		rSystemMatrix[i_global] += LocalMatrix(i_local, i_local);
	      }
	  }
	  KRATOS_CATCH_LEVEL_4("")
      }
      
      template<class TIndexArrayType>
      static void AssembleMatrix(SystemMatrixType& rSystemMatrix, LocalMatrixType const& LocalMatrix, TIndexArrayType const& EquationIds)
      {
	KRATOS_TRY_LEVEL_4
	
	SizeType global_size = TSystemSpace::Size(rSystemMatrix);
			
	SizeType local_size = TLocalSpace::Size1(LocalMatrix);
	for(SizeType i_local=0; i_local<local_size; i_local++)
	  {
	    SizeType i_global = EquationIds[i_local];
	    if ( i_global <  global_size)
	      {
		      rSystemMatrix[i_global] += LocalMatrix(i_local, i_local);
	      }
	  }
	  KRATOS_CATCH_LEVEL_4("")
      }
      
      static void AssembleVector(SystemVectorType& rSystemVector, LocalVectorType const& LocalVector, DofsArrayType const& Dofs)
      {
	KRATOS_TRY_LEVEL_4
	
	SizeType global_size = TSystemSpace::Size(rSystemVector);
			
	SizeType local_size = TLocalSpace::Size(LocalVector);

	for (SizeType i_local=0; i_local<local_size; i_local++)
	  {
	    SizeType i_global= Dofs[i_local].EquationId();
	    if ( i_global <  global_size) //on "free" DOFs
	      {	// ASSEMBLING THE SYSTEM VECTOR
		rSystemVector[i_global] += LocalVector[i_local];
	      }
	  }
	  KRATOS_CATCH_LEVEL_4("")
      }
      
      template<class TIndexArrayType>
      static void AssembleVector(SystemVectorType& rSystemVector, LocalVectorType const& LocalVector, TIndexArrayType const& EquationIds)
      {
	KRATOS_TRY_LEVEL_4
	
	SizeType global_size = TSystemSpace::Size(rSystemVector);
			
	SizeType local_size = TLocalSpace::Size(LocalVector);

	for (SizeType i_local=0; i_local<local_size; i_local++)
	  {
	    SizeType i_global=EquationIds[i_local];
	    if ( i_global < global_size) //on "free" DOFs
	      {	// ASSEMBLING THE SYSTEM VECTOR
		rSystemVector[i_global] += LocalVector[i_local];
	      }
	  }
	  KRATOS_CATCH_LEVEL_4("")
      }

      static void AssembleEquationSystem(SystemMatrixType& rSystemMatrix, 
					 SystemVectorType& rSystemVector,
					 LocalMatrixType const& LocalMatrix,
					 LocalVectorType const& LocalVector,
					 DofsArrayType const& Dofs)
      {
	KRATOS_TRY_LEVEL_4
	
	SizeType global_size = TSystemSpace::Size(rSystemMatrix);
			
	SizeType local_size = TLocalSpace::Size1(LocalMatrix);
			
 	for(SizeType i_local=0; i_local<local_size; i_local++)
 	  {
 	    SizeType i_global = Dofs[i_local].EquationId();
 	    if ( i_global <  global_size)
 	      {
		rSystemMatrix[i_global] += LocalMatrix(i_local, i_local);
 		rSystemVector[i_global] += LocalVector[i_local];
 	      }
 	  }
	  KRATOS_CATCH_LEVEL_4("")
      }
      
      template<class TIndexArrayType>
      static void AssembleEquationSystem(SystemMatrixType& rSystemMatrix, 
					 SystemVectorType& rSystemVector,
					 LocalMatrixType const& LocalMatrix, 
					 LocalVectorType const& LocalVector,
					 TIndexArrayType const& EquationIds)
      {
	KRATOS_TRY_LEVEL_4
	
	SizeType global_size = TSystemSpace::Size(rSystemMatrix);
			
	SizeType local_size = TLocalSpace::Size1(LocalMatrix);
			
 	for(SizeType i_local=0; i_local<local_size; i_local++)
 	  {
 	    SizeType i_global=EquationIds[i_local];
 	    if ( i_global <  global_size)
 	      {
		rSystemMatrix[i_global] += LocalMatrix(i_local, i_local);
 		rSystemVector[i_global] += LocalVector[i_local];
 	      }
 	  }
	  KRATOS_CATCH_LEVEL_4("")
      }
      
      void Resize(SizeType NewSize)
	{
	  mEquationSystemSize = NewSize;
	  mDirichletSize = (mDofs.size() > mEquationSystemSize) ? (mDofs.size() - mEquationSystemSize) : 0;
	  TSystemSpace::Resize(mSystemMatrix, mEquationSystemSize, mEquationSystemSize);
	  TSystemSpace::Resize(mDirichletMatrix, mEquationSystemSize, mDirichletSize);
	  TSystemSpace::Resize(mResults, mEquationSystemSize);
	  TSystemSpace::Resize(mRightHandSide, mEquationSystemSize);
	  std::cout << "Equation system size : " << mEquationSystemSize << std::endl;
	  std::cout << "Dirichlet system size : " << mDirichletSize << std::endl;
	}

      void Clear()
	{
	  mEquationSystemSize = 0;
	  mDirichletSize = 0;
	  mNonZeros =0;
	  mDirichletNonZeros =0;
	  
	  TSystemSpace::Clear(mSystemMatrix);
	  TSystemSpace::Clear(mDirichletMatrix);
	  TSystemSpace::Clear(mResults);
	  TSystemSpace::Clear(mRightHandSide);
	}

      void ClearAllData()
	{
	}

      void ClearData()
	{
	}

      ///@}
      ///@name Access
      ///@{

      SizeType Size()
	{
	  return mEquationSystemSize;
	}

      SizeType DirichletSize()
	{
	  return mDirichletSize;
	}

      SystemMatrixType& GetSystemMatrix()
      {
	return mSystemMatrix;
      }
      
      const SystemMatrixType& GetSystemMatrix() const
      {
	return mSystemMatrix;
      }

      SystemVectorType& GetResults()
      {
	return mResults;
      }

      const SystemVectorType& GetResults() const
      {
	return mResults;
      }
      
      SystemVectorType& GetRightHandSide()
      {
	return mRightHandSide;
      }

      const SystemVectorType& GetRightHandSide() const
      {
	return mRightHandSide;
      }

      SystemMatrixType& GetDirichletMatrix()
      {
	return mDirichletMatrix;
      }
      
      const SystemMatrixType& GetDirichletMatrix() const
      {
	return mDirichletMatrix;
      }

      void SetSystemMatrix(const SystemMatrixType& NewSystemMatrix)
      {
  	TSystemSpace::Copy(NewSystemMatrix, mSystemMatrix);
      }
      
      void SetResults(const SystemVectorType& NewResults)
      {
	TSystemSpace::Copy(NewResults, mResults);
      }

      void SetRightHandSide(const SystemVectorType& NewRightHandSide)
      {
	return TSystemSpace::Copy(NewRightHandSide, mRightHandSide);
      }

      DofsArrayType const& GetDofs() const
	{
	  return mDofs;
	}

      DofsArrayType& GetDofs()
	{
	  return mDofs;
	}

      void SetDofs(DofsArrayType const& ThisDofs)
	{
	  mDofs = ThisDofs;
	}

      void AddDofs(DofsArrayType const& ThisDofs)
	{
	  for(typename DofsArrayType::ptr_iterator i = ThisDofs.ptr_begin() ; i != ThisDofs.ptr_end() ; ++i)
	    mDofs.push_back(*i);
	}

      DofIterator DofsBegin() {return mDofs.begin();}

      DofConstantIterator DofsBegin() const {return mDofs.begin();}

      DofIterator DofsEnd() {return mDofs.end();}

      DofConstantIterator DofsEnd() const {return mDofs.end();}

      

      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	return "Equation system";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	rOStream << "Equation system";
      }


      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	rOStream << "    Equation system size :" << mEquationSystemSize << std::endl;
	rOStream << "    Dirichlet size       :" << mDirichletSize << std::endl;
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

      SizeType mEquationSystemSize;

      SizeType mDirichletSize;

      SizeType mNonZeros;

      SizeType mDirichletNonZeros;

      SystemMatrixType mSystemMatrix;

      SystemMatrixType mDirichletMatrix;

      SystemVectorType mResults;

      SystemVectorType mRightHandSide;

      DofsArrayType mDofs;

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
        
    }; // Class DiagonalEquationSystem 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSystemSpace, class TLocalSpace, class TDofType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    DiagonalEquationSystem<TSystemSpace, TLocalSpace, TDofType>& rThis);

  /// output stream function
  template<class TSystemSpace, class TLocalSpace, class TDofType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const DiagonalEquationSystem<TSystemSpace, TLocalSpace, TDofType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_DIAGONAL_EQUATION_SYSTEM_H_INCLUDED  defined 


