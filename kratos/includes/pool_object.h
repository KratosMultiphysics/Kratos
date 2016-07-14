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

#if !defined(KRATOS_POOL_OBJECT_H_INCLUDED )
#define  KRATOS_POOL_OBJECT_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/memory_pool.h"


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

  /// This redefines the new and delete operators of derived class to be allocated in pool
  /** The PoolObject is the base class for classes to be allocated in the pool
  */
  class PoolObject
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of PoolObject
      KRATOS_CLASS_POINTER_DEFINITION(PoolObject);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
	  PoolObject() {}

      /// Destructor.
	  virtual ~PoolObject() {}


      ///@}
      ///@name Operators
      ///@{

	  void* operator new(std::size_t Size){
		  return MemoryPool::Allocate(Size);
	  }

	  void operator delete(void* pPointerToRelease, std::size_t Size){
		  MemoryPool::Deallocate(pPointerToRelease, Size);
	  }

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
	  virtual std::string Info() const {
		  return "PoolObject";
	  }

      /// Print information about this object.
	  virtual void PrintInfo(std::ostream& rOStream) const {
		  rOStream << Info();
	  }

      /// Print object's data.
	  virtual void PrintData(std::ostream& rOStream) const {

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
      PoolObject& operator=(PoolObject const& rOther);

      /// Copy constructor.
      PoolObject(PoolObject const& rOther);


      ///@}

    }; // Class PoolObject

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    PoolObject& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const PoolObject& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POOL_OBJECT_H_INCLUDED  defined
