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

#if !defined(KRATOS_LOCK_OBJECT_H_INCLUDED )
#define  KRATOS_LOCK_OBJECT_H_INCLUDED

// System includes

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// This class defines stores a lock and gives interface to it.
  /** The class makes a tiny wrapper over the ahared memory locking mechanisms
  */
  class LockObject
    {
    public:
      ///@name Life Cycle
      ///@{

      /// Default constructor.
		LockObject() noexcept {
#ifdef KRATOS_SMP_OPENMP
			omp_init_lock(&mLock);
#endif
		}

		/// Copy constructor.
		LockObject(LockObject const& rOther) noexcept
#ifdef KRATOS_SMP_OPENMP
			: mLock(rOther.mLock)
#endif
		{
#ifdef KRATOS_SMP_OPENMP
			omp_init_lock(&mLock);
#endif
		}

		/// Move constructor.
		LockObject(LockObject&& rOther) noexcept
#ifdef KRATOS_SMP_OPENMP
			: mLock(rOther.mLock)
#endif
		{
#ifdef KRATOS_SMP_OPENMP
			omp_init_lock(&mLock);
#endif
		}

		/// Destructor.
		virtual ~LockObject() noexcept {
#ifdef KRATOS_SMP_OPENMP
			omp_destroy_lock(&mLock);
#endif
		}

      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator.
		LockObject& operator=(LockObject const& rOther) {
#ifdef KRATOS_SMP_OPENMP
			mLock = rOther.mLock;
#endif
			return *this;
		}

      ///@}
      ///@name Operations
      ///@{

      ///@}
      ///@name Access
      ///@{

		inline void SetLock() const
		{
			//does nothing if openMP is not present
#ifdef KRATOS_SMP_OPENMP
			omp_set_lock(&mLock);
#endif
		}

		inline void UnSetLock() const
		{
			//does nothing if openMP is not present
#ifdef KRATOS_SMP_OPENMP
			omp_unset_lock(&mLock);
#endif
		}

      ///@}

    private:
      ///@name Member Variables
      ///@{

#ifdef KRATOS_SMP_OPENMP
		mutable omp_lock_t mLock;
#endif

      ///@}

    }; // Class LockObject

  ///@}


  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LOCK_OBJECT_H_INCLUDED  defined
