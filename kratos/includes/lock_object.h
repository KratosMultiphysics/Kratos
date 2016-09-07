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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{
  
  /// This class defines stores a lock and gives interface to it.
  /** The class makes a tiny wrapper over the OpenMP Lock and also
  provides some auxiliary methods to obtain number of threads and 
  current thread number
  */
  class LockObject
    {
    public:
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
		LockObject() noexcept {
#ifdef _OPENMP
			omp_init_lock(&mLock);
#endif
		}

		/// Copy constructor.
		LockObject(LockObject const& rOther) noexcept
#ifdef _OPENMP
			: mLock(rOther.mLock)
#endif
		{
#ifdef _OPENMP
			omp_init_lock(&mLock);
#endif
		}

		/// Move constructor.
		LockObject(LockObject&& rOther) noexcept
#ifdef _OPENMP
			: mLock(rOther.mLock)
#endif
		{
#ifdef _OPENMP
			omp_init_lock(&mLock);
#endif
		}

		/// Destructor.
		virtual ~LockObject() noexcept {
#ifdef _OPENMP
			omp_destroy_lock(&mLock);
#endif
		}

      ///@}
      ///@name Operators 
      ///@{

	  /// Assignment operator.
		LockObject& operator=(LockObject const& rOther) {
#ifdef _OPENMP
			mLock = rOther.mLock;
#endif
			return *this;
		}

      
      ///@}
      ///@name Operations
      ///@{
      
		static inline int GetNumberOfThreads()
		{
#ifdef _OPENMP
			return omp_get_max_threads();
#else
			return 1;
#endif
		}

		static inline int GetThreadNumber()
		{
#ifdef _OPENMP
			return omp_get_thread_num();
#else
			return 0;
#endif
		}

      ///@}
      ///@name Access
      ///@{ 
      
#ifdef _OPENMP
		omp_lock_t& GetLock() const
		{
			return mLock;
		}
#endif

		inline void SetLock() const
		{
			//does nothing if openMP is not present
#ifdef _OPENMP
			omp_set_lock(&mLock);
#endif
		}

		inline void UnSetLock() const
		{
			//does nothing if openMP is not present
#ifdef _OPENMP
			omp_unset_lock(&mLock);
#endif
		}

      ///@}
      
    private:
      ///@name Member Variables 
      ///@{ 
        
#ifdef _OPENMP
		mutable omp_lock_t mLock;
#endif

      ///@}    
        
    }; // Class LockObject 

  ///@} 
  

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_LOCK_OBJECT_H_INCLUDED  defined 


