//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#if !defined(KRATOS_GLOBAL_POINTER_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_H_INCLUDED

#include <iostream>

#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/key_hash.h"
namespace Kratos {



template<class TDataType>
class GlobalPointer {
private:

  /// Pointer to the data
  TDataType * mDataPointer;
#ifdef KRATOS_USING_MPI
  int mRank;
#endif

public:

  typedef TDataType element_type;

 /** Default constructor
	* Default constructor
	* This should never be called as we need a local pointer to exists
	*/
	GlobalPointer() {
    mDataPointer = nullptr;
#ifdef KRATOS_USING_MPI
    this->mRank = 0;
#endif
  };

  /** Constructor by Data
   * Constructor by Data
   * @param Data the data.
   */
  GlobalPointer(TDataType Data) = delete;

  /** Constructor by DataPointer
   * Constructor by DataPointer
   * @param DataPointer Pointer to the data.
   */
  GlobalPointer(TDataType * DataPointer, int Rank = 0)
    : mDataPointer(DataPointer)
#ifdef KRATOS_USING_MPI
    , mRank(Rank)
#endif
     {
#ifndef KRATOS_USING_MPI
     KRATOS_DEBUG_ERROR_IF(Rank != 0) << "trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
     }

  /** Constructor by Kratos::shared_ptr
   * Constructor by Kratos::shared_ptr
   * @param DataPointer Boost Shared Pointer to the Data.
   */
  GlobalPointer(Kratos::shared_ptr<TDataType> DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
    , mRank(Rank) 
#endif
  {
#ifndef KRATOS_USING_MPI
     KRATOS_DEBUG_ERROR_IF(Rank != 0) << "trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
  }

  /** Constructor by Kratos::shared_ptr
   * Constructor by Kratos::shared_ptr
   * @param DataPointer Boost Shared Pointer to the Data.
   */
  GlobalPointer(Kratos::intrusive_ptr<TDataType>& DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
    , mRank(Rank) 
#endif
  {
#ifndef KRATOS_USING_MPI
     KRATOS_DEBUG_ERROR_IF(Rank != 0) << "trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
  }

  /** Constructor by Kratos::weak_ptr
   * Constructor by Kratos::weak_ptr
   * @param DataPointer Kratos Weak Pointer to the Data.
   */
  GlobalPointer(Kratos::weak_ptr<TDataType> DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.lock().get())
#ifdef KRATOS_USING_MPI
    , mRank(Rank) 
#endif
  {
#ifndef KRATOS_USING_MPI
     KRATOS_DEBUG_ERROR_IF(Rank != 0) << "trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
  }

  /** Constructor by std::unique_ptr
   * Constructor by std::unique_ptr
   * @param DataPointer Std Unique Pointer to the Data.
   */
  GlobalPointer(std::unique_ptr<TDataType> DataPointer, int Rank = 0) = delete;

  /** Copy constructor
   * Copy constructor
   * @ rOther: Copied GlobalPointer
   */
  GlobalPointer(const GlobalPointer & rOther)
    : mDataPointer(rOther.mDataPointer)
#ifdef KRATOS_USING_MPI
    , mRank(rOther.mRank)
#endif
  {
  }

  /** Move constructor
   * Move constructor
   * @ rOther: Copied GlobalPointer
   */
  GlobalPointer(const GlobalPointer && rOther)
    : mDataPointer(std::move(rOther.mDataPointer))
#ifdef KRATOS_USING_MPI
    , mRank(std::move(rOther.mRank))
#endif
  {
  }

  /** Assignment Operator
   * Assignment Operator
   */
  GlobalPointer & operator=(const GlobalPointer & rOther) {
    mDataPointer = rOther.mDataPointer;
#ifdef KRATOS_USING_MPI
    mRank = rOther.mRank;
#endif

    return *this;
  }

  /** Default Destructor
   * Default Destructor.
   */
  ~GlobalPointer() {
  }

  TDataType* get() {
	  return mDataPointer;
  }

  TDataType const* get() const {
	  return mDataPointer;
  }

  /** Pointer Operator
  * Pointer Operator
  */
  TDataType & operator*() {
	  return *mDataPointer;
  }

  /** Const Pointer Operator
  * Const Pointer Operator
  */
  TDataType const& operator*() const {
	  return *mDataPointer;
  }

  /** Arrow Operator
   * Arrow Operator
   */
  TDataType * operator->() {
    return mDataPointer;
  }

  TDataType const* operator->() const {
    return mDataPointer;
  }

  /** Fills buffer with the GlobalPoiter data
   * Fills buffer with the GlobalPoiter data
   * @param buffer Object data buffer
   */
  void save(char * buffer) const {
    memcpy(buffer, this, sizeof(GlobalPointer));
  }

  /** Restores the GlobalPoiter with the data from the buffer
   * Restores the GlobalPoiter with the data from the buffer
   * @param buffer Object data buffer
   */
  void load(char * buffer) {
    memcpy(this, buffer, sizeof(GlobalPointer));
  }

  /** Returns the rank of the global pointer
   * Returns the rank of the global pointer data is located or 0 if no mpi
   * @return rank of the global pointer data or 0
   */
  int GetRank() const {
#ifdef KRATOS_USING_MPI
    return this->mRank;
#else
    return 0;
#endif
  }

  private:

  friend class Serializer;



  void save(Serializer& rSerializer) const
  {
      if(rSerializer.Is(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION))
      {
        rSerializer.save("D", reinterpret_cast<const std::size_t>(mDataPointer));
      }
      else
      {
        rSerializer.save("D", mDataPointer);
      }
#ifdef KRATOS_USING_MPI
      rSerializer.save("R", mRank);
#endif
  }

  void load(Serializer& rSerializer)
  {
      if(rSerializer.Is(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION))
      {
        std::size_t tmp;
        rSerializer.load("D", tmp);
        mDataPointer = reinterpret_cast<TDataType*>(tmp);
      }
      else
      {
        rSerializer.load("D", mDataPointer);
      }
#ifdef KRATOS_USING_MPI
      rSerializer.load("R", mRank);
#endif
  }

};

template< class TDataType >
struct GlobalPointerHasher
{
    /**
     * @brief The () operator
     * @param pDoF The DoF pointer
     */
    std::size_t operator()(const GlobalPointer<TDataType>& pGp) const
    {
        std::size_t seed = 0;
        HashCombine(seed, &(*pGp) );
#ifdef KRATOS_USING_MPI
        HashCombine(seed, pGp.GetRank());
#endif
        return seed;
    }
};

/**
 * @brief This is a key comparer between two dof pointers
 * @details Used for example for the B&S
 */
template< class TDataType >
struct GlobalPointerComparor
{
    /**
     * @brief The () operator
     * @param pDoF1 The first DoF pointer
     * @param pDoF2 The second DoF pointer
     */
    bool operator()(const GlobalPointer<TDataType>& pGp1, const GlobalPointer<TDataType>& pGp2) const
    {
#ifdef KRATOS_USING_MPI
        return ( &(*pGp1) == &(*pGp2)  &&  pGp1.GetRank() == pGp2.GetRank()  );
#else
        return ( &(*pGp1) == &(*pGp2) );
#endif
    }
};

/// input stream function
template< class TDataType >
inline std::istream& operator >> (std::istream& rIStream,
                                  GlobalPointer<TDataType>& rThis)
                                  {return rIStream;};

/// output stream function
template< class TDataType >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GlobalPointer<TDataType>& rThis)
{

    rOStream << reinterpret_cast<const std::size_t>(&*rThis) << " : " << rThis.GetRank();

    return rOStream;
}

} // namespace Kratos

#endif // KRATOS_GLOBAL_POINTER_H_INCLUDED
