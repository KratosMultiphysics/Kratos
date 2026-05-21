//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#pragma once

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/key_hash.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
* @class GlobalPointer
* @ingroup KratosCore
* @brief This class is a wrapper for a pointer to a data that is located in a different rank
* @details This class is a wrapper for a pointer to a data that is located in a different rank
* @tparam TDataType The type of the data pointed by the GlobalPointer
* @see GlobalPointerHasher
* @author Carlos A. Roig
*/
template<class TDataType>
class GlobalPointer
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of element type
    using element_type = TDataType;

    ///@}
    ///@name Life Cycle
    ///@{

   /**
    * @brief Default constructor
    * @details This should never be called as we need a local pointer to exists
    */
    GlobalPointer() {
        mDataPointer = nullptr;
#ifdef KRATOS_USING_MPI
        this->mRank = 0;
#endif
    };

    /**
     * @brief Constructor by Data
     * @param Data the data.
     */
    GlobalPointer(TDataType Data) = delete;

    /**
     * @brief Constructor by DataPointer
     * @param DataPointer Pointer to the data.
     */
    GlobalPointer(
        TDataType* DataPointer,
        int Rank = 0
        ) : mDataPointer(DataPointer)
#ifdef KRATOS_USING_MPI
          , mRank(Rank)
#endif
    {
#ifndef KRATOS_USING_MPI
        KRATOS_DEBUG_ERROR_IF(Rank != 0) << "Trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
    }

    /**
     * @brief Constructor by Kratos::shared_ptr
     * @param DataPointer Std Shared Pointer to the Data.
     */
    GlobalPointer(
        Kratos::shared_ptr<TDataType> DataPointer,
        int Rank = 0
        ) : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
          , mRank(Rank)
#endif
    {
#ifndef KRATOS_USING_MPI
        KRATOS_DEBUG_ERROR_IF(Rank != 0) << "Trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
    }

    /**
     * @brief Constructor by Kratos::intrusive_ptr
     * @param DataPointer Boost intrusive Pointer to the Data.
     */
    GlobalPointer(
        Kratos::intrusive_ptr<TDataType>& DataPointer,
        int Rank = 0
        ) : mDataPointer(DataPointer.get())
    #ifdef KRATOS_USING_MPI
          , mRank(Rank)
    #endif
    {
#ifndef KRATOS_USING_MPI
        KRATOS_DEBUG_ERROR_IF(Rank != 0) << "Trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
#endif
    }

    /**
     * @brief Constructor by Kratos::weak_ptr
     * @param DataPointer Kratos Weak Pointer to the Data.
     */
    GlobalPointer(
        Kratos::weak_ptr<TDataType> DataPointer,
        int Rank = 0
        ) : mDataPointer(DataPointer.lock().get())
#ifdef KRATOS_USING_MPI
          , mRank(Rank)
#endif
    {
  #ifndef KRATOS_USING_MPI
        KRATOS_DEBUG_ERROR_IF(Rank != 0) << "Trying to construct a global pointer with rank different from zero when kratos is not in MPI mode " << std::endl;
  #endif
    }

    /**
     * @brief Constructor by std::unique_ptr
     * @param DataPointer Std Unique Pointer to the Data.
     */
    GlobalPointer(std::unique_ptr<TDataType> DataPointer, int Rank = 0) = delete;

    /**
     * @brief Copy constructor for GlobalPointer.
     * @param rOther The GlobalPointer to copy from.
     */
    GlobalPointer(const GlobalPointer& rOther)
        : mDataPointer(rOther.mDataPointer)
#ifdef KRATOS_USING_MPI
        , mRank(rOther.mRank)
#endif
    {
    }

    /**
     * @brief Move constructor for GlobalPointer.
     * @param rOther The GlobalPointer to move from.
     */
    GlobalPointer(const GlobalPointer&& rOther)
        : mDataPointer(std::move(rOther.mDataPointer))
#ifdef KRATOS_USING_MPI
        , mRank(std::move(rOther.mRank))
#endif
    {
    }

    /// Destructor.
    ~GlobalPointer() = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment operator for GlobalPointer.
     * @param rOther The GlobalPointer to assign from.
     * @return Reference to the assigned GlobalPointer.
     */
    GlobalPointer& operator=(const GlobalPointer& rOther) {
        mDataPointer = rOther.mDataPointer;
#ifdef KRATOS_USING_MPI
        mRank = rOther.mRank;
#endif
        return *this;
    }

    /**
     * @brief Dereference operator for non-const GlobalPointer.
     * @return Reference to the data pointed to by GlobalPointer.
     */
    TDataType& operator*() {
        return *mDataPointer;
    }

    /**
     * @brief Dereference operator for const GlobalPointer.
     * @return Const reference to the data pointed to by GlobalPointer.
     */
    TDataType const& operator*() const {
        return *mDataPointer;
    }

    /**
     * @brief Arrow operator for non-const GlobalPointer.
     * @return Pointer to the data pointed to by GlobalPointer.
     */
    TDataType* operator->() {
        return mDataPointer;
    }

    /**
     * @brief Arrow operator for const GlobalPointer.
     * @return Const pointer to the data pointed to by GlobalPointer.
     */
    TDataType const* operator->() const {
        return mDataPointer;
    }

    /**
     * @brief Overloads the '==' operator to compare two GlobalPointer objects of the same template type. Returns true if the underlying pointers are equal.
     * @param rOther The GlobalPointer object to be compared.
     * @return true if the underlying pointers are equal, false otherwise.
     */
    bool operator==(const GlobalPointer& rOther)
    {
#ifdef KRATOS_USING_MPI
      return this->get() == rOther.get() && this->GetRank() == rOther.GetRank();
#else
      return this->get() == rOther.get();
#endif
    }

    ///@}
    ///@name Access
    ///@{

    /**
    * @brief Get a non-const pointer to the data.
    * @return A pointer to the data.
    */
    TDataType* get() {
        return mDataPointer;
    }

    /**
    * @brief Get a const pointer to the data.
    * @return A const pointer to the data.
    */
    TDataType const* get() const {
        return mDataPointer;
    }

    /**
     * @brief Returns the rank of the global pointer
     * @details Returns the rank of the global pointer data is located or 0 if no mpi
     * @return rank of the global pointer data or 0
     */
    int GetRank() const
    {
#ifdef KRATOS_USING_MPI
      return this->mRank;
#else
      return 0;
#endif
    }

    /**
     * @brief Sets the rank of the global pointer
     * @details Sets the rank of the global pointer where data is located
     * @param Rank The rank to set
     */
    void SetRank(const int Rank)
    {
#ifdef KRATOS_USING_MPI
      this->mRank = Rank;
#endif
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Fills buffer with the GlobalPoiter data
     * @param buffer Object data buffer
     */
    void save(char * buffer) const {
      memcpy(buffer, this, sizeof(GlobalPointer));
    }

    /**
     * @brief Restores the GlobalPoiter with the data from the buffer
     * @param buffer Object data buffer
     */
    void load(char * buffer) {
      memcpy(this, buffer, sizeof(GlobalPointer));
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

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GlobalPointer" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointer";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointer from Rank: " << GetRank() << " contains: \n";
        mDataPointer->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    TDataType* mDataPointer; /// Pointer to the data
#ifdef KRATOS_USING_MPI
    int mRank;               /// Rank of the data
#endif

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
    ///@name Serialization
    ///@{

    /**
     * @brief Friend class Serializer for handling serialization.
     */
    friend class Serializer;

    /**
     * @brief Save function for serializing the data.
     * @param rSerializer The serializer to save data into.
     */
    void save(Serializer& rSerializer) const
    {
        if(rSerializer.Is(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION)) {
            rSerializer.save("D", reinterpret_cast<std::size_t>(mDataPointer));
        } else {
            rSerializer.save("D", mDataPointer);
        }
    #ifdef KRATOS_USING_MPI
        rSerializer.save("R", mRank);
    #endif
    }

    /**
     * @brief Load function for deserializing the data.
     * @param rSerializer The serializer to load data from.
     */
    void load(Serializer& rSerializer)
    {
        if(rSerializer.Is(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION)) {
            std::size_t tmp;
            rSerializer.load("D", tmp);
            mDataPointer = reinterpret_cast<TDataType*>(tmp);
        } else {
            rSerializer.load("D", mDataPointer);
        }
    #ifdef KRATOS_USING_MPI
        rSerializer.load("R", mRank);
    #endif
    }

    ///@}
};

/**
 * @brief Template struct for hashing GlobalPointer instances.
 * @tparam TDataType The type of data managed by the GlobalPointer.
 */
template <class TDataType>
struct GlobalPointerHasher
{
    /**
     * @brief The () operator for hashing a GlobalPointer.
     * @param pGp The GlobalPointer to be hashed.
     * @return The hash value of the GlobalPointer.
     */
    std::size_t operator()(const GlobalPointer<TDataType>& pGp) const
    {
        std::size_t seed = 0;
        HashCombine(seed, &(*pGp));
#ifdef KRATOS_USING_MPI
        HashCombine(seed, pGp.GetRank());
#endif
        return seed;
    }
};

/**
 * @brief This is a key comparer between two dof pointers checking for equal keys
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

/**
 * @brief This is a key compare between two pointers to the object object
 * @details This should be used in std::sort or any other algorithm requiring weak ordering
 * https://en.cppreference.com/w/cpp/named_req/Compare
 */
template< class TDataType >
struct GlobalPointerCompare
{
    /**
     * @brief The () operator
     * @param pDoF1 The first DoF pointer
     * @param pDoF2 The second DoF pointer
     */
    bool operator()(const GlobalPointer<TDataType>& pGp1, const GlobalPointer<TDataType>& pGp2) const
    {
#ifdef KRATOS_USING_MPI
        return (pGp1.GetRank() == pGp2.GetRank()) ? (pGp1.get() < pGp2.get()) : (pGp1.GetRank() < pGp2.GetRank());
#else
        return (pGp1.get() < pGp2.get());
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
