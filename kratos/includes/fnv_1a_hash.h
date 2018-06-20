//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_FNV_1A_HASH_H_INCLUDED)
#define KRATOS_FNV_1A_HASH_H_INCLUDED

namespace Kratos {
///@addtogroup Kratos Core
///@{
///@name Kratos Classes
///@{

/**
 * @class FNV1a32Hash
 * @ingroup KratosCore
 * @brief A constexpr version of FNV hash function. (32 bit version)
 * @details The algorithm is the FNV-1a version of Fowler–Noll–Vo (FNV) hash function as
 *  described in Wikipedia.
 *  https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
 * @author Pooyan Dadvand
*/
class FNV1a32Hash {

public:
    ///@name Type Definitions
    ///@{

    /// The hash to be employed is 32 bits this time
    typedef std::uint32_t HashType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// The class is unconstructable.
    FNV1a32Hash() = delete;

    /// Destructor.
    virtual ~FNV1a32Hash() = delete;

    /// Assignment operator.
    FNV1a32Hash &operator=(FNV1a32Hash const &rOther) = delete;

    /// Copy constructor.
    FNV1a32Hash(FNV1a32Hash const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static constexpr HashType CalculateHash(const char *const TheString) {
        return CalculateHash(mFNV32OfsetBasis, TheString);
    }

    ///@}

    private:
    ///@name Static Member Variables
    ///@{

    static constexpr HashType mFNV32OfsetBasis = 0x811c9dc5;
    static constexpr HashType mFNV32Prime = 0x1000193;

    ///@}
    ///@name Private Operations
    ///@{
    static constexpr HashType CalculateHash(const HashType Value,
                                                const char *const TheString) {
        return (TheString[0] == '\0')
                ? Value
                : CalculateHash((Value ^ static_cast<HashType>(TheString[0])) * mFNV32Prime,
                                TheString + 1);
    }

  ///@}

}; // Class FNV1a32Hash

/**
 * @class FNV1a64Hash
 * @ingroup KratosCore
 * @brief A constexpr version of FNV hash function. (64 bit version)
 * @details The algorithm is the FNV-1a version of Fowler–Noll–Vo (FNV) hash function as
 *  described in Wikipedia.
 *  https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
 * @author Vicente Mataix Ferrandiz
*/
class FNV1a64Hash {

public:
    ///@name Type Definitions
    ///@{

    /// The hash to be employed is 64 bits this time
    typedef std::uint64_t HashType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// The class is unconstructable.
    FNV1a64Hash() = delete;

    /// Destructor.
    virtual ~FNV1a64Hash() = delete;

    /// Assignment operator.
    FNV1a64Hash &operator=(FNV1a64Hash const &rOther) = delete;

    /// Copy constructor.
    FNV1a64Hash(FNV1a64Hash const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static constexpr HashType CalculateHash(const char *const TheString) {
        return CalculateHash(mFNV64OfsetBasis, TheString);
    }

    ///@}

    private:
    ///@name Static Member Variables
    ///@{

    static constexpr HashType mFNV64OfsetBasis = 0xcbf29ce484222325;
    static constexpr HashType mFNV64Prime = 0x100000001b3;

    ///@}
    ///@name Private Operations
    ///@{
    static constexpr HashType CalculateHash(const HashType Value,
                                                const char *const TheString) {
        return (TheString[0] == '\0')
                ? Value
                : CalculateHash((Value ^ static_cast<HashType>(TheString[0])) * mFNV64Prime,
                                TheString + 1);
    }

    ///@}

}; // Class FNV1a64Hash

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FNV_1A_HASH_H_INCLUDED  defined
