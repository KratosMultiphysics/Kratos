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
//
//

#if !defined(KRATOS_FNV_1A_HASH_H_INCLUDED)
#define KRATOS_FNV_1A_HASH_H_INCLUDED

namespace Kratos {
///@addtogroup Kratos Core
///@{
///@name Kratos Classes
///@{

/// A constexpr version of FNV hash function.
/** The algorithm is the FNV-1a version of Fowler–Noll–Vo (FNV) hash function as
 *  described in Wikipedia.
 *  https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
*/
class FNV1a32Hash {

public:
  ///@name Type Definitions
  ///@{

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

  static constexpr std::uint32_t CalculateHash(const char *const TheString) {
    return CalculateHash(mFNV32OfsetBasis, TheString);
  }

  ///@}

private:
  ///@name Static Member Variables
  ///@{

  static constexpr std::uint32_t mFNV32OfsetBasis = 0x811c9dc5;
  static constexpr std::uint32_t mFNV32Prime = 0x1000193;

  ///@}
  ///@name Private Operations
  ///@{
  static constexpr std::uint32_t CalculateHash(const std::uint32_t Value,
                                               const char *const TheString) {
    return (TheString[0] == '\0')
               ? Value
               : CalculateHash((Value ^ static_cast<std::uint32_t>(TheString[0])) * mFNV32Prime,
                               TheString + 1);
  }

  ///@}

}; // Class FNV1a32Hash

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FNV_1A_HASH_H_INCLUDED  defined
