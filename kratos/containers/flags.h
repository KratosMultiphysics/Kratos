//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"

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

class Serializer; /// Forward declaration

/**
 * @class Flags
 * @ingroup KratosCore
 * @brief Defines flags used in the Kratos framework for efficient storage and manipulation of state information.
 * @details Flags operate using bitwise logic, providing a crucial base class for numerous essential Kratos components.
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Flags
    KRATOS_CLASS_POINTER_DEFINITION(Flags);

#ifdef  _WIN32 // work around for windows int64_t error
    typedef __int64 int64_t;
#endif

    using BlockType = int64_t;

    using FlagType = int64_t;

    using IndexType = std::size_t;

    enum FlagsList
    {
        Flag0 = BlockType(1),
        Flag1 = BlockType(1) << 1,
        Flag2 = BlockType(1) << 2,
        Flag3 = BlockType(1) << 3,
        Flag4 = BlockType(1) << 4,
        Flag5 = BlockType(1) << 5,
        Flag6 = BlockType(1) << 6,
        Flag7 = BlockType(1) << 7,
        Flag8 = BlockType(1) << 8,
        Flag9 = BlockType(1) << 9,

        Flag10 = BlockType(1) << 10,
        Flag11 = BlockType(1) << 11,
        Flag12 = BlockType(1) << 12,
        Flag13 = BlockType(1) << 13,
        Flag14 = BlockType(1) << 14,
        Flag15 = BlockType(1) << 15,
        Flag16 = BlockType(1) << 16,
        Flag17 = BlockType(1) << 17,
        Flag18 = BlockType(1) << 18,
        Flag19 = BlockType(1) << 19,

        Flag20 = BlockType(1) << 20,
        Flag21 = BlockType(1) << 21,
        Flag22 = BlockType(1) << 22,
        Flag23 = BlockType(1) << 23,
        Flag24 = BlockType(1) << 24,
        Flag25 = BlockType(1) << 25,
        Flag26 = BlockType(1) << 26,
        Flag27 = BlockType(1) << 27,
        Flag28 = BlockType(1) << 28,
        Flag29 = BlockType(1) << 29,

        Flag30 = BlockType(1) << 30//,
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Flags() : mIsDefined(BlockType()), mFlags(BlockType()) {}

    /// Copy constructor.
    Flags(Flags const& rOther) : mIsDefined(rOther.mIsDefined), mFlags(rOther.mFlags)
    {
    }

    /// Destructor.
    virtual ~Flags() {}

    static Flags Create(IndexType ThisPosition, bool Value=true)
    {
        Flags flags;
        flags.SetPosition(ThisPosition, Value);
        return flags;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Flags& operator=(Flags const& rOther)
    {
        mIsDefined = rOther.mIsDefined;
        mFlags = rOther.mFlags;
        return *this;
    }

    /**
     * @brief Conversion operator to bool.
     * @return true if any flag is set, false otherwise.
     */
    operator bool() const
    {
        return mFlags;
    }

    /**
     * @brief Bitwise NOT operator.
     * @return Flags with each bit inverted.
     */
    Flags operator~() const
    {
        Flags results(*this);
        results.mFlags = ~mFlags;
        return  results;
    }

    /**
     * @brief Logical NOT operator.
     * @return Flags with logical negation applied.
     */
    bool operator!() const
    {
        Flags results(*this);
        results.mFlags = !mFlags;
        return  results;
    }

    /**
     * @brief Assign flags from another Flags object.
     * @param rOther The Flags object to assign from.
     */
    void AssignFlags(Flags const& rOther)
    {
        mIsDefined = rOther.mIsDefined;
        mFlags = rOther.mFlags;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set the specified flag.
     * @param ThisFlag The flag to set.
     */
    void Set(const Flags ThisFlag);

    /**
     * @brief Set the specified flag with a given value.
     * @param ThisFlag The flag to set.
     * @param Value The value to set the flag to.
     */
    void Set(const Flags ThisFlag, bool Value);

    /**
     * @brief Reset the specified flag.
     * @param ThisFlag The flag to reset.
     */
    void Reset(const Flags ThisFlag)
    {
        mIsDefined &= (~ThisFlag.mIsDefined);
        mFlags &= (~ThisFlag.mIsDefined); // I want to put to zero the ones are set regardless to their value. Pooyan.
    }

    /**
     * @brief Flip the specified flag.
     * @param ThisFlag The flag to flip.
     */
    void Flip(const Flags ThisFlag)
    {
        mIsDefined |= ThisFlag.mIsDefined;
        mFlags ^= (ThisFlag.mIsDefined); // I want to flip  the ones are set in this flags regardless to their value. Pooyan.
    }

    /**
     * @brief Set the value of the flag at the specified position.
     * @param Position The position of the flag to set.
     * @param Value The value to set the flag to (default is true).
     */
    void SetPosition(IndexType Position, const bool Value=true )
    {
        mIsDefined |= (BlockType(1) << Position);

        mFlags &= ~(BlockType(1) << Position);   // First reseting the position
        mFlags |= (BlockType(Value) << Position);
    }

    /**
     * @brief Get the value of the flag at the specified position.
     * @param Position The position of the flag to get.
     * @return True if the flag is set, false otherwise.
     */
    bool GetPosition(IndexType Position) const
    {
        return (mFlags & (BlockType(1) << Position));
    }

    /**
     * @brief Flip the value of the flag at the specified position.
     * @param Position The position of the flag to flip.
     */
    void FlipPosition(IndexType Position)
    {
        mIsDefined |= (BlockType(1) << Position);
        mFlags ^= (BlockType(1) << Position);
    }

    /**
     * @brief Clear the flag at the specified position.
     * @param Position The position of the flag to clear.
     */
    void ClearPosition(IndexType Position)
    {
        mIsDefined &= ~((BlockType(1) << Position));
        mFlags &= ~(BlockType(1) << Position);
    }

    /**
     * @brief Clear all flags.
     */
    void Clear()
    {
        mIsDefined = BlockType();
        mFlags = BlockType();
    }

    /**
     * @brief Return a Flags object with all flags set to false.
     * @return Flags object with all flags set to false.
     */
    Flags AsFalse() const
    {
        Flags this_flag_false(*this);
        this_flag_false.mFlags = !((*this).mFlags);
        return this_flag_false;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Returns a Flags object with all flags defined.
     * @return Flags object with all flags defined.
     */
    static const Flags AllDefined()
    {
        return Flags(~0, 0);
    }

    /**
     * @brief Returns a Flags object with all flags set to true. 
     * @return Flags object with all flags set to true.
     */
    static const Flags AllTrue()
    {
        return Flags(~0, ~0);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if the flags in this object match those in the given Flags object.
     * @param rOther The Flags object to compare against.
     * @return True if the flags match, false otherwise.
     */
    bool Is(Flags const & rOther) const
    {
        return (mFlags & rOther.mFlags) | ((rOther.mIsDefined ^ rOther.mFlags) & (~mFlags));
    }

    /**
     * @brief Checks if all flags in this object are defined (set).
     * @param rOther The Flags object to compare against.
     * @return True if all flags are defined, false otherwise.
     */
    bool IsDefined(Flags const & rOther) const
    {
        return (mIsDefined & rOther.mIsDefined);
    }

    /**
     * @brief Checks if the flags in this object do not match those in the given Flags object.
     * @param rOther The Flags object to compare against.
     * @return True if the flags do not match, false otherwise.
     */
    bool IsNot(Flags const& rOther) const
    {
        return !((mFlags & rOther.mFlags) | ((rOther.mIsDefined ^ rOther.mFlags) & (~mFlags)));
    }

    /**
     * @brief Checks if any flags in this object are not defined (unset).
     * @param rOther The Flags object to compare against.
     * @return True if any flags are not defined, false otherwise.
     */
    bool IsNotDefined(Flags const& rOther) const
    {
        return !(mIsDefined & rOther.mIsDefined);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Flags" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Flags";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for(std::size_t i = sizeof(BlockType) * 8 ; i > 0 ; i--)
            rOStream << bool(mFlags & (BlockType(1) << i));
    }

    ///@}
    ///@name Friends
    ///@{

    /**
     * @brief Equality comparison operator for Flags.
     * @param Left The left-hand side Flags object.
     * @param Right The right-hand side Flags object.
     * @return True if the Flags objects are equal, false otherwise.
     */
    friend bool KRATOS_API(KRATOS_CORE) operator==(const Flags& Left, const Flags& Right );

    /**
     * @brief Inequality comparison operator for Flags.
     * @param Left The left-hand side Flags object.
     * @param Right The right-hand side Flags object.
     * @return True if the Flags objects are not equal, false otherwise.
     */
    friend bool KRATOS_API(KRATOS_CORE) operator!=(const Flags& Left, const Flags& Right );

    /**
     * @brief Bitwise OR operator for Flags.
     * @param Left The left-hand side Flags object.
     * @param Right The right-hand side Flags object.
     * @return Flags object resulting from the bitwise OR operation.
     */
    friend Flags KRATOS_API(KRATOS_CORE) operator|(const Flags& Left, const Flags& Right );

    /**
     * @brief Bitwise AND operator for Flags.
     * @param Left The left-hand side Flags object.
     * @param Right The right-hand side Flags object.
     * @return Flags object resulting from the bitwise AND operation.
     */
    friend Flags KRATOS_API(KRATOS_CORE) operator&(const Flags& Left, const Flags& Right );

    /**
     * @brief Compound assignment operator for bitwise OR.
     * @param Other The Flags object to perform bitwise OR with.
     * @return Reference to the modified Flags object after bitwise OR.
     */
    const Flags& operator|=(const Flags& Other );

    /**
     * @brief Compound assignment operator for bitwise AND.
     * @param Other The Flags object to perform bitwise AND with.
     * @return Reference to the modified Flags object after bitwise AND.
     */
    const Flags& operator&=(const Flags& Other );

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

    BlockType mIsDefined; /// Bitmask representing defined flags.
    BlockType mFlags; /// Bitmask representing flag values.

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    friend class MPIDataCommunicator; ///< Friend class for accessing private members.

    /**
     * @brief Get the bitmask representing defined flags.
     * @return BlockType representing defined flags.
     */
    BlockType GetDefined() const;

    /**
     * @brief Set the bitmask representing defined flags.
     * @param rDefined The BlockType to set as defined flags.
     */
    void SetDefined(const BlockType& rDefined);

    /**
     * @brief Get the bitmask representing flag values.
     * @return BlockType representing flag values.
     */
    BlockType GetFlags() const;

    /**
     * @brief Set the bitmask representing flag values.
     * @param rValues The BlockType to set as flag values.
     */
    void SetFlags(const BlockType& rValues);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    Flags(BlockType DefinedFlags, BlockType SetFlags):
        mIsDefined(DefinedFlags), mFlags(SetFlags)
    {}

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class Flags

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Flags& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Flags& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
