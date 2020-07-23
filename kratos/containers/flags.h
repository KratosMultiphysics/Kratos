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
//                   Riccardo Rossi
//
//

#if !defined(KRATOS_FLAGS_H_INCLUDED )
#define  KRATOS_FLAGS_H_INCLUDED



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

/// Short class definition.
/** Detail class definition.
*/
class Serializer;
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

    typedef int64_t BlockType;

    typedef int64_t FlagType;

    typedef std::size_t IndexType;

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

//    Flags(BlockType rOther) : mIsDefined(rOther), mFlags(rOther)
//    {
//    }

//    template<class TFlagsType>
//    Flags(TFlagsType rOther) : mIsDefined(static_cast<BlockType>(rOther)), mFlags(static_cast<BlockType>(rOther))
//    {
//    }

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


    operator bool() const
    {
        return mFlags;
    }

    Flags operator~() const
    {
        Flags results(*this);
        results.mFlags = ~mFlags;
        return  results;
    }

    bool operator!() const
    {
        Flags results(*this);
        results.mFlags = !mFlags;
        return  results;
    }


    void AssignFlags(Flags const& rOther)
    {
        mIsDefined = rOther.mIsDefined;
        mFlags = rOther.mFlags;
    }

    ///@}
    ///@name Operations
    ///@{

    void Set(Flags ThisFlag);

    void Set(Flags ThisFlag, bool Value);

    void Reset(Flags ThisFlag)
    {
        mIsDefined &= (~ThisFlag.mIsDefined);
        mFlags &= (~ThisFlag.mIsDefined); // I want to put to zero the ones are set regardless to their value. Pooyan.
    }

    void Flip(Flags ThisFlag)
    {
        mIsDefined |= ThisFlag.mIsDefined;
        mFlags ^= (ThisFlag.mIsDefined); // I want to flip  the ones are set in this flags regardless to their value. Pooyan.
    }

    void SetPosition(IndexType Position, bool Value=true )
    {
        mIsDefined |= (BlockType(1) << Position);


        mFlags &= ~(BlockType(1) << Position);   // First reseting the position
        mFlags |= (BlockType(Value) << Position);

    }

    bool GetPosition(IndexType Position) const
    {
        return (mFlags & (BlockType(1) << Position));
    }


   void FlipPosition(IndexType Position)
    {
        mIsDefined |= (BlockType(1) << Position);
        mFlags ^= (BlockType(1) << Position);
    }


    void ClearPosition(IndexType Position)
    {
        mIsDefined &= ~((BlockType(1) << Position));
        mFlags &= ~(BlockType(1) << Position);
    }


    void Clear()
    {
        mIsDefined = BlockType();
        mFlags = BlockType();
    }

    Flags AsFalse() const
    {
        Flags this_flag_false(*this);
        this_flag_false.mFlags = !((*this).mFlags);
        return this_flag_false;
    }

    ///@}
    ///@name Access
    ///@{

    static const Flags AllDefined()
    {
        return Flags(~0,0);
    }

    static const Flags AllTrue()
    {
        return Flags(~0,~0);
    }

    ///@}
    ///@name Inquiry
    ///@{


//    template<class TFlagsType>
//    bool Is(TFlagsType Flag)
//    {
//        return (mFlags & Flag);
//    }


    bool Is(Flags const & rOther) const
    {
        return (mFlags & rOther.mFlags) | ((rOther.mIsDefined ^ rOther.mFlags) & (~mFlags));
    }

    bool IsDefined(Flags const & rOther) const
    {
        return (mIsDefined & rOther.mIsDefined);
    }


//    template<class TFlagsType>
//    bool IsNot(TFlagsType const& Flag )
//    {
//        return !(mFlags & Flag);
//    }

    bool IsNot(Flags const& rOther) const
    {
        return !((mFlags & rOther.mFlags) | ((rOther.mIsDefined ^ rOther.mFlags) & (~mFlags)));
    }

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


    friend bool KRATOS_API(KRATOS_CORE) operator==(const Flags& Left, const Flags& Right );

    friend bool KRATOS_API(KRATOS_CORE) operator!=(const Flags& Left, const Flags& Right );

    friend Flags KRATOS_API(KRATOS_CORE) operator|(const Flags& Left, const Flags& Right );

    friend Flags KRATOS_API(KRATOS_CORE) operator&(const Flags& Left, const Flags& Right );

    const Flags& operator|=(const Flags& Other );

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

    BlockType mIsDefined;
    BlockType mFlags;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    friend class MPIDataCommunicator;

    BlockType GetDefined() const;

    void SetDefined(const BlockType& rDefined);

    BlockType GetFlags() const;

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

#endif // KRATOS_FLAGS_H_INCLUDED  defined
