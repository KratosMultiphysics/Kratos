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

// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "includes/serializer.h"

namespace Kratos {

void Flags::save(Serializer& rSerializer) const
{
    rSerializer.save("IsDefined",  mIsDefined);
    rSerializer.save("Flags",  mFlags);
}

void Flags::load(Serializer& rSerializer)
{
    rSerializer.load("IsDefined",  mIsDefined);
    rSerializer.load("Flags",  mFlags);
}

void Flags::Set(const Flags ThisFlag)
{
    mIsDefined |= ThisFlag.mIsDefined;
    mFlags = (mFlags & ~ThisFlag.mIsDefined) | (ThisFlag.mIsDefined & ThisFlag.mFlags);
}

void Flags::Set(const Flags ThisFlag, bool value)
{
    mIsDefined |= ThisFlag.mIsDefined;
    mFlags = (mFlags & ~ThisFlag.mIsDefined) | (ThisFlag.mIsDefined * BlockType(value));
}

bool operator==(const Flags& Left, const Flags& Right )
{
    return (Left.mFlags == Right.mFlags);
}

bool operator!=(const Flags& Left, const Flags& Right )
{
    return (Left.mFlags != Right.mFlags);
}

Flags KRATOS_API(KRATOS_CORE) operator|(const Flags& Left, const Flags& Right )
{
    Flags results(Left);
    results |= Right;
    return results;
}

Flags KRATOS_API(KRATOS_CORE) operator&(const Flags& Left, const Flags& Right )
{
    Flags results(Left);
    results &= Right;
    return results;
}

const Flags& Flags::operator|=(const Flags& Other )
{
    mIsDefined |= Other.mIsDefined;
    mFlags |= Other.mFlags;
    return *this;
}

const Flags& Flags::operator&=(const Flags& Other )
{
    mIsDefined |= Other.mIsDefined;
    mFlags &= Other.mFlags;
    return *this;
}


Flags::BlockType Flags::GetDefined() const
{
    return mIsDefined;
}

void Flags::SetDefined(const BlockType& rDefined)
{
    mIsDefined = rDefined;
}

Flags::BlockType Flags::GetFlags() const
{
    return mFlags;
}

void Flags::SetFlags(const BlockType& rValues)
{
    mFlags = rValues;
}

}
