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

#include "containers/flags.h"

namespace Kratos {


Flags operator|(const Flags& Left, const Flags& Right )
{
    Flags results(Left);
    results |= Right;
    return results;
}

Flags operator&(const Flags& Left, const Flags& Right )
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
    mFlags |= ~mIsDefined;
    mFlags &= Other.mFlags | ~Other.mIsDefined;
    mIsDefined |= Other.mIsDefined;
    mFlags &= mIsDefined;
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

Flags::BlockType Flags::GetValues() const
{
    return mFlags;
}

void Flags::SetValues(const BlockType& rValues)
{
    mFlags = rValues;
}

}