//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "includes/schema.hpp"
#include "includes/kratos_parameters.h"


namespace Kratos {


struct Schema::Impl
{
    Parameters mDefinition;
}; // struct Schema::Impl


Schema::Schema()
    : Schema(Parameters())
{
}


Schema::Schema(Schema&&) noexcept = default;


Schema::Schema(Parameters&& rDefinition)
    : mpImpl(new Impl {/*mDefinition:*/std::move(rDefinition)})
{
    if (not mpImpl->mDefinition.Has("$schema"))
        mpImpl->mDefinition.AddString("$schema", this->GetJSONSchemaStandard());
}


Schema::Schema(const Schema& rOther)
    : mpImpl(new Impl {/*mDefinition:*/Parameters(rOther.mpImpl->mDefinition)})
{
}


Schema::~Schema() = default;


Schema& Schema::operator=(Schema&&) noexcept = default;


Schema& Schema::operator=(const Schema& rOther)
{
    mpImpl->mDefinition = rOther.mpImpl->mDefinition;
    return *this;
}


Parameters& Schema::GetDefinition() noexcept
{
    return mpImpl->mDefinition;
}

const Parameters& Schema::GetDefinition() const noexcept
{
    return mpImpl->mDefinition;
}


Schema::operator Parameters& () noexcept
{
    return mpImpl->mDefinition;
}


Schema::operator const Parameters& () const noexcept
{
    return mpImpl->mDefinition;
}


Schema::operator std::string () const
{
    return mpImpl->mDefinition.PrettyPrintJsonString();
}


std::string Schema::GetJSONSchemaStandard()
{
    return "https://json-schema.org/draft/2020-12/schema";
}


std::ostream& operator<<(std::ostream& rStream, const Schema& rSchema)
{
    return rStream << (std::string)rSchema;
}


} // namespace Kratos
