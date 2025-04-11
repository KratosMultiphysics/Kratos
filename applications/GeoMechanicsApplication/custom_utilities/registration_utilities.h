// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include "includes/kratos_export_api.h"
#include "includes/serializer.h"

#include <string>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ScopedSerializerRegistration
{
public:
    template <typename... Pairs>
    explicit ScopedSerializerRegistration(Pairs... PackedNameAndInstancePairs)
    {
        Register(PackedNameAndInstancePairs...);
    }

    ~ScopedSerializerRegistration();

    ScopedSerializerRegistration(const ScopedSerializerRegistration&)            = delete;
    ScopedSerializerRegistration& operator=(const ScopedSerializerRegistration&) = delete;

    ScopedSerializerRegistration(ScopedSerializerRegistration&&) noexcept            = default;
    ScopedSerializerRegistration& operator=(ScopedSerializerRegistration&&) noexcept = default;

private:
    std::vector<std::string> mNames;

    template <typename Pair>
    void Register(const Pair& rNameAndInstance)
    {
        const auto& [name, instance] = rNameAndInstance;
        Serializer::Register(name, instance);
        mNames.push_back(name);
    }

    template <typename Pair, typename... Pairs>
    void Register(const Pair& rNameAndInstance, Pairs... PackedNameAndInstancePairs)
    {
        Register(rNameAndInstance);
        Register(PackedNameAndInstancePairs...);
    }
};

} // namespace Kratos