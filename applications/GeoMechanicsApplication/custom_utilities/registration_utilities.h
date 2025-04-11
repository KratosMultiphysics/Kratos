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
    template <typename... Ts>
    explicit ScopedSerializerRegistration(Ts... args)
    {
        Register(args...);
    }

    ~ScopedSerializerRegistration();

    ScopedSerializerRegistration(const ScopedSerializerRegistration&)            = delete;
    ScopedSerializerRegistration& operator=(const ScopedSerializerRegistration&) = delete;

    ScopedSerializerRegistration(ScopedSerializerRegistration&&) noexcept            = default;
    ScopedSerializerRegistration& operator=(ScopedSerializerRegistration&&) noexcept = default;

private:
    std::vector<std::string> mNames;

    template <class T>
    void Register(const T& rNameAndInstance)
    {
        const auto& [name, instance] = rNameAndInstance;
        Serializer::Register(name, instance);
        mNames.push_back(name);
    }

    template <class T, typename... Ts>
    void Register(const T& rNameAndInstance, Ts... args)
    {
        Register(rNameAndInstance);
        Register(args...);
    }
};

} // namespace Kratos