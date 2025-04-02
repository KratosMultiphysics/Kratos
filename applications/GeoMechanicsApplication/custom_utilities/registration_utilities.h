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
    template <typename T>
    ScopedSerializerRegistration(const std::string& rName, const T& rInstance) : mName(rName)
    {
        Serializer::Register(rName, rInstance);
    }

    ~ScopedSerializerRegistration();

    ScopedSerializerRegistration(const ScopedSerializerRegistration&)            = delete;
    ScopedSerializerRegistration& operator=(const ScopedSerializerRegistration&) = delete;

    ScopedSerializerRegistration(ScopedSerializerRegistration&&) noexcept            = default;
    ScopedSerializerRegistration& operator=(ScopedSerializerRegistration&&) noexcept = default;

private:
    std::string mName;
};

} // namespace Kratos