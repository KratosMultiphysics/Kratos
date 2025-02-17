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

#pragma once

// Project includes
#include "includes/smart_pointers.h" // KRATOS_CLASS_POINTER_DEFINITION
#include "includes/kratos_export_api.h" // KRATOS_API

// System includes
#include <memory> // std::unique_ptr
#include <string> // std::string
#include <iosfwd> // std::ostream


namespace Kratos {


class Parameters;


class KRATOS_API(KRATOS_CORE) Schema
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(Schema);

    Schema();

    Schema(Parameters&& rDefinition);

    Schema(Schema&&) noexcept;

    Schema(const Schema&);

    ~Schema();

    Schema& operator=(Schema&&) noexcept;

    Schema& operator=(const Schema&);

    Parameters& GetDefinition() noexcept;

    const Parameters& GetDefinition() const noexcept;

    explicit operator Parameters& () noexcept;

    explicit operator const Parameters& () const noexcept;

    explicit operator std::string () const;

    static std::string GetJSONSchemaStandard();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class Schema


std::ostream& operator<<(std::ostream& rStream, const Schema& rSchema);


} // namespace Kratos
