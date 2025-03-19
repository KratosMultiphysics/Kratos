// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include <filesystem>

namespace Kratos::Testing
{

namespace Defaults
{

constexpr auto absolute_tolerance = 1.0e-12;
constexpr auto relative_tolerance = 1.0e-6;

} // namespace Defaults

class TestUtilities
{
public:
    static bool CompareFiles(const std::filesystem::path& rPath1, const std::filesystem::path& rPath2);
};

class ScopedSerializerRegistrationOfAllStressStatePolicies
{
public:
    ScopedSerializerRegistrationOfAllStressStatePolicies();
    ~ScopedSerializerRegistrationOfAllStressStatePolicies();
    ScopedSerializerRegistrationOfAllStressStatePolicies(const ScopedSerializerRegistrationOfAllStressStatePolicies&) = delete;
    ScopedSerializerRegistrationOfAllStressStatePolicies& operator=(
        const ScopedSerializerRegistrationOfAllStressStatePolicies&) = delete;
    ScopedSerializerRegistrationOfAllStressStatePolicies(ScopedSerializerRegistrationOfAllStressStatePolicies&&) noexcept = default;
    ScopedSerializerRegistrationOfAllStressStatePolicies& operator=(
        ScopedSerializerRegistrationOfAllStressStatePolicies&&) noexcept = default;
};

class ScopedSerializerRegistrationOfAllYieldSurfaces
{
public:
    ScopedSerializerRegistrationOfAllYieldSurfaces();
    ~ScopedSerializerRegistrationOfAllYieldSurfaces();
    ScopedSerializerRegistrationOfAllYieldSurfaces(const ScopedSerializerRegistrationOfAllYieldSurfaces&) = delete;
    ScopedSerializerRegistrationOfAllYieldSurfaces& operator=(const ScopedSerializerRegistrationOfAllYieldSurfaces&) = delete;
    ScopedSerializerRegistrationOfAllYieldSurfaces(ScopedSerializerRegistrationOfAllYieldSurfaces&&) noexcept = default;
    ScopedSerializerRegistrationOfAllYieldSurfaces& operator=(ScopedSerializerRegistrationOfAllYieldSurfaces&&) noexcept = default;
};

} // namespace Kratos::Testing
