//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <string_view>

// External includes
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes
#include "meshioplusplus_application.h"

#ifndef KRATOS_MESHIOPLUSPLUS_VERSION
#define KRATOS_MESHIOPLUSPLUS_VERSION "unknown"
#endif

namespace Kratos
{

KratosMeshioPlusPlusApplication::KratosMeshioPlusPlusApplication()
    : KratosApplication("MeshioPlusPlusApplication")
{
    // Logging information about the used library, for debugging purposes.
    // Note: the detail severity must be explicitly enabled, it is not shown by default.
    KRATOS_DETAIL("MeshioPlusPlusApplication")
        << "Version of the meshio++ library used during compilation: "
        << KRATOS_MESHIOPLUSPLUS_VERSION << std::endl;
    KRATOS_DETAIL("MeshioPlusPlusApplication")
        << "meshio++ mesh backend: " << meshioplusplus::mesh_backend_name()
        << " - parallel backend: " << meshioplusplus::parallel_backend_name() << std::endl;

    // The bridge to Kratos (meshioplusplus::ModelPart, kratos_bridge.hpp) only exists in
    // the KRATOS mesh backend. Linking meshioplusplus::core_kratos is what selects it, so
    // this is a build-configuration guard rather than a runtime condition.
    static_assert(std::string_view(meshioplusplus::mesh_backend_name()) == "kratos",
                  "MeshioPlusPlusApplication must be linked against meshioplusplus::core_kratos "
                  "(the KRATOS mesh backend); link that target instead of meshioplusplus::core.");
}

void KratosMeshioPlusPlusApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosMeshioPlusPlusApplication..." << std::endl;
}

} // namespace Kratos.
