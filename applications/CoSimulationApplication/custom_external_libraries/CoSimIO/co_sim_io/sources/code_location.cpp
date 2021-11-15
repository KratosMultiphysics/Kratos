//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// Project includes
#include "includes/code_location.hpp"
#include "includes/filesystem_inc.hpp"

namespace CoSimIO {
namespace Internals {

std::string CodeLocation::GetCleanFileName() const
{
    return fs::canonical(fs::path(GetFileName())).lexically_relative(fs::absolute(".")).string();
}

} // namespace Interals
} // namespace CoSimIO
