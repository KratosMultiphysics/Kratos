//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher (https://github.com/philbucher)
//

// Project includes
#include "includes/serializer.hpp"

namespace CoSimIO {
namespace Internals {

Serializer::RegisteredObjectsContainerType Serializer::msRegisteredObjects;
Serializer::RegisteredObjectsNameContainerType Serializer::msRegisteredObjectsName;

} // namespace Internals
}
