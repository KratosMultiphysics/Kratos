
//
//  Main authors:   Aron Noordam
//

// Application includes
#include "railway_application.h"

namespace Kratos
{
// We define the node type
using NodeType = Node;

KratosRailwayApplication::KratosRailwayApplication()
    : KratosApplication("KratosRailwayApplication")
{
}

void KratosRailwayApplication::Register()
{
    KratosApplication::Register();
}
} // namespace Kratos.
