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
#include "utilities/dynamic_proxies.h" // DynamicEntityProxy
#include "includes/global_variables.h"


namespace Kratos {


namespace {

std::string GetDataLocationEnumName(Globals::DataLocation Location)
{
    std::string output;
    switch (Location) {
        case Globals::DataLocation::NodeHistorical: output = "NodeHistorical";
        case Globals::DataLocation::NodeNonHistorical: output = "NodeNonHistorical";
        case Globals::DataLocation::Element: output = "Element";
        case Globals::DataLocation::Condition: output = "Condition";
        default: KRATOS_ERROR << "Unsupported DataLocation enum";
    };
    return output;
}

} // unnamed namespace


DynamicEntityProxy::DynamicEntityProxy(Globals::DataLocation Location, Node& rNode)
{
    if (Location == Globals::DataLocation::NodeHistorical) {
        mProxy = EntityProxy<Globals::DataLocation::NodeHistorical,true>(rNode);
    } else if (Location == Globals::DataLocation::NodeNonHistorical) {
        mProxy = EntityProxy<Globals::DataLocation::NodeNonHistorical,true>(rNode);
    } else {
        KRATOS_ERROR
        << "Constructing a nodal proxy requires 'NodeHistorical' or 'NodeNonHistorical' DataLocation, not '" << GetDataLocationEnumName(Location) << "'";
    }
}


DynamicEntityProxy::DynamicEntityProxy(Globals::DataLocation Location, Element& rElement)
    : mProxy(EntityProxy<Globals::DataLocation::Element,true>(rElement))
{
    KRATOS_ERROR_IF_NOT(Location == Globals::DataLocation::Element)
    << "Constructing an element proxy requires 'Element' DataLocation, not '" << GetDataLocationEnumName(Location) << "'";
}


DynamicEntityProxy::DynamicEntityProxy(Globals::DataLocation Location, Condition& Condition)
    : mProxy(EntityProxy<Globals::DataLocation::Condition,true>(Condition))
{
    KRATOS_ERROR_IF_NOT(Location == Globals::DataLocation::Condition)
    << "Constructing a condition proxy requires 'Condition' DataLocation, not '" << GetDataLocationEnumName(Location) << "'";
}


} // namespace Kratos
