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
#include "includes/global_variables.h" // Globals::DataLocation
#include "includes/model_part.h" // ModelPart::NodesContainerType, ModelPart::ElementsContainerType, ModelPart::ConditionsContainerType
#include "utilities/proxies.h" // EntityProxy, ContainerProxy
#include "includes/kratos_export_api.h" // KRATOS_API

// System includes
#include <variant> // variant


namespace Kratos {


class KRATOS_API(KRATOS_CORE) DynamicEntityProxy
{
public:
    DynamicEntityProxy() noexcept = default;

    template <Globals::DataLocation TLocation>
    DynamicEntityProxy(EntityProxy<TLocation,true> Proxy) noexcept : mProxy(Proxy) {}

    DynamicEntityProxy(Globals::DataLocation Location, Node& rNode);

    DynamicEntityProxy(Globals::DataLocation Location, Element& rElement);

    DynamicEntityProxy(Globals::DataLocation Location, Condition& Condition);

    template <class TVariable>
    bool HasValue(const TVariable& rVariable) const
    {
        KRATOS_TRY
        return std::visit(
            [&rVariable](auto Proxy){
                return Proxy.HasValue(rVariable);
            },
            mProxy
        );
        KRATOS_CATCH("")
    }

    template <class TVariable>
    const typename TVariable::Type& GetValue(const TVariable& rVariable) const
    {
        KRATOS_TRY
        return std::visit(
            [&rVariable](auto Proxy){
                return Proxy.GetValue(rVariable);
            },
            mProxy
        );
        KRATOS_CATCH("")
    }

    template <class TVariable>
    typename TVariable::Type& GetValue(const TVariable& rVariable)
    {
        KRATOS_TRY
        return std::visit(
            [&rVariable](auto Proxy){
                return Proxy.GetValue(rVariable);
            },
            mProxy
        );
        KRATOS_CATCH("")
    }

    template <class TVariable>
    void SetValue(const TVariable& rVariable, const typename TVariable::Type& rValue) const
    {
        KRATOS_TRY
        return std::visit(
            [&rVariable, &rValue](auto Proxy){
                return Proxy.SetValue(rVariable, rValue);
            },
            mProxy
        );
        KRATOS_CATCH("")
    }

private:
    std::variant<
        EntityProxy<Globals::DataLocation::NodeHistorical,true>,
        EntityProxy<Globals::DataLocation::NodeNonHistorical,true>,
        EntityProxy<Globals::DataLocation::Element,true>,
        EntityProxy<Globals::DataLocation::Condition,true>
    > mProxy;
}; // class DynamicEntityProxy


} // namespace Kratos
