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

// System includes
#include <type_traits>

// Project includes
#include "includes/kratos_export_api.h"
#include "includes/model_part.h"
#include "includes/global_variables.h"

namespace Kratos {


/// @brief Class collecting a set of free-floating utility functions for querying and mutating @ref ModelPart s.
class ModelPartUtils
{
public:
    /// @brief Templated interface for getting nodes, elements, conditions or @ref ProcessInfo from a @ref ModelPart.
    template <Globals::DataLocation TLocation>
    static const auto& GetContainer(const ModelPart& rModelPart)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical) {
            return rModelPart.Nodes();
        } else if constexpr (TLocation == Globals::DataLocation::Element) {
            return rModelPart.Elements();
        } else if constexpr (TLocation == Globals::DataLocation::Condition) {
            return rModelPart.Conditions();
        } else if constexpr (TLocation == Globals::DataLocation::ProcessInfo) {
            return rModelPart.GetProcessInfo();
        } else if constexpr (TLocation == Globals::DataLocation::ModelPart) {
            return rModelPart;
        }
    }

    /// @brief Templated interface for getting nodes, elements, conditions or @ref ProcessInfo from a @ref ModelPart.
    template <Globals::DataLocation TLocation>
    static auto& GetContainer(ModelPart& rModelPart)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical) {
            return rModelPart.Nodes();
        } else if constexpr (TLocation == Globals::DataLocation::Element) {
            return rModelPart.Elements();
        } else if constexpr (TLocation == Globals::DataLocation::Condition) {
            return rModelPart.Conditions();
        } else if constexpr (TLocation == Globals::DataLocation::ProcessInfo) {
            return rModelPart.GetProcessInfo();
        } else if constexpr (TLocation == Globals::DataLocation::ModelPart) {
            return rModelPart;
        }
    }

    /// @brief Templated interface to get nodes, elements and conditions from a @ref ModelPart
    template<class TContainerType>
    static const auto& GetContainer(const ModelPart& rModelPart)
    {
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            return GetContainer<Globals::DataLocation::NodeNonHistorical>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            return GetContainer<Globals::DataLocation::Condition>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            return GetContainer<Globals::DataLocation::Element>(rModelPart);
        } else {
            static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
            return 0;
        }
    }

    /// @brief Templated interface to get nodes, elements and conditions from a @ref ModelPart
    template<class TContainerType>
    static auto& GetContainer(ModelPart& rModelPart)
    {
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            return GetContainer<Globals::DataLocation::NodeNonHistorical>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            return GetContainer<Globals::DataLocation::Condition>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            return GetContainer<Globals::DataLocation::Element>(rModelPart);
        } else {
            static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
            return 0;
        }
    }

    /// @brief Add nodes to ModelPart from an ordered container.
    /** By assuming that the input is ordered (by increasing Id), the nodes can be added more efficiently.
     *  Note that the function makes no check of the ordering, it is the responsability of the caller to ensure that it is correct.
     *  @tparam TIteratorType Iterator type for the nodes to add.
     *  @param rTargetModelPart ModelPart the nodes will be added to.
     *  @param iNodesBegin First element of the container of nodes to be added to rTargetModelPart.
     *  @param iNodesEnd End position for the container of nodes to be added to rTargetModelPart.
     */
    template<class TIteratorType >
    static void AddNodesFromOrderedContainer(ModelPart& rTargetModelPart, TIteratorType iNodesBegin,  TIteratorType iNodesEnd)
    {
        KRATOS_TRY
        ModelPart::NodesContainerType aux;
        ModelPart* root_model_part = &rTargetModelPart.GetRootModelPart();

        for(TIteratorType it = iNodesBegin; it!=iNodesEnd; it++) {
            aux.push_back( *(it.base()) );
        }
        // Add the nodes to the root modelpart
        root_model_part->Nodes() = JoinOrderedNodesContainerType(root_model_part->Nodes().begin(), root_model_part->Nodes().end(), aux.begin(), aux.end());

        // Add to all of the leaves
        ModelPart* current_part = &rTargetModelPart;
        while(current_part->IsSubModelPart()) {
            current_part->Nodes() = JoinOrderedNodesContainerType(current_part->Nodes().begin(), current_part->Nodes().end(), aux.begin(), aux.end());
            current_part = &(current_part->GetParentModelPart());
        }

        KRATOS_CATCH("")
    }

private:

    /// @brief Helper function for AddNodesFromOrderedContainer
    template<class TIteratorType >
    static typename ModelPart::NodesContainerType JoinOrderedNodesContainerType(
        TIteratorType iC1Begin,  TIteratorType iC1End, TIteratorType iC2Begin,  TIteratorType iC2End)
    {
        std::size_t c1length = std::distance(iC1Begin,iC1End);
        std::size_t c2length = std::distance(iC2Begin,iC2End);
        TIteratorType blong, elong, bshort, eshort;
        ModelPart::NodesContainerType aux;
        aux.reserve(c1length + c2length);

        // We order c1 and c2 to long and short
        if(c1length>c2length){
            blong = iC1Begin; elong = iC1End;
            bshort = iC2Begin; eshort = iC2End;
        } else {
            blong = iC2Begin; elong = iC2End;
            bshort = iC1Begin; eshort = iC1End;
        }

        // If short is empty we return long. If both empty it returns empty aux
        if(c2length == 0 || c1length == 0){
            for(TIteratorType it1=blong; it1!=elong; it1++){
                aux.push_back(*(it1.base()) );
            }
            return aux;
        }
        TIteratorType it2 = blong;
        for(TIteratorType it1=bshort; it1!=eshort; it1++) {
            while(it2!=elong && it2->Id()<it1->Id()) {
                aux.push_back(*(it2.base()));
                it2++;
            }
            aux.push_back(*(it1.base()) );
            if(it2!=elong && (it1->Id() == it2->Id())) { //If both are the same, then we need to skip
                it2++;
            }
        }
        while(it2 != elong) {
            aux.push_back(*(it2.base()));
            it2++;
        }

        return aux;
    }


}; // class ModelPartUtils




} // namespace Kratos
