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
}; // class ModelPartUtils


} // namespace Kratos
