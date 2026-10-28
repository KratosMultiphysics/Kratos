//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

#include "includes/model_part.h"

namespace Kratos {
/** @brief Build additive gap quadrature and connector entities from surrogate faces. */
class KRATOS_API(KRATOS_CORE) GapConstructionUtility
{
public:
    static void AssignParents(ModelPart& rRoot, const std::vector<std::size_t>& rConditions,
        const std::vector<std::size_t>& rElements);
    static void Build(ModelPart& rOutput, ModelPart& rSurrogate,
        ModelPart& rGapElements, ModelPart& rGapInterfaces,
        const std::string& rElementName, const std::string& rInterfaceConditionName,
        double Tolerance, int EchoLevel = 0);
};
} // namespace Kratos
