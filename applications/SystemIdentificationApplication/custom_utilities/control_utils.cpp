//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <algorithm>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes

// Include base h
#include "control_utils.h"

namespace Kratos {

template<class TContainerType>
void ControlUtils::AssignEquivalentProperties(
        TContainerType& rSourceContainer,
        TContainerType& rDestinationContainer)
{
    KRATOS_TRY

    // TODO: To be removed once the PointerVectorSet mutable find is fixed.
    rSourceContainer.Sort();
    rDestinationContainer.Sort();

    const IndexType number_of_entities = rSourceContainer.size();

    KRATOS_ERROR_IF_NOT(number_of_entities == rDestinationContainer.size())
        << "Number of entities mismatch [ rSourceContainer.size() = " << number_of_entities
        << ", rDestinationContainer.size() = " << rDestinationContainer.size() << " ].\n";

    IndexPartition<IndexType>(number_of_entities).for_each([&rSourceContainer, &rDestinationContainer](const auto Index) {
        auto& r_destination_entity = *(rDestinationContainer.begin() + Index);
        auto p_itr = rSourceContainer.find(r_destination_entity.Id());

        KRATOS_ERROR_IF(p_itr == rSourceContainer.end()) <<
            "The entity with id = " << r_destination_entity.Id() << " not found in the source container.";

        r_destination_entity.SetProperties(p_itr->pGetProperties());
    });

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::AssignEquivalentProperties(ModelPart::ConditionsContainerType&, ModelPart::ConditionsContainerType&);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::AssignEquivalentProperties(ModelPart::ElementsContainerType&, ModelPart::ElementsContainerType&);

} /* namespace Kratos.*/