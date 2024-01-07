//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(FLUID_MODEL_PART_PREPROCESSING_UTILITIES)
#define FLUID_MODEL_PART_PREPROCESSING_UTILITIES

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidModelPartPreProcessingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public Static Operations
    ///@{

    bool static CreateModelPartForCommenInterface(
        ModelPart& rModelPart,
        const std::string& rCommonInterfaceModelPartName,
        const std::vector<std::string>& rListOfInterfaceModelPartNames
    );

    std::vector<IndexType> static GetElementIdsWithAllNodesOnBoundaries(
        ModelPart& rModelPart,
        const std::vector<std::string>& rListOfBoundaryModelPartNames
    );

    void static BreakElements(
        ModelPart& rModelPart,
        const std::string& rNewElementName,
        const std::vector<IndexType>& rElementIds
    );

    ///@}

private:
    ///@name Private member functions
    ///@{

    IndexType static BreakElement(
        ModelPart& rModelPart,
        ModelPart::ElementType& rElement,
        IndexType& NewNodeId,
        IndexType& NewElementId,
        const std::string& rNewElementName
    );

    ///@}
};

} // namespace Kratos

#endif // FLUID_MODEL_PART_PREPROCESSING_UTILITIES
