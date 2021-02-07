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

#if !defined(KRATOS_FLUID_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_FLUID_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes

namespace Kratos
{
///@name Classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using PropertiesType = ModelPart::PropertiesType;

    ///@}
    ///@name Static operations
    ///@{

    template<class TDataType>
    static void AssignRandomValues(
        TDataType& rValue,
        const std::string& rSeed,
        const int DomainSize,
        const double MinValue = 0.0,
        const double MaxValue = 1.0);

    static ModelPart& CreateTestModelPart(
        Model& rModel,
        const std::string& rModelPartName,
        const std::string& rElementName,
        const std::string& rConditionName,
        const std::function<void(PropertiesType&)>& rSetElementProperties,
        const std::function<void(PropertiesType&)>& rSetConditionProperties,
        const std::function<void(ModelPart&)>& rAddNodalSolutionStepVariablesFuncion,
        const std::function<void(NodeType&)>& rAddDofsFunction,
        const int BufferSize = 2);

    template<class TDataType>
    static void RandomFillNodalHistoricalVariable(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const double MinValue = 0.0,
        const double MaxValue = 1.0,
        const int Step = 0);

    template<class TContainerType, class TDataType>
    static void RandomFillContainerNonHistoricalVariable(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize,
        const double MinValue = 0.0,
        const double MaxValue = 1.0);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_TEST_UTILITIES_H_INCLUDED