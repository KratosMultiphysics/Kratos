//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#pragma once

// System includes
#include <type_traits>
#include <variant>

// Project includes
#include "includes/model_part.h"
#include "containers/variable.h"
#include "includes/global_variables.h"
#include "expression/expression.h"
#include "expression/variable_expression_data_io.h"

namespace Kratos
{

namespace DataContainers
{

using IndexType = std::size_t;

using DataLocation = Globals::DataLocation;

template<class TDataType, DataLocation TDataLocation>
class ModelPartVariableDataContainer
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using DataType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    ModelPartVariableDataContainer(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable)
        : mrMesh(rModelPart.GetCommunicator().LocalMesh()),
          mrVariable(rVariable)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    IndexType Size() const
    {
        if constexpr(TDataLocation == DataLocation::NodeHistorical) {
            return mrMesh.NumberOfNodes();
        } else if constexpr(TDataLocation == DataLocation::NodeNonHistorical) {
            return mrMesh.NumberOfNodes();
        } else if constexpr(TDataLocation == DataLocation::Condition) {
            return mrMesh.NumberOfConditions();
        } else if constexpr(TDataLocation == DataLocation::Element) {
            return mrMesh.NumberOfElements();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data location provided.");
        }

        return 0;
    }

    void GetValue(
        TDataType& rOutputValue,
        const IndexType EntityIndex) const
    {
        if constexpr(TDataLocation == DataLocation::NodeHistorical) {
            rOutputValue = (mrMesh.NodesBegin() + EntityIndex)->FastGetSolutionStepValue(mrVariable);
        } else if constexpr(TDataLocation == DataLocation::NodeNonHistorical) {
            rOutputValue = (mrMesh.NodesBegin() + EntityIndex)->GetValue(mrVariable);
        } else if constexpr(TDataLocation == DataLocation::Condition) {
            rOutputValue = (mrMesh.ConditionsBegin() + EntityIndex)->GetValue(mrVariable);
        } else if constexpr(TDataLocation == DataLocation::Element) {
            rOutputValue = (mrMesh.ElementsBegin() + EntityIndex)->GetValue(mrVariable);
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data location provided.");
        }
    }

    unsigned int GetId(const IndexType EntityIndex) const
    {
        if constexpr(TDataLocation == DataLocation::NodeHistorical || TDataLocation == DataLocation::NodeNonHistorical) {
            return (mrMesh.NodesBegin() + EntityIndex)->GetId();
        } else if constexpr(TDataLocation == DataLocation::Condition) {
            return (mrMesh.ConditionsBegin() + EntityIndex)->GetId();
        } else if constexpr(TDataLocation == DataLocation::Element) {
            return (mrMesh.ElementsBegin() + EntityIndex)->GetId();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data location provided.");
            return 0;
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const ModelPart::MeshType& mrMesh;

    const Variable<TDataType>& mrVariable;

    ///@}
};

template<DataLocation TDataLocation>
class ModelPartFlagContainer
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using DataType = unsigned int;

    ///@}
    ///@name Life cycle
    ///@{

    ModelPartFlagContainer(
        const ModelPart& rModelPart,
        const Flags& rFlag)
        : mrMesh(rModelPart.GetCommunicator().LocalMesh()),
          mrFlag(rFlag)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    IndexType Size() const
    {
        if constexpr(TDataLocation == DataLocation::NodeHistorical) {
            return mrMesh.NumberOfNodes();
        } else if constexpr(TDataLocation == DataLocation::NodeNonHistorical) {
            return mrMesh.NumberOfNodes();
        } else if constexpr(TDataLocation == DataLocation::Condition) {
            return mrMesh.NumberOfConditions();
        } else if constexpr(TDataLocation == DataLocation::Element) {
            return mrMesh.NumberOfElements();
        } else {
            static_assert(TDataLocation != TDataLocation, "Unsupported data location provided.");
        }

        return 0;
    }

    void GetValue(
        DataType& rOutputValue,
        const IndexType EntityIndex) const
    {
        if constexpr(TDataLocation == DataLocation::NodeHistorical) {
            rOutputValue = (mrMesh.NodesBegin() + EntityIndex)->Is(mrFlag);
        } else if constexpr(TDataLocation == DataLocation::NodeNonHistorical) {
            rOutputValue = (mrMesh.NodesBegin() + EntityIndex)->Is(mrFlag);
        } else if constexpr(TDataLocation == DataLocation::Condition) {
            rOutputValue = (mrMesh.ConditionsBegin() + EntityIndex)->Is(mrFlag);
        } else if constexpr(TDataLocation == DataLocation::Element) {
            rOutputValue = (mrMesh.ElementsBegin() + EntityIndex)->Is(mrFlag);
        } else {
            static_assert(TDataLocation != TDataLocation, "Unsupported data location provided.");
        }
    }

    unsigned int GetId(const IndexType EntityIndex) const
    {
        if constexpr(TDataLocation == DataLocation::NodeHistorical || TDataLocation == DataLocation::NodeNonHistorical) {
            return (mrMesh.NodesBegin() + EntityIndex)->GetId();
        } else if constexpr(TDataLocation == DataLocation::Condition) {
            return (mrMesh.ConditionsBegin() + EntityIndex)->GetId();
        } else if constexpr(TDataLocation == DataLocation::Element) {
            return (mrMesh.ElementsBegin() + EntityIndex)->GetId();
        } else {
            static_assert(TDataLocation != TDataLocation, "Unsupported data location provided.");
            return 0;
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const ModelPart::MeshType& mrMesh;

    const Flags& mrFlag;

    ///@}
};

template<class TDataType>
class ExpressionDataContainer
{
public:
    ///@name Type definitions
    ///@{

    using DataType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    ExpressionDataContainer(const Expression& rExpression)
        : mrExpression(rExpression),
          mVariableExpressionDataIO(rExpression.GetItemShape())
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    IndexType Size() const { return mrExpression.NumberOfEntities(); }

    void GetValue(
        TDataType& rOutputValue,
        const IndexType EntityIndex) const
    {
        mVariableExpressionDataIO.Assign(rOutputValue, mrExpression, EntityIndex);
    }

    unsigned int GetId(const IndexType EntityIndex) const
    {
        return EntityIndex;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Expression& mrExpression;

    const VariableExpressionDataIO<TDataType> mVariableExpressionDataIO;

    ///@}
};

template<class TDataType>
std::variant<
    ModelPartVariableDataContainer<TDataType, DataLocation::NodeHistorical>,
    ModelPartVariableDataContainer<TDataType, DataLocation::NodeNonHistorical>,
    ModelPartVariableDataContainer<TDataType, DataLocation::Condition>,
    ModelPartVariableDataContainer<TDataType, DataLocation::Element>> GetDataContainer(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rDataLocation)
{
    switch (rDataLocation) {
        case DataLocation::NodeHistorical:
            return ModelPartVariableDataContainer<TDataType, DataLocation::NodeHistorical>(rModelPart, rVariable);
        case DataLocation::NodeNonHistorical:
            return ModelPartVariableDataContainer<TDataType, DataLocation::NodeNonHistorical>(rModelPart, rVariable);
        case DataLocation::Condition:
            return ModelPartVariableDataContainer<TDataType, DataLocation::Condition>(rModelPart, rVariable);
        case DataLocation::Element:
            return ModelPartVariableDataContainer<TDataType, DataLocation::Element>(rModelPart, rVariable);
        default:
            KRATOS_ERROR << "Only supports following data containers of the model parts:"
                         << "\n\t" << "NodeHistorical"
                         << "\n\t" << "NodeNonHistorical"
                         << "\n\t" << "Condition"
                         << "\n\t" << "Element";
            return ModelPartVariableDataContainer<TDataType, DataLocation::NodeHistorical>(rModelPart, rVariable);
    }
}

std::variant<
    ModelPartFlagContainer<DataLocation::NodeHistorical>,
    ModelPartFlagContainer<DataLocation::NodeNonHistorical>,
    ModelPartFlagContainer<DataLocation::Condition>,
    ModelPartFlagContainer<DataLocation::Element>> GetDataContainer(
    const ModelPart& rModelPart,
    const Flags& rFlag,
    const DataLocation& rDataLocation)
{
    switch (rDataLocation) {
        case DataLocation::NodeHistorical:
            return ModelPartFlagContainer<DataLocation::NodeHistorical>(rModelPart, rFlag);
        case DataLocation::NodeNonHistorical:
            return ModelPartFlagContainer<DataLocation::NodeNonHistorical>(rModelPart, rFlag);
        case DataLocation::Condition:
            return ModelPartFlagContainer<DataLocation::Condition>(rModelPart, rFlag);
        case DataLocation::Element:
            return ModelPartFlagContainer<DataLocation::Element>(rModelPart, rFlag);
        default:
            KRATOS_ERROR << "Only supports following data containers of the model parts:"
                         << "\n\t" << "NodeHistorical"
                         << "\n\t" << "NodeNonHistorical"
                         << "\n\t" << "Condition"
                         << "\n\t" << "Element";
            return ModelPartFlagContainer<DataLocation::NodeHistorical>(rModelPart, rFlag);
    }
}

std::variant<
    ExpressionDataContainer<double>,
    ExpressionDataContainer<array_1d<double, 3>>,
    ExpressionDataContainer<array_1d<double, 4>>,
    ExpressionDataContainer<array_1d<double, 6>>,
    ExpressionDataContainer<array_1d<double, 9>>,
    ExpressionDataContainer<Vector>,
    ExpressionDataContainer<Matrix>> GetDataContainer(const Expression& rExpression)
{
    const auto& r_shape = rExpression.GetItemShape();
    if (r_shape.size() == 0) {
        return ExpressionDataContainer<double>(rExpression);
    } else if (r_shape.size() == 1) {
        switch (r_shape[0]) {
            case 3:
                return ExpressionDataContainer<array_1d<double, 3>>(rExpression);
            case 4:
                return ExpressionDataContainer<array_1d<double, 4>>(rExpression);
            case 6:
                return ExpressionDataContainer<array_1d<double, 6>>(rExpression);
            case 9:
                return ExpressionDataContainer<array_1d<double, 9>>(rExpression);
            default:
                return ExpressionDataContainer<Vector>(rExpression);
        }
    } else if (r_shape.size() == 2) {
        return ExpressionDataContainer<Matrix>(rExpression);
    } else {
        KRATOS_ERROR << "Unsupported expression.";
        return ExpressionDataContainer<double>(rExpression);
    }
}

} // namespace DataContainers

} // namespace Kratos
