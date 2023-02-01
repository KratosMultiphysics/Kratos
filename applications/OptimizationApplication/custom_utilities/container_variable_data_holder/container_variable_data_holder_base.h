//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class HistoricalContainerDataIO
{
public:
    ///@name Public operations
    ///@{

    template<class TDataType>
    static TDataType& GetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable)
    {
        return rNode.FastGetSolutionStepValue(rVariable);
    }

    template<class TDataType>
    static void SetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.FastGetSolutionStepValue(rVariable) = rValue;
    }

    ///@}
};

class NonHistoricalContainerDataIO
{
public:
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.SetValue(rVariable, rValue);
    }

    ///@}
};

class PropertiesContainerDataIO
{
public:
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetProperties().GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.GetProperties().SetValue(rVariable, rValue);
    }

    ///@}
};

template <class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolderBase {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableDataHolderBase);

    ///@}
    ///@name Public operations
    ///@{

    void CopyDataFrom(const ContainerVariableDataHolderBase<TContainerType>& rOther);

    ///@}
    ///@name Input and output
    ///@{

    IndexType GetDataDimension() const;

    Vector& GetData();

    const Vector& GetData() const;

    ModelPart& GetModelPart();

    const ModelPart& GetModelPart() const;

    TContainerType& GetContainer();

    const TContainerType& GetContainer() const;

    virtual std::string Info() const;

    ///@}
protected:
    ///@name Life cycle
    ///@{

    /// Constructor with the model part
    ContainerVariableDataHolderBase(ModelPart& rModelPart);

    /// Copy constructor
    ContainerVariableDataHolderBase(const ContainerVariableDataHolderBase& rOther);

    ///@}
    ///@name Protected member variables
    ///@{

    Vector mData;

    IndexType mDataDimension;

    ModelPart* const mpModelPart;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableDataHolderBase<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

} // namespace Kratos