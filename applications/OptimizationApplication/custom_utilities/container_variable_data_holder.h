//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
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

///@name Forward declarations
///@{

class HistoricalDataValueContainer;

template<class TContainerType>
class NonHistoricalDataValueContainer;

template<class TContainerType>
class PropertiesDataValueContainer;

template<class TContainerVariableDataHolderType>
class ContainerVariableDataHolder;

///@}
///@name Class declarations
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolderBase {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition of ContainerVariableDataHolderBase
    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableDataHolderBase);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~ContainerVariableDataHolderBase() = default;

    ///@}
    ///@name Input and output
    ///@{

    IndexType GetDataDimension() const;

    Vector& GetData();

    const Vector& GetData() const;

    void CopyDataFrom(const ContainerVariableDataHolderBase& rOther);

    ModelPart& GetModelPart();

    const ModelPart& GetModelPart() const;

    bool IsCompatibleWithContainerVariableDataHolder(const ContainerVariableDataHolderBase& rOther) const;

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    ///@}
protected:
    ///@name Enums
    ///@{

    enum ContainerVariableDataHolderType
    {
        HistoricalContainerVariableDataHolder          = 0,
        NodalContainerVariableDataHolder               = 1,
        ConditionContainerVariableDataHolder           = 2,
        ElementContainerVariableDataHolder             = 3,
        ConditionPropertiesContainerVariableDataHolder = 4,
        ElementPropertiesContainerVariableDataHolder   = 5
    };

    ///@}
    ///@name Protected life cycle
    ///@{

    /// Constructor
    ContainerVariableDataHolderBase(
        ModelPart& rModelPart,
        const ContainerVariableDataHolderType& rContainerVariableDataHolderType);

    /// Copy constructor
    ContainerVariableDataHolderBase(const ContainerVariableDataHolderBase& rOther);

    ///@}
    ///@name Protected member variables
    ///@{

    ModelPart& mrModelPart;

    Vector mData;

    IndexType mDataDimension;

    ///@}
private:
    ///@name Private member variables
    ///@{

    const ContainerVariableDataHolderType mContainerVariableDataHolderType;

    ///@}
    ///@name Friend classes
    ///@{

    friend class HistoricalDataValueContainer;

    template<class TContainerType>
    friend class NonHistoricalDataValueContainer;

    template<class TContainerType>
    friend class PropertiesDataValueContainer;

    ///@}
};

class HistoricalDataValueContainer
{
private:
    ///@name Private type definitions
    ///@{

    using ContainerType = ModelPart::NodesContainerType;

    static constexpr ContainerVariableDataHolderBase::ContainerVariableDataHolderType ContainerVariableDataHolderType
        = ContainerVariableDataHolderBase::ContainerVariableDataHolderType::HistoricalContainerVariableDataHolder;

    ///@}
    ///@name Private static operations
    ///@{

    template<class TDataType>
    static TDataType& GetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable);

    template<class TDataType>
    static void SetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable, const TDataType& rValue);

    ///@}
    ///@name Friend classes
    ///@{

    friend class ContainerVariableDataHolder<HistoricalDataValueContainer>;

    ///@}
};

template<class TContainerType>
class NonHistoricalDataValueContainer
{
private:
    ///@name Private type definitions
    ///@{

    using ContainerType = TContainerType;

    static constexpr ContainerVariableDataHolderBase::ContainerVariableDataHolderType ContainerVariableDataHolderType
        = std::is_same_v<TContainerType, ModelPart::NodesContainerType>
            ? ContainerVariableDataHolderBase::ContainerVariableDataHolderType::NodalContainerVariableDataHolder
            : std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                ? ContainerVariableDataHolderBase::ContainerVariableDataHolderType::ConditionContainerVariableDataHolder
                : ContainerVariableDataHolderBase::ContainerVariableDataHolderType::ElementContainerVariableDataHolder;

    ///@}
    ///@name Private static operations
    ///@{

    template<class TDataType>
    static TDataType& GetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable);

    template<class TDataType>
    static void SetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable, const TDataType& rValue);

    ///@}
    ///@name Friend classes
    ///@{

    friend class ContainerVariableDataHolder<NonHistoricalDataValueContainer<TContainerType>>;

    ///@}
};

template<class TContainerType>
class PropertiesDataValueContainer
{
private:
    ///@name Private type definitions
    ///@{

    using ContainerType = TContainerType;

    static constexpr ContainerVariableDataHolderBase::ContainerVariableDataHolderType ContainerVariableDataHolderType
        = std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
            ? ContainerVariableDataHolderBase::ContainerVariableDataHolderType::ConditionPropertiesContainerVariableDataHolder
            : ContainerVariableDataHolderBase::ContainerVariableDataHolderType::ElementPropertiesContainerVariableDataHolder;

    ///@}
    ///@name Private static operations
    ///@{

    template<class TDataType>
    static TDataType& GetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable);

    template<class TDataType>
    static void SetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable, const TDataType& rValue);

    ///@}
    ///@name Friend classes
    ///@{

    friend class ContainerVariableDataHolder<PropertiesDataValueContainer<TContainerType>>;

    ///@}
};

template<class TContainerVariableDataHolderType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolder : public ContainerVariableDataHolderBase
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerVariableDataHolderBase;

    using IndexType = std::size_t;

    /// Pointer definition of ContainerVariableDataHolder
    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableDataHolder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ContainerVariableDataHolder(ModelPart& rModelPart) : BaseType(rModelPart, TContainerVariableDataHolderType::ContainerVariableDataHolderType) {}

    /// Copy constructor
    ContainerVariableDataHolder(const ContainerVariableDataHolder& rOther) : BaseType(rOther) {}

    /// Destructor.
    ~ContainerVariableDataHolder() override = default;

    ///@}
    ///@name Public operations
    ///@{

    ContainerVariableDataHolder<TContainerVariableDataHolderType> Clone();

    template<class TDataType>
    void ReadDataFromContainerVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void AssignDataToContainerVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue);

    ///@}
    ///@name Operators
    ///@{

    ContainerVariableDataHolder operator+(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator+=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator+(const double Value) const;

    ContainerVariableDataHolder& operator+=(const double Value);

    ContainerVariableDataHolder operator-(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator-=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator-(const double Value) const;

    ContainerVariableDataHolder& operator-=(const double Value);

    ContainerVariableDataHolder operator*(const double Value) const;

    ContainerVariableDataHolder& operator*=(const double Value);

    ContainerVariableDataHolder operator/(const double Value) const;

    ContainerVariableDataHolder& operator/=(const double Value);

    ContainerVariableDataHolder operator^(const double Value) const;

    ContainerVariableDataHolder& operator^=(const double Value);

    ContainerVariableDataHolder& operator=(const ContainerVariableDataHolder& rOther);

    ///@}
    ///@name Input and output
    ///@{

    typename TContainerVariableDataHolderType::ContainerType& GetContainer();

    const typename TContainerVariableDataHolderType::ContainerType& GetContainer() const;

    std::string Info() const override;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableDataHolderBase& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

template<class TContainerVariableDataHolderType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableDataHolder<TContainerVariableDataHolderType>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

///@}
} // namespace Kratos