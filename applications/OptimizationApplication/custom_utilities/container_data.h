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

template<class TContainerDataType>
class ContainerData;

///@}
///@name Class declarations
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerDataBase {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition of ContainerDataBase
    KRATOS_CLASS_POINTER_DEFINITION(ContainerDataBase);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~ContainerDataBase() = default;

    ///@}
    ///@name Input and output
    ///@{

    IndexType GetDataDimension() const;

    Vector& GetData();

    const Vector& GetData() const;

    void CopyDataFrom(const ContainerDataBase& rOther);

    ModelPart& GetModelPart();

    const ModelPart& GetModelPart() const;

    bool IsCompatibleWithContainerData(const ContainerDataBase& rOther) const;

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    ///@}
protected:
    ///@name Enums
    ///@{

    enum ContainerDataType
    {
        HistoricalContainerData          = 0,
        NodalContainerData               = 1,
        ConditionContainerData           = 2,
        ElementContainerData             = 3,
        ConditionPropertiesContainerData = 4,
        ElementPropertiesContainerData   = 5
    };

    ///@}
    ///@name Protected life cycle
    ///@{

    /// Constructor
    ContainerDataBase(
        ModelPart& rModelPart,
        const ContainerDataType& rContainerDataType);

    /// Copy constructor
    ContainerDataBase(const ContainerDataBase& rOther);

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

    const ContainerDataType mContainerDataType;

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

    static constexpr ContainerDataBase::ContainerDataType ContainerDataType
        = ContainerDataBase::ContainerDataType::HistoricalContainerData;

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

    friend class ContainerData<HistoricalDataValueContainer>;

    ///@}
};

template<class TContainerType>
class NonHistoricalDataValueContainer
{
private:
    ///@name Private type definitions
    ///@{

    using ContainerType = TContainerType;

    static constexpr ContainerDataBase::ContainerDataType ContainerDataType
        = std::is_same_v<TContainerType, ModelPart::NodesContainerType>
            ? ContainerDataBase::ContainerDataType::NodalContainerData
            : std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                ? ContainerDataBase::ContainerDataType::ConditionContainerData
                : ContainerDataBase::ContainerDataType::ElementContainerData;

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

    friend class ContainerData<NonHistoricalDataValueContainer<TContainerType>>;

    ///@}
};

template<class TContainerType>
class PropertiesDataValueContainer
{
private:
    ///@name Private type definitions
    ///@{

    using ContainerType = TContainerType;

    static constexpr ContainerDataBase::ContainerDataType ContainerDataType
        = std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
            ? ContainerDataBase::ContainerDataType::ConditionPropertiesContainerData
            : ContainerDataBase::ContainerDataType::ElementPropertiesContainerData;

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

    friend class ContainerData<PropertiesDataValueContainer<TContainerType>>;

    ///@}
};

template<class TContainerDataType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerData : public ContainerDataBase
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerDataBase;

    using IndexType = std::size_t;

    /// Pointer definition of ContainerData
    KRATOS_CLASS_POINTER_DEFINITION(ContainerData);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ContainerData(ModelPart& rModelPart) : BaseType(rModelPart, TContainerDataType::ContainerDataType) {}

    /// Copy constructor
    ContainerData(const ContainerData& rOther) : BaseType(rOther) {}

    /// Destructor.
    ~ContainerData() override = default;

    ///@}
    ///@name Public operations
    ///@{

    ContainerData<TContainerDataType> Clone();

    template<class TDataType>
    void ReadDataFromContainerVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void AssignDataToContainerVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue);

    ///@}
    ///@name Operators
    ///@{

    ContainerData operator+(const ContainerData& rOther) const;

    ContainerData& operator+=(const ContainerData& rOther);

    ContainerData operator+(const double Value) const;

    ContainerData& operator+=(const double Value);

    ContainerData operator-(const ContainerData& rOther) const;

    ContainerData& operator-=(const ContainerData& rOther);

    ContainerData operator-(const double Value) const;

    ContainerData& operator-=(const double Value);

    ContainerData operator*(const double Value) const;

    ContainerData& operator*=(const double Value);

    ContainerData operator/(const double Value) const;

    ContainerData& operator/=(const double Value);

    ContainerData operator^(const double Value) const;

    ContainerData& operator^=(const double Value);

    ContainerData& operator=(const ContainerData& rOther);

    ///@}
    ///@name Input and output
    ///@{

    typename TContainerDataType::ContainerType& GetContainer();

    const typename TContainerDataType::ContainerType& GetContainer() const;

    std::string Info() const override;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerDataBase& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

template<class TContainerDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerData<TContainerDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

///@}
} // namespace Kratos