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

class HistoricalDataValueContainer
{
public:
    using ContainerType = ModelPart::NodesContainerType;

    template<class TDataType>
    static TDataType& GetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable);

    template<class TDataType>
    static void SetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable, const TDataType& rValue);
};

template<class TContainerType>
class NonHistoricalDataValueContainer
{
public:
    using ContainerType = TContainerType;

    template<class TDataType>
    static TDataType& GetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable);

    template<class TDataType>
    static void SetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable, const TDataType& rValue);
};

template<class TContainerType>
class PropertiesDataValueContainer
{
public:
    using ContainerType = TContainerType;

    template<class TDataType>
    static TDataType& GetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable);

    template<class TDataType>
    static void SetValue(typename ContainerType::data_type& rEntity, const Variable<TDataType>& rVariable, const TDataType& rValue);
};

template<class TContainerDataType>
class ContainerData;

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
    ///@name Public operations
    ///@{

    double NormInf() const;

    ///@}
    ///@name Input and output
    ///@{

    IndexType GetDataDimension() const;

    Vector& GetData();

    const Vector& GetData() const;

    void CopyData(const ContainerDataBase& rOther);

    ModelPart& GetModelPart();

    const ModelPart& GetModelPart() const;

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    ///@}
private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart;

    Vector mData;

    IndexType mDataDimension;

    ///@}
    ///@name Private life cycle
    ///@{

    /// Constructor
    ContainerDataBase(ModelPart& rModelPart);

    /// Copy constructor
    ContainerDataBase(const ContainerDataBase& rOther);

    ///@}
    ///@name Friend classes
    ///@{

    template<class TContainerDataType>
    friend class ContainerData;

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
    ContainerData(ModelPart& rModelPart) : BaseType(rModelPart) {}

    /// Copy constructor
    ContainerData(const ContainerData& rOther) : BaseType(rOther) {}

    /// Destructor.
    virtual ~ContainerData() override = default;

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

    double InnerProduct(const ContainerData& rOther) const;

    ///@}
    ///@name Operators
    ///@{

    ContainerData operator+(const ContainerData& rOther) const;

    ContainerData& operator+=(const ContainerData& rOther);

    ContainerData operator-(const ContainerData& rOther) const;

    ContainerData& operator-=(const ContainerData& rOther);

    ContainerData operator*(const double Value) const;

    ContainerData& operator*=(const double Value);

    ContainerData operator/(const double Value) const;

    ContainerData& operator/=(const double Value);

    ContainerData& operator=(const ContainerData& rOther);

    ///@}
    ///@name Input and output
    ///@{

    typename TContainerDataType::ContainerType& GetContainer();

    const typename TContainerDataType::ContainerType& GetContainer() const;

    std::string Info() const override;

    ///@}
};

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