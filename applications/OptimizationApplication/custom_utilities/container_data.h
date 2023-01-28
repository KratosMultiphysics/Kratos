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
#include <variant>
#include <functional>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerData {
public:
    ///@name  Enum's
    ///@{

    enum ContainerDataType {
        NodalHistorical        = 0,
        NodalNonHistorical     = 1,
        ConditionNonHistorical = 2,
        ConditionProperties    = 3,
        ElementNonHistorical   = 4,
        ElementProperties      = 5
    };

    ///@}
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition of ContainerData
    KRATOS_CLASS_POINTER_DEFINITION(ContainerData);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ContainerData(
        ModelPart& rModelPart,
        const ContainerDataType& rContainerDataType);

    /// Copy constructor
    ContainerData(const ContainerData& rOther);

    /// Destructor.
    ~ContainerData() = default;

    ///@}
    ///@name Operations
    ///@{

    template<class TDataType>
    void AssignDataToContainerVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void ReadDataFromContainerVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue);

    ContainerData Clone();

    bool IsSameContainer(const ContainerData& rOther) const;

    double NormInf() const;

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

    Vector& GetData();

    const Vector& GetData() const;

    const ContainerDataType& GetContainerDataType() const;

    const ModelPart& GetModelPart() const;

    ModelPart& GetModelPart();

    std::variant<
        std::reference_wrapper<ModelPart::NodesContainerType>,
        std::reference_wrapper<ModelPart::ConditionsContainerType>,
        std::reference_wrapper<ModelPart::ElementsContainerType>> GetContainer();

    std::variant<
        std::reference_wrapper<const ModelPart::NodesContainerType>,
        std::reference_wrapper<const ModelPart::ConditionsContainerType>,
        std::reference_wrapper<const ModelPart::ElementsContainerType>> GetContainer() const;

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Private members
    ///@{

    ModelPart& mrModelPart;
    Vector mData;
    const ContainerDataType& mrContainerDataType;

    ///@}
};

///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const ContainerData& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

///@}
} // namespace Kratos