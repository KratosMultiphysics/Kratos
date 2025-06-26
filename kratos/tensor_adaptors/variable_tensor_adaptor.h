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

#pragma once

// System includes
#include <string>
#include <variant>

// External includes

// Project includes
#include "tensor_adaptor.h"
#include "utilities/container_io_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the tensor adaptor types.
 */
class KRATOS_API(KRATOS_CORE) VariableTensorAdaptor: public TensorAdaptor {
private:
    template <class TArgs, class... TExtensionArgs> struct variant_concat;

    template <class... TArgs, class... TExtensionArgs>
    struct variant_concat<std::variant<TArgs...>, TExtensionArgs...> {
        using type = std::variant<TArgs..., TExtensionArgs...>;
    };


    template<class... TDataTypes>
    struct VariableInfo
    {
        using VariableType = std::variant<Variable<TDataTypes> const *...>;
        using __io_historical = std::variant<HistoricalIO<TDataTypes> const *...>;
        using __io_non_historical = typename variant_concat<__io_historical, NonHistoricalIO<TDataTypes> const * ...>::type;
        using ContainerIOType = typename variant_concat<__io_non_historical, GaussPointIO<TDataTypes> const * ...>::type;
    };

public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<VariableTensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const VariableTensorAdaptor>;

    using BaseType = TensorAdaptor;

    using ContainerType = typename BaseType::ContainerType;

    using VariableInfoType = VariableInfo<
                                    int,
                                    double,
                                    array_1d<double, 3>,
                                    array_1d<double, 4>,
                                    array_1d<double, 6>,
                                    array_1d<double, 9>,
                                    Vector,
                                    Matrix
                                >;

    using VariableType = typename VariableInfoType::VariableType;

    using ContainerIOType = typename VariableInfoType::ContainerIOType;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Variable Tensor Adaptor for nodal historical data value container.
     *
     * @param pContainer        Nodes container
     * @param pVariable         Variable
     * @param StepIndex         Step index of the historical data value container.
     */
    VariableTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pContainer,
        VariableType pVariable,
        const int StepIndex);

    /**
     * @brief Construct a new Variable Tensor Adaptor for nodal historical data value container.
     *
     * @param pContainer        Nodes container
     * @param pVariable         Variable
     * @param StepIndex         Step index of the historical data value container.
     * @param rShape            Shape of the data to be stored.
     */
    VariableTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pContainer,
        VariableType pVariable,
        const int StepIndex,
        const std::vector<int>& rShape);

    /**
     * @brief Construct a new Variable Tensor Adaptor object for non-historical data value container in Nodes, Conditions, ...
     *
     * @tparam TContainerPointerType        Pointer type of the container.
     * @param pContainer                    Pointer to the container.
     * @param pVariable                     Variable
     */
    template<class TContainerPointerType>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariableType pVariable);

    /**
     * @brief Construct a new Variable Tensor Adaptor object for non-historical data value container in Nodes, Conditions, ...
     *
     * @tparam TContainerPointerType        Pointer type of the container.
     * @param pContainer                    Pointer to the container.
     * @param pVariable                     Variable
     * @param rShape                        Shape of the data to be stored.
     */
    template<class TContainerPointerType>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariableType pVariable,
        const std::vector<int>& rShape);

    /**
     * @brief Construct a new Variable Tensor Adaptor object for calculating gauss point values in Conditions and Elements
     *
     * @tparam TContainerPointerType        Pointer type of the container.
     * @param pContainer                    Pointer to the container.
     * @param pVariable                     Variable
     * @param rProcessInfo                  Process info
     */
    template<class TContainerPointerType>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariableType pVariable,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Construct a new Variable Tensor Adaptor object for calculating gauss point values in Conditions and Elements
     *
     * @tparam TContainerPointerType        Pointer type of the container.
     * @param pContainer                    Pointer to the container.
     * @param pVariable                     Variable
     * @param rProcessInfo                  Process info
     * @param rShape                        Shape of the data to be stored.
     */
    template<class TContainerPointerType>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariableType pVariable,
        const ProcessInfo& rProcessInfo,
        const std::vector<int>& rShape);

    // Destructor
    ~VariableTensorAdaptor() override;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data from Kratos data structures
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given TContainerType container.
     *
     */
    void StoreData() override;

    ContainerType GetContainer() const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ContainerIOType mpContainerIO;

    ContainerType mpContainer;

    ///@}
};

/// @}
} // namespace Kratos