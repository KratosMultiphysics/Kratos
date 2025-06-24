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
#include "tensor_adaptor_utils.h"
#include "utilities/container_io_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the tensor adaptor types.
 */
template<class TContainerType, template<class> class TContainerIOType, class... TArgs>
class VariantVariableTensorAdaptor: public TensorAdaptor<TContainerType, double> {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<VariantVariableTensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const VariantVariableTensorAdaptor>;

    using BaseType = TensorAdaptor<TContainerType, double>;

    using ContainerType = typename BaseType::ContainerType;

    using VariableType = std::variant<
                                    const Variable<double>*,
                                    const Variable<array_1d<double, 3>>*,
                                    const Variable<array_1d<double, 4>>*,
                                    const Variable<array_1d<double, 6>>*,
                                    const Variable<array_1d<double, 9>>*,
                                    const Variable<Vector>*,
                                    const Variable<Matrix>*
                                >;

    using ContainerIOType = std::variant<
                                    TContainerIOType<double> const *,
                                    TContainerIOType<array_1d<double, 3>> const *,
                                    TContainerIOType<array_1d<double, 4>> const *,
                                    TContainerIOType<array_1d<double, 6>> const *,
                                    TContainerIOType<array_1d<double, 9>> const *,
                                    TContainerIOType<Vector> const *,
                                    TContainerIOType<Matrix> const *
                                >;

    ///@}
    ///@name Life cycle
    ///@{

    VariantVariableTensorAdaptor(
        typename ContainerType::Pointer pContainer,
        VariableType pVariable,
        TArgs&&... rArgs)
        : BaseType(pContainer)
    {
        std::visit([&](const auto pVariable){
            using data_type = typename std::remove_cv_t<std::decay_t<decltype(*pVariable)>>::Type;
            auto p_container_io = new TContainerIOType<data_type>(*pVariable, rArgs...);
            TensorAdaptorUtils::GetShape(this->mShape, *(this->mpContainer), *p_container_io);
            this->mpContainerIO = p_container_io;
        }, pVariable);

        this->mData.resize(TensorAdaptorUtils::GetFlatLength(this->mShape.data(), this->mShape.data() + this->mShape.size()));
    }

    VariantVariableTensorAdaptor(
        typename ContainerType::Pointer pContainer,
        VariableType pVariable,
        const std::vector<int>& rShape,
        TArgs&&... rArgs)
        : BaseType(pContainer)
    {
        std::visit([&](const auto pVariable){
            using data_type = typename std::remove_cv_t<std::decay_t<decltype(*pVariable)>>::Type;

            auto p_container_io = new TContainerIOType<data_type>(*pVariable, rArgs...);
            this->mpContainerIO = p_container_io;

            KRATOS_ERROR_IF_NOT(DataTypeTraits<data_type>::Dimension + 1 == rShape.size())
                << "Dimensions mismatch for " << pVariable->Name() << " [ Required dimensions by variable = "
                << DataTypeTraits<data_type>::Dimension + 1 << ", shape = " << rShape
                << ", TensorAdaptor = " << *this << " ].\n";

            if constexpr(!DataTypeTraits<data_type>::IsDynamic) {
                // now we know the shape exactly, hence we can check whether
                // user has provided the shape correctly. These types are
                // double, array3, array4, ...
                TensorAdaptorUtils::GetShape(this->mShape, *(this->mpContainer), *p_container_io);

                for (IndexType i = 0; i < rShape.size(); ++i) {
                    KRATOS_ERROR_IF_NOT(rShape[i] == this->mShape[i])
                                << "Shape mismatch for " << pVariable->Name()
                                << " [ Required variable shape = " << this->mShape << ", shape = " << rShape
                                << ", TensorAdaptor = " << *this << " ].\n";
                }
            }

        }, pVariable);

        KRATOS_ERROR_IF_NOT(rShape[0] == pContainer->size())
            << "First dimension of the shape should represent the container size [ shape = "
            << rShape << ", size of the container = " << pContainer->size() << ", Tensor = " << *this << " ].\n";

        this->mShape = rShape;
        this->mData.resize(TensorAdaptorUtils::GetFlatLength(this->mShape.data(), this->mShape.data() + this->mShape.size()));
    }

    ~VariantVariableTensorAdaptor()
    {
        std::visit([](const auto pContainerIO){ delete pContainerIO; }, this->mpContainerIO);
    }

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data from Kratos data structures
     */
    void CollectData() override
    {
        KRATOS_TRY

        // sanity checks
        KRATOS_ERROR_IF_NOT(this->mShape[0] == static_cast<int>(this->mpContainer->size()))
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << "Tensor adapter shape = " << this->mShape << ", container size = " << this->mpContainer->size()
            << ", TensorAdaptor = " << *this << " ].\n";

        std::visit([this](auto pContainerIO) {
            CopyToContiguousArray(*(this->mpContainer), *pContainerIO, this->mData.data().begin(), this->mData.size());
        }, this->mpContainerIO);

        KRATOS_CATCH("");
    }

    /**
     * @brief Store internal data to the given TContainerType container.
     *
     */
    void StoreData() override
    {
        KRATOS_TRY

        // sanity checks
        KRATOS_ERROR_IF_NOT(this->mShape[0] == static_cast<int>(this->mpContainer->size()))
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << "Tensor adapter shape = " << this->mShape << ", container size = " << this->mpContainer->size()
            << ", TensorAdaptor = " << *this << " ].\n";

        std::visit([this](auto pContainerIO) {
            std::vector<unsigned int> shape;
            shape.resize(this->mShape.size() - 1);
            std::copy(this->mShape.begin() + 1, this->mShape.end(), shape.begin());
            CopyFromContiguousDataArray(*(this->mpContainer), *pContainerIO, this->mData.data().begin(), shape);
        }, this->mpContainerIO);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return std::visit([this](auto pContainer) {
            std::stringstream msg;
            msg << "VariableTensorAdaptor: " << pContainer->Info() << " with " << this->mpContainer->size()
                << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s).";
            return msg.str();
        }, mpContainerIO);
        return "";
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    ContainerIOType mpContainerIO;

    ///@}
};

/// @}
} // namespace Kratos