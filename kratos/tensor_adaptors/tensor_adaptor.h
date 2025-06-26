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
#include <atomic>
#include <variant>
#include <type_traits>

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the tensor adaptor types.
 */
class KRATOS_API(KRATOS_CORE) TensorAdaptor {
private:
    ///@name Private classes
    ///@{

    template<class... TPrimitiveDataTypes>
    struct PrimitiveDataTypes
    {
        ///@name Type definitions
        ///@{

        using DataType = std::variant<DenseVector<TPrimitiveDataTypes>...>;

        using ViewType = std::variant<Kratos::span<TPrimitiveDataTypes>...>;

        using ConstViewType = std::variant<Kratos::span<const TPrimitiveDataTypes>...>;

        ///@}
    };

    ///@}
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<TensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const TensorAdaptor>;

    using PrimitiveDataTypeInfo = PrimitiveDataTypes<bool, int, double>;

    using DataType = PrimitiveDataTypeInfo::DataType;

    using ViewType = PrimitiveDataTypeInfo::ViewType;

    using ConstViewType = PrimitiveDataTypeInfo::ConstViewType;

    using ContainerType = std::variant<
                                ModelPart::NodesContainerType::Pointer,
                                ModelPart::ConditionsContainerType::Pointer,
                                ModelPart::ElementsContainerType::Pointer,
                                ModelPart::PropertiesContainerType::Pointer,
                                // ModelPart::MasterSlaveConstraintContainerType::Pointer,
                                ModelPart::GeometryContainerType::GeometriesMapType::Pointer
                            >;

    ///@}
    ///@name Life cycle
    ///@{

    TensorAdaptor() = default;

    virtual ~TensorAdaptor() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data from Kratos data structures
     */
    virtual void CollectData() = 0;

    /**
     * @brief Store internal data to the given TContainerType container.
     *
     */
    virtual void StoreData() = 0;

    virtual ContainerType GetContainer() const = 0;

    DataType MoveData() { return std::move(mData); }

    ConstViewType  ViewData() const
    {
        return std::visit([](const auto& rData) -> ConstViewType {
            using data_type = typename std::remove_cv_t<std::decay_t<decltype(rData)>>::value_type;
            return Kratos::span<const data_type>(rData.data().begin(), rData.data().end());
        }, this->mData);
    }

    ViewType ViewData()
    {
        return std::visit([](auto& rData) -> ViewType {
            using data_type = typename std::remove_cv_t<std::decay_t<decltype(rData)>>::value_type;
            return Kratos::span<data_type>(rData.data().begin(), rData.data().end());
        }, this->mData);
    }

    /**
     * @brief Get the Shape of the tensor
     */
    DenseVector<int> Shape() const { return mShape; };

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string Info() const = 0;

    ///@}

protected:
    ///@name Protected member variables
    ///@{

    DataType mData;

    DenseVector<int> mShape;

    ///@}

private:
    ///@name Private operations
    ///@{

    //*********************************************
    // this block is needed for refcounting in the @ref intrusive ptr
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const TensorAdaptor* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const TensorAdaptor* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }

    //*********************************************

    ///@}
};

/// @}
/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TensorAdaptor& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos