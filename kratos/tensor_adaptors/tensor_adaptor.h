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

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the tensor adaptor types.
 */
template<class TContainerType, class TPrimitiveDataType>
class KRATOS_API(KRATOS_CORE) TensorAdaptor {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<TensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const TensorAdaptor>;

    using PrimitiveDataType = TPrimitiveDataType;

    using ContainerType = TContainerType;

    ///@}
    ///@name Life cycle
    ///@{

    TensorAdaptor(typename TContainerType::Pointer pContainer) : mpContainer(pContainer) {};

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

    typename TContainerType::Pointer GetContainer() const { return mpContainer; }

    DenseVector<TPrimitiveDataType> MoveData() { return std::move(mData); }

    Kratos::span<const TPrimitiveDataType> ViewData() const { return Kratos::span<const TPrimitiveDataType>(mData.data().begin(), mData.data().end()); }

    Kratos::span<TPrimitiveDataType> ViewData() { return Kratos::span<TPrimitiveDataType>(mData.data().begin(), mData.data().end()); }

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

    typename TContainerType::Pointer mpContainer;

    DenseVector<TPrimitiveDataType> mData;

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
template<class TContainerType, class TPrimitiveDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TensorAdaptor<TContainerType, TPrimitiveDataType>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos