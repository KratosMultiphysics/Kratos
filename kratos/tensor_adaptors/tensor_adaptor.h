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
#include <numeric>

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
template<class TDataType>
class KRATOS_API(KRATOS_CORE) TensorAdaptor {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<TensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const TensorAdaptor>;

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

    DenseVector<TDataType> MoveData() { return std::move(mData); }

    Kratos::span<const TDataType>  ViewData() const  { return Kratos::span<const TDataType>(mData.data().begin(), mData.data().end()); }

    Kratos::span<TDataType> ViewData() { return Kratos::span<TDataType>(mData.data().begin(), mData.data().end()); }

    /**
     * @brief Get the Shape of the tensor
     */
    DenseVector<unsigned int> Shape() const { return mShape; };

    DenseVector<unsigned int> GetDataShape() const {
        DenseVector<unsigned int> data_shape(mShape.size() - 1);
        std::copy(mShape.begin() + 1, mShape.end(), data_shape.begin());
        return data_shape;
     }

    int Size() const { return std::accumulate(mShape.data().begin(), mShape.data().end(), 1, std::multiplies<int>{}); }

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string Info() const = 0;

    ///@}

protected:
    ///@name Protected member variables
    ///@{

    DenseVector<TDataType> mData;

    DenseVector<unsigned int> mShape;

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
template<class TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TensorAdaptor<TDataType>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos