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
     * @brief Fill the internal data from Kratos data structures.
     */
    virtual void CollectData() = 0;

    /**
     * @brief Store internal data to the given Kratos data structure.
     */
    virtual void StoreData() = 0;

    /**
     * @brief Get the data container which is associated with the TensorAdaptor.
     */
    virtual ContainerType GetContainer() const = 0;

    /**
     * @brief Moves the internal data.
     * @warning The TensorAdaptor should not be used after the move is called.
     * @return DenseVector<TDataType>   Dense vector containing the internal data.
     */
    DenseVector<TDataType> MoveData() { return std::move(mData); }

    /**
     * @brief Return a view of the internal data structure.
     */
    Kratos::span<const TDataType>  ViewData() const  { return Kratos::span<const TDataType>(mData.data().begin(), mData.data().end()); }

    /**
     * @brief Return a view of the internal data structure.
     */
    Kratos::span<TDataType> ViewData() { return Kratos::span<TDataType>(mData.data().begin(), mData.data().end()); }

    /**
     * @brief Get the Shape of the tensor adaptor.
     * @details The first dimension of the tensor adaptor shape will represent the number of entities in the
     *          container. The rest of the dimensions will represent the shape of the data it caries for each of the
     *          entities.
     * @return DenseVector<unsigned int>    Shape of the tensor adaptor.
     */
    DenseVector<unsigned int> Shape() const { return mShape; };

    /**
     * @brief Get the data which the tensor carries for each of the entities.
     * @details This method returns the shape of the data which the tensor adaptor carries for the entities
     *          of the container. This will be always one dimension less than what is given in \ref Shape, because
     *          this will not contain the dimension representing number of entities in the container.
     *
     * @return DenseVector<unsigned int>    Shape of the data which is for one entity.
     */
    DenseVector<unsigned int> GetDataShape() const {
        DenseVector<unsigned int> data_shape(mShape.size() - 1);
        std::copy(mShape.begin() + 1, mShape.end(), data_shape.begin());
        return data_shape;
    }

    /**
     * @brief Total size of the tensor.
     */
    unsigned int Size() const { return std::accumulate(mShape.data().begin(), mShape.data().end(), 1, std::multiplies<unsigned int>{}); }

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