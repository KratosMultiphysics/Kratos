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
#include "includes/model_part.h"
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the tensor adaptor types.
 */
template<class TDataType>
class TensorAdaptor {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TensorAdaptor);

    using ContainerPointerType = std::variant<
                                        ModelPart::NodesContainerType::Pointer,
                                        ModelPart::ConditionsContainerType::Pointer,
                                        ModelPart::ElementsContainerType::Pointer,
                                        ModelPart::PropertiesContainerType::Pointer,
                                        ModelPart::MasterSlaveConstraintContainerType::Pointer,
                                        ModelPart::GeometryContainerType::Pointer
                                    >;

    ///@}
    ///@name Life cycle
    ///@{

    TensorAdaptor() = default;

    /**
     * @brief Destroy the Tensor Adaptor
     * @details This method destroys the TensorAdaptor and
     *          the internal data as well if it is still owned by the TensorAdaptor.
     *
     */
    virtual ~TensorAdaptor()
    {
        if (mpData) {
            delete[] mpData;
        }
    }

    ///@}
    ///@name Public operations
    ///@{

    virtual TensorAdaptor::Pointer Clone() const = 0;

    /**
     * @brief Fill the internal data from Kratos data structures.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only gather data from Kratos data structures to the TensorAdaptor.
     */
    virtual void CollectData() = 0;

    /**
     * @brief Store internal data to the given Kratos data structure.
     */
    virtual void StoreData() = 0;

    /**
     * @brief Get the data container which is associated with the TensorAdaptor.
     */
    virtual ContainerPointerType GetContainer() const = 0;

    /**
     * @brief Moves the internal data.
     * @warning The TensorAdaptor should not be used after the move is called.
     *          The management of the data should be carried out by the owner of the
     *          returned Kratos::span. Otherwise, there will be memory leaks.
     * @throws If the internal data is already moved.
     * @return Kratos::span<TDataType>  Returns a span containing the internal data.
     */
    Kratos::span<TDataType> MoveData()
    {
        KRATOS_ERROR_IF_NOT(mpData) << "The data is already moved [ " << *this << " ].\n";
        auto p_data = mpData;
        mpData = nullptr;
        return Kratos::span<TDataType>(p_data, p_data + this->Size());
    }

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved.
     */
    Kratos::span<const TDataType>  ViewData() const
    {
        KRATOS_ERROR_IF_NOT(mpData) << "The data is already moved [ " << *this << " ].\n";
        return Kratos::span<const TDataType>(mpData, mpData + this->Size());
    }

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved.
     */
    Kratos::span<TDataType> ViewData()
    {
        KRATOS_ERROR_IF_NOT(mpData) << "The data is already moved [ " << *this << " ].\n";
        return Kratos::span<TDataType>(mpData, mpData + this->Size());
    }

    /**
     * @brief Get the Shape of the tensor adaptor.
     * @details The first dimension of the tensor adaptor shape will represent the number of entities in the
     *          container. The rest of the dimensions will represent the shape of the data it caries for each of the
     *          entities.
     * @return DenseVector<unsigned int>    Shape of the tensor adaptor.
     */
    DenseVector<unsigned int> Shape() const
    {
        return mShape;
    };

    /**
     * @brief Get the shape of the data which the tensor carries for each of the entities.
     * @details This method returns the shape of the data which the tensor adaptor carries for the entities
     *          of the container. This will be always one dimension less than what is given in \ref Shape, because
     *          this will not contain the dimension representing number of entities in the container.
     *
     * @return DenseVector<unsigned int>    Shape of the data which is for one entity.
     */
    DenseVector<unsigned int> DataShape() const
    {
        DenseVector<unsigned int> data_shape(mShape.size() - 1);
        std::copy(mShape.begin() + 1, mShape.end(), data_shape.begin());
        return data_shape;
    }

    /**
     * @brief Total size of the tensor adaptor.
     */
    unsigned int Size() const
    {
        return std::accumulate(mShape.data().begin(), mShape.data().end(), 1, std::multiplies<unsigned int>{});
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Provides some information about the TensorAdaptor.
     */
    virtual std::string Info() const = 0;

    ///@}

protected:
    ///@name Protected operations
    ///@{

    void SetShape(const DenseVector<unsigned int>& rNewShape)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rNewShape.empty())
            << "The tensors cannot have an empty shape. It should have atleast "
            << "first dimension which represents number of items in the container [ "
            << *this << " ].\n";

        // checks whether the mShape is empty, that means
        // the SetShape method is not yet called. This is because, once this method is called
        // the mShape size cannot be zero. It will at least have size 1 representing the size of the container.
        KRATOS_ERROR_IF_NOT(this->Shape().empty())
            << "The tensor is already initialized with a shape [ new shape = " << rNewShape
            << ", current shape = " << this->Shape() << ", " << *this << " ].\n";

        // set the shape
        mShape = rNewShape;

        // get the new size
        const auto new_size = this->Size();

        // allocate new memory
        mpData = new TDataType[new_size];

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    TDataType * mpData = nullptr; // nullptr is used to indicate there is no owning data.

    DenseVector<unsigned int> mShape;

    ///@}
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