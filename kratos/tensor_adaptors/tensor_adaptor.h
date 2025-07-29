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
#include "utilities/parallel_utilities.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
class TensorData
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TensorData);

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

    TensorData(
        ContainerPointerType pContainer,
        const DenseVector<unsigned int>& rShape)
        : mpContainer(pContainer),
          mShape(rShape)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mShape.empty())
            << "The tensor data shape cannot be empty. It atleast needs one dimension representing the number of items in the container [ tensor data = "
            << *this << " ].\n";

        std::visit([this](auto pContainer){
            KRATOS_ERROR_IF_NOT(mShape[0] == pContainer->size())
                << "The value of the first dimension should be equal to the number of items in the pContainer [ container size = "
                << pContainer->size() << ", tensor data = " << *this << " ].\n";
        }, mpContainer);

        // allocate new memory
        mpData = new TDataType[this->Size()];

        KRATOS_CATCH("");
    }

    /**
     * @brief Destroy the Tensor Adaptor storage
     * @details This method destroys the TensorAdaptorStorage and
     *          the internal data as well if it is still owned by the TensorAdaptorStorage.
     *
     */
    ~TensorData()
    {
        if (mpData) {
            delete[] mpData;
        }
    }

    ///@}
    ///@name Public operations
    ///@{

    TensorData::Pointer Clone() const
    {
        auto p_tensor_data = Kratos::make_intrusive<TensorData<TDataType>>(this->GetContainer(), this->Shape());
        const auto&  cloned_span = p_tensor_data->ViewData();
        const auto& origin_span = this->ViewData();

        IndexPartition<IndexType>(cloned_span.size()).for_each([&cloned_span, &origin_span](const auto Index) {
            cloned_span[Index] = origin_span[Index];
        });

        return p_tensor_data;
    }

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

    DenseVector<unsigned int> Shape() const
    {
        return mShape;
    };

    DenseVector<unsigned int> DataShape() const
    {
        const auto& shape = this->Shape();
        DenseVector<unsigned int> data_shape(shape.size() - 1);
        std::copy(shape.begin() + 1, shape.end(), data_shape.begin());
        return data_shape;
    }


    /**
     * @brief Total size of the tensor adaptor.
     */
    unsigned int Size() const
    {
        return std::accumulate(mShape.data().begin(), mShape.data().end(), 1, std::multiplies<unsigned int>{});
    }

    ContainerPointerType GetContainer() const
    {
        return mpContainer;
    }

    std::string Info() const
    {
        std::stringstream info;
        std::visit([&info, this](auto pContainer) {
            using container_type = std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;
            info << "TensorData with " << pContainer->size() << " " << ModelPart::Container<container_type>::GetEntityName() << "(s) with shape = " << this->Shape();
        }, mpContainer);
        return info.str();
    }

    ///@}


private:
    ///@name private member variables
    ///@{

    const ContainerPointerType mpContainer;

    const DenseVector<unsigned int> mShape;

    TDataType * mpData = nullptr; // nullptr is used to indicate there is no owning data.

    ///@}
    ///@name Private operations
    ///@{

    //*********************************************
    // this block is needed for refcounting in the @ref intrusive ptr
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const TensorData* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const TensorData* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }

    //*********************************************

    ///@}
};

template<class TDataType>
class TensorAdaptor {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TensorAdaptor);

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
    typename TensorData<TDataType>::ContainerPointerType GetContainer() const
    {
        return mpStorage->GetContainer();
    }

    typename TensorData<TDataType>::Pointer GetTensorData()
    {
        return mpStorage;
    }

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
        return mpStorage->MoveData();
    }

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved.
     */
    Kratos::span<const TDataType>  ViewData() const
    {
        const auto& storage = *mpStorage;
        return storage.ViewData();
    }

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved.
     */
    Kratos::span<TDataType> ViewData()
    {
        return mpStorage->ViewData();
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
        return mpStorage->Shape();
    }

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
        return mpStorage->DataShape();
    }

    /**
     * @brief Total size of the tensor adaptor.
     */
    unsigned int Size() const
    {
        return mpStorage->Size();
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
    ///@name Protected member variables
    ///@{

    typename TensorData<TDataType>::Pointer mpStorage;

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
    const TensorData<TDataType>& rThis)
{
    return rOStream << rThis.Info();
}

template<class TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TensorAdaptor<TDataType>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos