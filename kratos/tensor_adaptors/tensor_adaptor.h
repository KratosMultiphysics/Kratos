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

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
/**
 * @class TensorAdaptor
 * @brief Provides an abstraction for managing tensor data associated with Kratos ModelPart containers.
 *
 * @details The TensorAdaptor class is designed to facilitate the storage, retrieval, and manipulation of tensor data
 *          linked to various Kratos ModelPart containers (such as nodes, elements, conditions, etc.). It encapsulates
 *          mechanisms for safe memory management, including reference counting via intrusive pointers, and provides
 *          interfaces for copying, moving, and viewing internal tensor data. The class also supports querying the shape
 *          and size of the tensor, as well as collecting and storing data from/to Kratos data structures.
 *
 * @tparam TDataType The type of the data stored in the tensor adaptor.
 */
class KRATOS_API(KRATOS_CORE) TensorAdaptor {
protected:
    ///@name Class definitions
    ///@{

    /**
     * @class Storage
     * @brief Manages the storage and lifetime of tensor adaptor data associated with various Kratos ModelPart containers.
     *
     * The Storage class encapsulates a pointer to a container (such as nodes, elements, conditions, etc.) and manages
     * the associated tensor data. It provides mechanisms for copying, moving, and viewing the internal data, as well as
     * querying the shape and size of the tensor. Reference counting is implemented for safe memory management using
     * intrusive pointers.
     *
     * @tparam TDataType The type of the data stored in the tensor adaptor.
     */
    class Storage
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Storage);

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

        Storage(
            ContainerPointerType pContainer,
            const DenseVector<unsigned int>& rShape);

        /**
         * @brief Destroy the Tensor Adaptor storage
         * @details This method destroys the TensorAdaptorStorage and
         *          the internal data as well if it is still owned by the TensorAdaptorStorage.
         *
         */
        ~Storage();

        ///@}
        ///@name Public operations
        ///@{

        /**
         * @brief Creates a copy of the current storage object.
         * @details This method returns a pointer to a new Storage object that is having copied data and copied
         *          pointer to the container.
         * @return Storage::Pointer A pointer to the copied Storage object.
         */
        Storage::Pointer Copy() const;

        /**
         * @brief Moves the internal data.
         * @warning The TensorAdaptor should not be used after the move is called.
         *          The management of the data should be carried out by the owner of the
         *          returned Kratos::span. Otherwise, there will be memory leaks.
         * @throws If the internal data is already moved.
         * @return Kratos::span<TDataType>  Returns a span containing the internal data.
         */
        Kratos::span<TDataType> MoveData();

        /**
         * @brief Return a view of the internal data structure.
         * @throws If the internal data is already moved.
         */
        Kratos::span<const TDataType>  ViewData() const;

        /**
         * @brief Return a view of the internal data structure.
         * @throws If the internal data is already moved.
         */
        Kratos::span<TDataType> ViewData();

        /**
         * @brief Returns the shape of the tensor as a dense vector of unsigned integers.
         * @details This function provides the dimensions of the tensor, where each element in the returned
         *          DenseVector corresponds to the size of the tensor in that particular dimension. First
         *          dimension always represents how many entities are in the container.
         *
         * @return DenseVector<unsigned int> A vector containing the size of each dimension of the tensor.
         */
        DenseVector<unsigned int> Shape() const;

        /**
         * @brief Returns the shape of the underlying tensor data as a dense vector.
         * @details This method returns the shape of the data which may be collected and stored
         *          in each entity of the specified container. This is always one dimension less
         *          than the @ref Storage::Shape method.
         * @return DenseVector<unsigned int> A vector containing the dimensions of the tensor.
         */
        DenseVector<unsigned int> DataShape() const;


        /**
         * @brief Returns the number of elements in the tensor adaptor.
         * @return The size as an unsigned integer.
         */
        unsigned int Size() const;

        /**
         * @brief Returns a pointer to the underlying container.
         * @details This method provides access to the internal container used by the tensor adaptor.
         *          The returned pointer allows read-only operations on the container.
         *
         * @return ContainerPointerType Pointer to the underlying container.
         */
        ContainerPointerType GetContainer() const;


        /**
         * @brief Returns a string containing information about the tensor adaptor.
         * @return A std::string with descriptive information about the current state or properties of the tensor adaptor.
         */
        std::string Info() const;
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

        friend void intrusive_ptr_add_ref(const Storage* x)
        {
            x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
        }

        friend void intrusive_ptr_release(const Storage* x)
        {
            if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
                std::atomic_thread_fence(std::memory_order_acquire);
                delete x;
            }
        }

        //*********************************************

        ///@}
    };

    ///@}
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TensorAdaptor);

    ///@}
    ///@name Life cycle
    ///@{

    TensorAdaptor(
        const TensorAdaptor& rOther,
        const bool Copy = false);

    virtual ~TensorAdaptor() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Check if necessary data is present in the underlying pointer vector sets.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only check data from Kratos data structures whether necessary information is there to perform
     *          CollectData and StoreData without any errors.
     */
    virtual void Check() const;

    /**
     * @brief Fill the internal data from Kratos data structures.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only gather data from Kratos data structures to the TensorAdaptor.
     */
    virtual void CollectData();

    /**
     * @brief Store internal data to the given Kratos data structure.
     */
    virtual void StoreData();

    /**
     * @brief Get the data container which is associated with the TensorAdaptor.
     */
    typename Storage::ContainerPointerType GetContainer() const;

    /**
     * @brief Moves the internal data.
     * @warning The TensorAdaptor should not be used after the move is called.
     *          The management of the data should be carried out by the owner of the
     *          returned Kratos::span. Otherwise, there will be memory leaks.
     * @throws If the internal data is already moved.
     * @return Kratos::span<TDataType>  Returns a span containing the internal data.
     */
    Kratos::span<TDataType> MoveData();

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved.
     */
    Kratos::span<const TDataType>  ViewData() const;

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved.
     */
    Kratos::span<TDataType> ViewData();

    /**
     * @brief Get the Shape of the tensor adaptor.
     * @details The first dimension of the tensor adaptor shape will represent the number of entities in the
     *          container. The rest of the dimensions will represent the shape of the data it caries for each of the
     *          entities.
     * @return DenseVector<unsigned int>    Shape of the tensor adaptor.
     */
    DenseVector<unsigned int> Shape() const;

    /**
     * @brief Get the shape of the data which the tensor carries for each of the entities.
     * @details This method returns the shape of the data which the tensor adaptor carries for the entities
     *          of the container. This will be always one dimension less than what is given in \ref Shape, because
     *          this will not contain the dimension representing number of entities in the container.
     *
     * @return DenseVector<unsigned int>    Shape of the data which is for one entity.
     */
    DenseVector<unsigned int> DataShape() const;

    /**
     * @brief Total size of the tensor adaptor.
     */
    unsigned int Size() const;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Provides some information about the TensorAdaptor.
     */
    virtual std::string Info() const;

    ///@}

protected:
    ///@name Protected life cycle constructors
    ///@{

    TensorAdaptor() = default;

    ///@}
    ///@name Protected member variables
    ///@{

    typename Storage::Pointer mpStorage;

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