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

/**
 * @defgroup TensorAdaptors "Tensor Adaptors"
 * @brief Tensor adaptors are used to read or write data from/to entities (nodes, elements, conditions, constraints) of a @ref ModelPart.
 * @details The @ref TensorAdaptor class is designed to handle reading from variables, info, data from each entity
 *          in the containers of @p ModelPart to a @p TensorAdaptor instance or vice-versa. It encapsulates
 *          interfaces for copying, moving, and viewing internal data.
 *
 *          Example in C++:
 *          @code
 *              auto model = Model();
 *              auto& model_part = model.CreateModelPart("test");
 *              for (IndexType i = 0; i < 10; ++i) {
 *                  model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0)->SetValue(PRESSURE, i + 1);
 *              }
 *
 *              VariableTensorAdaptor var_ta(model_part.pNodes(), PRESSURE);
 *              var_ta.Check();       // Checks if the PRESSURE variable is available in each of the nodes.
 *              var_ta.CollectData(); // Reads the PRESSURE variable from each of the nodes and puts them into an internal storage.
 *
 *              auto data = var_ta.ViewData(); // Create a span viewing the data
 *
 *              // here we modify the internal data storage of the TensorAdaptor.
 *              // This will not modify the PRESSURE of the nodes.
 *              for (IndexType i = 0; i < data.size(); ++i) {
 *                  data[i] += 1.0;
 *              }
 *
 *              var_ta.StoreData(); // Writes back the modified tensor adaptor data back to each node's PRESSURE variable.
 *          @endcode
 *
 *          Equivalent python code:
 *          @code{.py}
 *              import KratosMultiphysics as Kratos
 *
 *              model = Kratos.Model()
 *              model_part = model.CreateModelPart("test")
 *              for i in range(10):
 *                  model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0).SetValue(Kratos.PRESSURE, i + 1)
 *
 *              var_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes(), Kratos.PRESSURE)
 *              var_ta.Check()
 *              var_ta.CollectData()
 *              var_ta.data += 1.0
 *              var_ta.StoreData()
 *          @endcode
 */

 /**
 * @class TensorAdaptor
 * @ingroup TensorAdaptors
 *
 * @section TensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::DofsArrayType
 * - @ref ModelPart::NodesContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 * - @ref ModelPart::PropertiesContainerType
 * - @ref ModelPart::GeometryContainerType
 * - @ref ModelPart::MasterSlaveConstraintContainerType
 *
 * @author Suneth Warnakulasuriya
 * @brief Provides an abstraction for managing tensor data associated with Kratos ModelPart containers.
 * @tparam TDataType The type of the data stored in the tensor adaptor.
 */
template<class TDataType>
class KRATOS_API(KRATOS_CORE) TensorAdaptor {
protected:

    /**
     * @class Storage
     * @brief Manages the storage and lifetime of tensor adaptor data associated with various containers in @ref ModelPart.
     *
     * @p Storage encapsulates a pointer to a container (such as nodes, elements, conditions, etc.) and manages
     * the associated tensor data. It provides mechanisms for copying, moving, and viewing the internal data, as well as
     * querying the shape and size of the tensor.
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
                                            ModelPart::DofsArrayType::Pointer,
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
         *          the internal data as well if it is still owned by @p Storage.
         *
         */
        ~Storage();

        ///@}
        ///@name Public operations
        ///@{

        /**
         * @brief Creates a copy of the current storage object.
         * @details This method returns a pointer to a new @p Storage which is having the copied internal data and the copied
         *          pointer to the container.
         */
        Storage::Pointer Copy() const;

        /**
         * @brief Moves the internal data.
         * @warning The TensorAdaptor should not be used after the move is called.
         *          The management of the data should be carried out by the owner of the
         *          returned Kratos::span. Otherwise, there will be memory leaks.
         * @throws If the internal data is already moved.
         * @return Returns a span containing the internal data.
         */
        Kratos::span<TDataType> MoveData();

        /**
         * @brief Return a view of the internal data structure.
         * @throws If the internal data was moved via @ref Storage::MoveData.
         */
        Kratos::span<const TDataType>  ViewData() const;

        /**
         * @brief Return a view of the internal data structure.
         * @throws If the internal data is already moved.
         */
        Kratos::span<TDataType> ViewData();

        /**
         * @brief Returns the shape of the tensor.
         * @details This function provides the dimensions of the tensor, where each element in the returned
         *          array corresponds to the size of the tensor in that particular dimension. The first
         *          dimension always represents the number of entities stored in the container.
         *
         * @return A vector containing the size of each dimension of the tensor.
         */
        DenseVector<unsigned int> Shape() const;

        /**
         * @brief Returns the shape of the underlying tensor data.
         * @details This method returns the shape of the data which may be collected and stored
         *          in each entity of the specified container. This is always one dimension less
         *          than the @ref Storage::Shape method.
         * @return A vector containing the dimensions of the tensor.
         */
        DenseVector<unsigned int> DataShape() const;


        /**
         * @brief Returns the number of elements in the tensor adaptor.
         */
        unsigned int Size() const;

        /**
         * @brief Returns a pointer to the underlying container.
         * @details This method provides access to the internal container used by the tensor adaptor.
         *          The returned pointer allows read-only operations on the container.
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

public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TensorAdaptor);

    ///@}
    ///@name Life cycle
    ///@{

    TensorAdaptor(
        const TensorAdaptor& rOther,
        const bool Copy = true);

    virtual ~TensorAdaptor() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Check if the necessary data is present in the underlying container.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only make sure that @ref TensorAdaptor::CollectData "CollectData" and @ref TensorAdaptor::StoreData "StoreData"
     *          can be performed without errors.
     */
    virtual void Check() const;

    /**
     * @brief Fill the internal data from Kratos data structures.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only gather data from Kratos data structures to the @ref TensorAdaptor instance.
     */
    virtual void CollectData();

    /**
     * @brief Store internal data to the given Kratos data structure.
     */
    virtual void StoreData();

    /**
     * @brief Get the entity container which is associated with the @p TensorAdaptor instance.
     */
    typename Storage::ContainerPointerType GetContainer() const;

    /**
     * @brief Moves the internal data.
     * @warning The TensorAdaptor should not be used after the move is called.
     *          The management of the data should be carried out by the owner of the
     *          returned @ref Kratos::span "span". Otherwise, there will be memory leaks.
     * @throws If the internal data is already moved as defined in @ref Storage::MoveData.
     * @return Kratos::span<TDataType>  Returns a span containing the internal data.
     */
    Kratos::span<TDataType> MoveData();

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved as defined in @ref Storage::ViewData.
     */
    Kratos::span<const TDataType>  ViewData() const;

    /**
     * @brief Return a view of the internal data structure.
     * @throws If the internal data is already moved as defined in @ref Storage::ViewData.
     */
    Kratos::span<TDataType> ViewData();

    /**
     * @brief Get the Shape of the tensor adaptor.
     * @details The first dimension of the tensor adaptor shape will represent the number of entities in the
     *          container. The rest of the dimensions will represent the shape of the data it caries for each of the
     *          entities.
     */
    DenseVector<unsigned int> Shape() const;

    /**
     * @brief Get the shape of the data which the tensor carries for each of the entities.
     * @details This method returns the shape of the data which the tensor adaptor carries for the entities
     *          of the container. This will be always one dimension less than what is given in \ref Shape, because
     *          this will not contain the dimension representing number of entities in the container.
     *
     * @return DenseVector  Shape of the data which is for one entity.
     */
    DenseVector<unsigned int> DataShape() const;

    /**
     * @brief Total size of the tensor adaptor.
     */
    unsigned int Size() const;

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string Info() const;

    ///@}

protected:
    ///@name Protected life cycle
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

///@}
template<class TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TensorAdaptor<TDataType>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos