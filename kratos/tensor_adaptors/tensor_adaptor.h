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
#include <optional>

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/nd_data.h"
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @defgroup TensorAdaptors Tensor Adaptors
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
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(TensorAdaptor);

    using ContainerPointerType = std::variant<
                                        ModelPart::DofsArrayType::Pointer,
                                        ModelPart::NodesContainerType::Pointer,
                                        ModelPart::ConditionsContainerType::Pointer,
                                        ModelPart::ElementsContainerType::Pointer,
                                        ModelPart::PropertiesContainerType::Pointer,
                                        ModelPart::MasterSlaveConstraintContainerType::Pointer,
                                        ModelPart::GeometryContainerType::Pointer
                                    >;

    using Storage = NDData<TDataType>;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Tensor Adaptor from a given @p pContainer and a @p pData
     * @details This constructor creates an instance of @ref TensorAdaptor using the @p pContainer
     *          as the container and the data in the @p pData .
     *              - Copying the internal storage of @p pData if @p Copy is true.
     *              - Sharing the internal storage of @p pData if @p Copy is false.
     * @param pContainer    Container to use in the @ref TensorAdaptor.
     * @param pData         Data to be used in the @ref TensorAdaptor.
     * @param Copy          If true, the internal storage will be copied. Otherwise, internal storage will be shared.
     */
    TensorAdaptor(
        ContainerPointerType pContainer,
        typename NDData<TDataType>::Pointer pData,
        const bool Copy = true);

    /**
     * @brief Copy constructor doubling up as an optional deep copy.
     * @details This constructor will construct an instance of a @ref TensorAdaptor by:
     *              - Copying the internal storage of @p rOther if @p Copy is true, and assigning the container pointer from @p rOther to the new instance.
     *              - Sharing the internal storage of @p rOther if @p Copy is false, and assigning the container pointer from @p rOther to the new instance.
     *
     * @param rOther    Other @ref TensorAdaptor instance.
     * @param Copy      If true, the internal storage will be copied. Otherwise, internal storage will be shared.
     */
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
     * @throws std::runtime_error if the @ref HasContainer() method returns false.
     */
    ContainerPointerType GetContainer() const;

    /**
     * @brief Returns whether this tensor adaptor has an associated valid underlying container.
     * @details If the underlying tensor adaptor does not have a valid container representation, then
     *          this method will return false. Example is @ref CombinedTensorAdaptor, where
     *          there is no valid representation of the container, so in this case this method
     *          returns false.
     */
    bool HasContainer() const;

    /**
     * @brief Return a view of the internal data structure.
     * @throws std::runtime_error If the internal data is already moved as defined in @ref Storage::ViewData.
     */
    Kratos::span<const TDataType>  ViewData() const;

    /**
     * @brief Return a view of the internal data structure.
     * @throws std::runtime_error If the internal data is already moved as defined in @ref Storage::ViewData.
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

    /**
     * @brief Returns the internal data storage.
     * @details This returns the internal storage which is a @ref NDData. This
     *          does not allow changing the internal storage, and the shape of the internal storage
     *          since @ref NDData does not have interfaces to change the shape or
     *          to change the array which it points data to.
     */
    typename Storage::Pointer pGetStorage();

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

    std::optional<ContainerPointerType> mpContainer;

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