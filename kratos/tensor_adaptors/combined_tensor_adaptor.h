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
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

 /**
 * @class CombinedTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Provides an abstraction for managing tensor data associated with Kratos ModelPart containers.
 * @details This class provides mechanisms to combine multiple @ref TensorAdaptor instances to one
 *          @ref CombinedTensorAdaptor instance. So, this will enable working on multiple tensor adaptors
 *          as a single @ref TensorAdaptor, and interfacing with numpy or C arrays.
 *
 * @note    @ref CombinedTensorAdaptor always keeps a copy of all the internal data of the sub @ref TensorAdaptor instances
 *          which are passed as @p rTensorAdaptorVector.
 *
 *          Example use case 1: (Extending tensor adaptors by raveling)
 *              In this example, PRESSURE of nodes in one model part with VELOCITY of nodes of another model part are combined together.
 *              @code{.py}
 *                  import KratosMultiphysics as Kratos
 *                  model = Kratos.Model()
 *                  model_part_1 = model.CreateModelPart("test1")
 *                  model_part_2 = model.CreateModelPart("test2")
 *                  for i in range(10):
 *                      node_1 = model_part_1.CreateNewNode(i + 1, 0, 0, 0)
 *                      node_1.SetValue(Kratos.PRESSURE, i)
 *                      if (i % 2 == 0):
 *                          node_2 = model_part_2.CreateNewNode(i + 1, 0, 0, 0)
 *                          node_2.SetValue(Kratos.VELOCITY, Kratos.Array3D([i, i, i]))
 *
 *                  # following will be tensor adaptor with shape = [10, 1]
 *                  tensor_adaptor_1 = Kratos.TensorAdaptor.VariableTensorAdaptor(self.model_part_1.Nodes, Kratos.PRESSURE)
 *
 *                  # following will be a tensor adaptor with shape = [5, 3]
 *                  tensor_adaptor_2 = Kratos.TensorAdaptor.VariableTensorAdaptor(self.model_part_2.Nodes, Kratos.VELOCITY)
 *
 *                  # following will be a tensor adaptor with shape = [25, 1]
 *                  combined_tensor_adaptor = Kratos.TensorAdaptor.DoubleCombinedTensorAdaptor([tensor_adaptor_1, tensor_adaptor_2])
 *              @endcode
 *
 *          Example use case 2: (Extending tensor adaptors in axis = 1)
 *              In this example, PRESSURE, VELOCITY of nodes are combined together.
 *              @code{.py}
 *                  import KratosMultiphysics as Kratos
 *                  model = Kratos.Model()
 *                  model_part = model.CreateModelPart("test")
 *                  for i in range(10):
 *                      node = model_part.CreateNewNode(i + 1, 0, 0, 0)
 *                      node.SetValue(Kratos.PRESSURE, i)
 *                      node.SetValue(Kratos.VELOCITY, Kratos.Array3D([i, i, i]))
 *
 *                  # following will be tensor adaptor with shape = [10]
 *                  pressure_tensor_adaptor = Kratos.TensorAdaptor.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
 *
 *                  # following will be a tensor adaptor with shape = [10, 3]
 *                  velocity_tensor_adaptor = Kratos.TensorAdaptor.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
 *
 *                  # following will be a tensor adaptor with shape = [10, 4]
 *                  combined_tensor_adaptor = Kratos.TensorAdaptor.DoubleCombinedTensorAdaptor([pressure_tensor_adaptor, velocity_tensor_adaptor], axis=1)
 *              @endcode
 *
 *          Example use case 3: (Extending tensor adaptors in axis = 0)
 *              In this example, VELOCITY of nodes from two different model parts are combined together
 *              @code{.py}
 *                  import KratosMultiphysics as Kratos
 *                  model = Kratos.Model()
 *                  model_part_1 = model.CreateModelPart("test1")
 *                  model_part_2 = model.CreateModelPart("test2")
 *                  for i in range(10):
 *                      node_1 = model_part_1.CreateNewNode(i + 1, 0, 0, 0)
 *                      node_1.SetValue(Kratos.VELOCITY, Kratos.Array3D([i, i, i]))
 *                      if (i % 2 == 0):
 *                          node_2 = model_part_2.CreateNewNode(i + 1, 0, 0, 0)
 *                          node_2.SetValue(Kratos.VELOCITY, Kratos.Array3D([i, i, i]))
 *
 *                  # following will be tensor adaptor with shape = [10, 3]
 *                  tensor_adaptor_1 = Kratos.TensorAdaptor.VariableTensorAdaptor(self.model_part_1.Nodes, Kratos.VELOCITY)
 *
 *                  # following will be a tensor adaptor with shape = [5, 3]
 *                  tensor_adaptor_2 = Kratos.TensorAdaptor.VariableTensorAdaptor(self.model_part_2.Nodes, Kratos.VELOCITY)
 *
 *                  # following will be a tensor adaptor with shape = [15, 3]
 *                  combined_tensor_adaptor = Kratos.TensorAdaptor.DoubleCombinedTensorAdaptor([tensor_adaptor_1, tensor_adaptor_2], axis=0)
 *              @endcode
 *
 * @section CombinedTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::DofsArrayType
 * - @ref ModelPart::NodesContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 * - @ref ModelPart::PropertiesContainerType
 * - @ref ModelPart::GeometryContainerType
 * - @ref ModelPart::MasterSlaveConstraintContainerType
 *
 * @section CombinedTensorAdaptor_Usage Usage
 * - @ref CollectData() - To read in the data from sub tensor adaptors to the @ref CombinedTensorAdaptor. If @p PerformCollectDataRecursively is true, then this will first call @ref CollectData() of the sub tensor adaptors.
 * - @ref StoreData() - To store back the data from @ref CombinedTensorAdaptor to sub tensor adaptors. If @p PerformStoreDataRecursively is true, then this will call sub tensor adaptors' @ref StoreData() afterwards.
 * - @ref GetContainer() - Not used. throws an error
 * - @ref HasContainer() - Can be used to check if the container is there. Returns false for CombinedTensorAdaptor.
 *
 * @author Suneth Warnakulasuriya
 */
template<class TDataType>
class KRATOS_API(KRATOS_CORE) CombinedTensorAdaptor : public TensorAdaptor<TDataType> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CombinedTensorAdaptor);

    using BaseType = TensorAdaptor<TDataType>;

    using TensorAdaptorVectorType = std::vector<typename BaseType::Pointer>;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Combined Tensor Adaptor given list of @ref TensorAdaptor instances by raveling them.
     * @details This will construct a @ref CombinedTensorAdaptor by reveling all the tensor adaptors given in @p rTensorAdaptorVector.
     *          Final shape of the CombinedTensorAdaptor will be summation of @ref TensorAdaptor::Size() values from each tensor adaptor
     *          in @p rTensorAdaptorVector.
     *
     * @note    This always return a flat tensor adaptor having only one dimensionality.
     *
     * @param rTensorAdaptorVector                          List of tensor adaptors.
     * @param PerformCollectDataRecursively    If true, @ref CollectData() method will call combined tensor adaptors' @ref CollectData.
     * @param PerformStoreDataRecursively      If true, @ref StoreData() method will call combined tensor adaptors' @ref StoreData.
     */
    CombinedTensorAdaptor(
        const TensorAdaptorVectorType& rTensorAdaptorVector,
        const bool PerformCollectDataRecursively = true,
        const bool PerformStoreDataRecursively = true,
        const bool Copy = true);

    /**
     * @brief Construct a new Combined Tensor Adaptor  given list of @ref TensorAdaptor instances along the specified @p Axis.
     * @details This constructor will construct a @ref TensorAdaptor, by concatenating the given
     *          list of tensor adaptors on the specified @p Axis.
     *
     * @note   If there are tensor adaptors with different number of dimensionalities, then the final number of dimensions of the resulting
     *         CombineTensorAdaptor will be having the same number of dimensions of a tensor adaptor which has the maximum number of dimensions.
     *         All the tensor adaptors which has less number of dimensions will have a modifed shape appending @p 1 to its original shape making all
     *         the tensor adaptors having the same number of dimensions at the end.
     *
     *          Example:
     *          @code{.py}
     *            a = VariableTensorAdaptor(nodes, PRESSURE) -> shape = [10]
     *            b = VariableTensorAdaptor(nodes, VELOCITY) -> shape = [10, 3]
     *            # When combining, a's shape will be modified to [10, 1].
     *          @endcode
     *
     * @throws std::runtime_error if the @p rTensorAdaptorVector is empty.
     * @throws std::runtime_error if the @p Axis is less than zero.
     * @throws std::runtime_error if the @p Axis is equal or larger than the largest number of dimensions found in the tensor adaptors of @p rTensorAdaptorVector.
     * @throws std::runtime_error if the number of components in each dimension except for axis dimension is not equal on all the tensor adaptors in @p rTensorAdaptorVector.
     *
     * @param rTensorAdaptorVector                          List of tensor adaptors.
     * @param Axis                                          Axis to concatenate the tensor adaptors.
     * @param PerformCollectDataRecursively    If true, @ref CollectData() method will call combined tensor adaptors' @ref CollectData.
     * @param PerformStoreDataRecursively      If true, @ref StoreData() method will call combined tensor adaptors' @ref StoreData.
     */
    CombinedTensorAdaptor(
        const TensorAdaptorVectorType& rTensorAdaptorVector,
        const unsigned int Axis,
        const bool PerformCollectDataRecursively = true,
        const bool PerformStoreDataRecursively = true,
        const bool Copy = true);

    /**
     * @brief Construct a new Combined Tensor Adaptor based on an existing @p rOther.
     *
     * @param rOther                                        Other @ref CombinedTensorAdaptor to copy/refer from.
     * @param PerformCollectDataRecursively    If true, @ref CollectData() method will call combined tensor adaptors' @ref CollectData.
     * @param PerformStoreDataRecursively      If true, @ref StoreData() method will call combined tensor adaptors' @ref StoreData.
     * @param Copy                                          If true, the underlying internal data will be copied along with the pointers to the containers.
     */
    CombinedTensorAdaptor(
        const CombinedTensorAdaptor& rOther,
        const bool PerformCollectDataRecursively = true,
        const bool PerformStoreDataRecursively = true,
        const bool Copy = true);

    virtual ~CombinedTensorAdaptor() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones the existing tensor adaptor.
     */
    TensorAdaptor<TDataType>::Pointer Clone() const override;

    /**
     * @brief Check if the necessary data is present in the underlying container.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only make sure that @ref CombinedTensorAdaptor::CollectData "CollectData" and @ref CombinedTensorAdaptor::StoreData "StoreData"
     *          can be performed without errors.
     */
    void Check() const override;

    /**
     * @brief Fill the internal data from Kratos data structures.
     * @details This method should not change anything in the underlying Kratos data structures. It should
     *          only gather data from Kratos data structures to the @ref CombinedTensorAdaptor instance.
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given Kratos data structure.
     */
    void StoreData() override;

    /**
     * @brief Get the list of tensor adaptors stored in this CombinedTensorAdaptor
     */
    TensorAdaptorVectorType GetTensorAdaptors() const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const bool mPerformCollectDataRecursively;

    const bool mPerformStoreDataRecursively;

    const int mAxis;

    TensorAdaptorVectorType mTensorAdaptors;

    ///@}
};

///@}

} // namespace Kratos