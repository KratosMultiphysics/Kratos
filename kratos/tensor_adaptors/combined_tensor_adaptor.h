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
 *          Example use case 1: (Extending tensor adaptors in axis = 1)
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
 *          Example use case 2: (Extending tensor adaptors in axis = 0)
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
 *                  combined_tensor_adaptor = Kratos.TensorAdaptor.DoubleCombinedTensorAdaptor([tensor_adaptor_1, tensor_adaptor_2])
 *              @endcode
 *
 *          Example use case 3: (Extending tensor adaptors by raveling)
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
 * @section CombinedTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::DofsArrayType
 * - @ref ModelPart::NodesContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 * - @ref ModelPart::PropertiesContainerType
 * - @ref ModelPart::GeometryContainerType
 * - @ref ModelPart::MasterSlaveConstraintContainerType
 *
 * @author Suneth Warnakulasuriya
 * @tparam TDataType The type of the data stored in the tensor adaptor.
 */
template<class TDataType>
class KRATOS_API(KRATOS_CORE) CombinedTensorAdaptor : public TensorAdaptor<TDataType> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CombinedTensorAdaptor);

    using BaseType = TensorAdaptor<TDataType>;

    using TensorAdaptorVectorType = std::vector<typename BaseType::Pointer>;

    ///@}
    ///@name Life cycle
    ///@{

    CombinedTensorAdaptor(
        const TensorAdaptorVectorType& rTensorAdaptorVector,
        const bool CollectAndStoreRecursively = true);

    CombinedTensorAdaptor(
        const TensorAdaptorVectorType& rTensorAdaptorVector,
        const unsigned int Axis,
        const bool CollectAndStoreRecursively = true);

    CombinedTensorAdaptor(
        const TensorAdaptorVectorType& rTensorAdaptorVector,
        const bool Ravel,
        const bool CollectAndStoreRecursively = true);

    CombinedTensorAdaptor(
        const CombinedTensorAdaptor& rOther,
        const bool CollectAndStoreRecursively = true,
        const bool Copy = true);

    virtual ~CombinedTensorAdaptor() = default;

    ///@}
    ///@name Public operations
    ///@{

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

    const bool mCollectAndStoreRecursively;

    const int mAxis;

    TensorAdaptorVectorType mTensorAdaptors;

    ///@}
};

///@}

} // namespace Kratos