//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "tensor_adaptor.h"
#include "includes/ublas_interface.h"
#include "utilities/data_type_traits.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class FixityTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling fixity of nodal DOFs
 *
 * @details This class provides an interface to collect and store fixity from
 * the nodal DOFs of a ModelPart. It extends TensorAdaptor<bool>.
 * The fixity is stored in a boolean tensor where each entry is true if the
 * corresponding DOF is fixed and false otherwise.
 * The data size is given by the number of nodes in the container times the
 * number of DOF variables in the fixity variables list.
 * Note that thte fixity data is sorted according to the order of the DOF
 * variables in the fixity variables list.
 *
 * @section FixityTensorAdaptor_supported_container Supported container
 * types
 * - @ref ModelPart::NodesContainerType
 *
 * @section FixityTensorAdaptor_usage Usage
 * - Use @ref Check to verify the container is valid.
 * - Use @ref CollectData to read fixity from the nodal entities.
 * - Use @ref StoreData to set fixity to the nodal entities.
 *
 * @author Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) FixityTensorAdaptor : public TensorAdaptor<bool>
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FixityTensorAdaptor);

    using BaseType = TensorAdaptor<bool>;

    ///@}
    ///@name Life cycle
    ///@{

    FixityTensorAdaptor(
        ContainerPointerType pContainer,
        const std::vector<const Variable<double>*>& rDofsVarPointerList);

    FixityTensorAdaptor(
        const BaseType& rOther,
        const bool Copy = true);

    FixityTensorAdaptor(const FixityTensorAdaptor& rOther) = default;

    // Destructor
    ~FixityTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    TensorAdaptor::Pointer Clone() const override;

    void Check() const override;

    void CollectData() override;

    void StoreData() override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}
    private:

    ///@name Private members
    ///@{

    const std::vector<const Variable<double>*> mDofsVarPointerList;

    ///@}
    ///@name Private static operations
    ///@{

    /**
     * @brief Get the shape of the tensor adaptor data
     * The shape is determined by the number of nodes and the number of DOFs in the fixity variables list
     * @param rContainer Nodal container
     * @return DenseVector<unsigned int> Shape of the tensor adaptor data
     */
    DenseVector<unsigned int> GetShape(const ModelPart::NodesContainerType& rContainer);

    ///@}
    ///@name Private static operations
    ///@{

    /**
     * @brief Collect fixity data from the nodal entities
     * @param Span Span of the tensor adaptor data
     * @param rDofsVarPointerList List of DOF variables pointers
     * @param rNodes Nodal container
     */
    static void CollectFixityData(
        Kratos::span<bool> Span,
        const std::vector<const Variable<double>*>& rDofsVarPointerList,
        const ModelPart::NodesContainerType& rNodes);

    /**
     * @brief Store fixity data to the nodal entities
     * @param Span Span of the tensor adaptor data
     * @param rDofsVarPointerList List of DOF variables pointers
     * @param rNodes Nodal container
     */
    static void StoreFixityData(
        Kratos::span<bool> Span,
        const std::vector<const Variable<double>*>& rDofsVarPointerList,
        const ModelPart::NodesContainerType& rNodes);

    /**
     * @brief Check the validity of the tensor adaptor data
     * @details This method checks if the DOF variables of the tensor adaptor
     * have been addded to as a DOF to each of the nodes in the container
     * @param rDofsVarPointerList List of DOF variables pointers
     * @param rNodes Nodal container
     */
    static void CheckContainer(
        const std::vector<const Variable<double>*>& rDofsVarPointerList,
        const ModelPart::NodesContainerType& rNodes);

    ///@}
};

/// @}
} // namespace Kratos
