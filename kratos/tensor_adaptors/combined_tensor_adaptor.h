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