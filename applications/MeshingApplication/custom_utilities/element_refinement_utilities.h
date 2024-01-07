// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_ELEMENT_REFINEMENT_UTILITIES_H)
#define KRATOS_ELEMENT_REFINEMENT_UTILITIES_H

// System includes
#include <vector>
#include <utility>
#include <unordered_map>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/data_communicator.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class ElementRefinementUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using ConditionType = ModelPart::ConditionType;

    using ElementType = ModelPart::ElementType;

    using Array3D = array_1d<double, 3>;

    template<unsigned int N>
    using SurfaceIndicesArray = array_1d<IndexType, N>;

    template<unsigned int N>
    using BVectorN = BoundedVector<double, N>;

    template<unsigned int N>
    using BMatrixNN = BoundedMatrix<double, N, N>;

    template<unsigned int N, unsigned int M>
    using BMatrixNM = BoundedMatrix<double, N, M>;

    ///@}
    ///@name Life Cycle
    ///@{

    ElementRefinementUtilities(
        ModelPart& rModelPart,
        const std::string& rRefinedModelPartName,
        const ElementType& rReferenceElement,
        const IndexType RefinementLevel,
        const std::string rRefinedElementName,
        const std::string rRefinedConditionName,
        const std::vector<std::string>& rHistoricalVariableNamesList,
        const std::vector<std::string>& rNonHistoricalVariableNamesList,
        const std::vector<std::string>& rNodalFlagNamesList);

    ///@}
    ///@name Operations
    ///@{

    void ComputeSurfaceMap();

    void InterpolateToRefinedMeshFromCoarseElement(const ElementType& rCoarseElement);

    void SetIds(const Variable<int>& rVariable);

    void SetConditionParentIds(const Variable<int>& rVariable);

    ModelPart& GetRefinedModelPart();

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    ModelPart& mrModelPart;

    ModelPart* mpRefinedModelPart;

    std::vector<ModelPart*> mSurfaceModelParts;

    std::vector<std::vector<IndexType>> mSurfacesNodeIndexOrdering;

    // <CoarseElementID, <CoarseConditionId, RefinedSurfaceIndex>>
    std::unordered_map<IndexType, std::unordered_map<IndexType, IndexType>> mElementSurfaceIndicesMap;

    std::vector<std::vector<std::pair<ConditionType*, ElementType*>>> mRefinedModelPartSurfaceParentElements;

    const ElementType& mrReferenceElement;

    const IndexType mRefinementLevel;

    const std::string& mrRefinedElementName;

    const std::string& mrRefinedConditionName;

    std::vector<const Variable<double>*> mNodalHistoricalVariablesList;

    std::vector<const Variable<double>*> mNodalNonHistoricalVariablesList;

    std::vector<const Flags*> mNodalFlagsList;

    ///@}
    ///@name Private Operations
    ///@{

    template<unsigned int TDim, unsigned int TNumNodes, unsigned int TNumSurfaces>
    void RefineElement(
        const IndexType RefinementLevel,
        IndexType& rCurrentNodeId,
        IndexType& rCurrentConditionId,
        IndexType& rCurrentElementId,
        const BMatrixNN<TNumNodes>& rNodalInterpolationValues,
        const SurfaceIndicesArray<TNumSurfaces>& rSurfaceIndices);

    void CreateSurfaceCondition(
        const IndexType SurfaceIndex,
        Properties::Pointer pProperties,
        Element* pParentElement,
        const IndexType ConditionId,
        const std::vector<IndexType>& rNodeIds);

    void CreateSurfaceAuxiliaries(
        const std::vector<std::vector<IndexType>>& rSurfaceNodalIndicesList);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_ELEMENT_REFINEMENT_UTILITIES_H
