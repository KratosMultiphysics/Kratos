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

#if !defined(KRATOS_ELEMENT_REFINEMENT_PROCESS_H)
#define KRATOS_ELEMENT_REFINEMENT_PROCESS_H

// System includes
#include <string>
#include <vector>
#include <unordered_map>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup MeshingApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(MESHING_APPLICATION) ElementRefinementProcess : public Process
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using Array3D = array_1d<double, 3>;

    template<unsigned int N>
    using SurfaceIndicesArray = array_1d<IndexType, N>;

    template<unsigned int N, unsigned int M>
    using BMatrixNM = BoundedMatrix<double, N, M>;

    using ConditionType = ModelPart::ConditionType;

    using ElementType = ModelPart::ElementType;

    /// Pointer definition of ElementRefinementProcess
    KRATOS_CLASS_POINTER_DEFINITION(ElementRefinementProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ElementRefinementProcess(
        Model& rModel,
        Parameters rParameters);

    ~ElementRefinementProcess() override = default;

    ///@}

    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void InterpolateThreadLocalRefinedMeshFromCoarseElement(const ElementType& rCoarseElement);

    void InterpolateAllRefinedMeshesFromCoarseElement(const ElementType& rCoarseElement);

    void SetEntityIds(const Variable<int>& rVariable);

    void SetConditionParentIds(const Variable<int>& rVariable);

    ModelPart& GetThreadLocalModelPart();

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private classes
    ///@{

    struct ThreadLocalStorage
    {
        ModelPart* pRefinedModelPart;
        std::vector<ModelPart*> mSurfaceModelPartList;
        std::vector<std::vector<std::pair<ConditionType::Pointer, ElementType::Pointer>>> mRefinedModelPartSurfaceParentElements;
    };

    ///@}
    ///@name Private members
    ///@{

    Model& mrModel;
    int mEchoLevel;

    std::string mModelPartName;
    std::string mRefinedModelPartNamePrefix;

    IndexType mRefinementLevel;

    ElementType::Pointer mReferenceCoarseElement;
    std::string mRefinedElementName;
    std::string mRefinedConditionName;

    std::vector<const Variable<double>*> mNodalHistoricalVariablesList;
    std::vector<const Variable<double>*> mNodalNonHistoricalVariablesList;
    std::vector<const Flags*> mNodalFlagsList;
    IndexType mHistoricalVariableInterpolationBufferSize;

    std::vector<ThreadLocalStorage> mThreadLocalStorage;

    std::unordered_map<IndexType, std::unordered_map<IndexType, IndexType>> mElementSurfaceIndicesMap;
    std::vector<std::vector<IndexType>> mSurfacesNodeIndexOrdering;

    ///@}
    ///@name Private  non thread local operations (Must not call in private threads)
    ///@{

    void ComputeSurfaceMap();

    ///@}
    ///@name Private thread local operations (may be called in private threads)
    ///@{

    void CreateSurfaceCondition(
        ThreadLocalStorage& rTLS,
        const IndexType SurfaceIndex,
        Properties::Pointer pProperties,
        ElementType::Pointer pParentElement,
        const IndexType ConditionId,
        const std::vector<IndexType>& rNodeIds);

    void CreateSurfaceAuxiliaries(ThreadLocalStorage& rTLS);

    template<unsigned int TDim, unsigned int TNumNodes, unsigned int TNumSurfaces>
    void RefineElement(
        ThreadLocalStorage& rTLS,
        const IndexType RefinementLevel,
        IndexType& rCurrentNodeId,
        IndexType& rCurrentConditionId,
        IndexType& rCurrentElementId,
        const BMatrixNM<TNumNodes, TNumNodes>& rNodalInterpolationValues,
        const SurfaceIndicesArray<TNumSurfaces>& rSurfaceIndices);

    const std::string GetThreadLocalModelPartName(const IndexType ThreadId) const;

    ThreadLocalStorage& GetThreadLocalStorage();

    ///@}
};

}

///@}

///@}

#endif // KRATOS_ELEMENT_REFINEMENT_PROCESS_H