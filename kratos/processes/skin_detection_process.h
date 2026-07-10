//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes
#include <unordered_set>
#include <unordered_map>

// Project includes
#include "processes/process.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // General geometry type definitions
    using GeometryType = Geometry<Node>;

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The definition of the sizetype
    using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SkinDetectionProcess
 * @ingroup KratosCore
 * @brief An algorithm that looks for neighbour elements in a mesh and creates a submodelpart containing the skin of the disconnected elements (interface elements)
 * @details For that pourpose if builds an unordered map of the surrounding elements and nodes and performs different checks.
 * @tparam TDim The dimension where the problem is computed
 * @author Vicente Mataix Ferrandiz
*/
template<SizeType TDim>
class KRATOS_API(KRATOS_CORE) SkinDetectionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SkinDetectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(SkinDetectionProcess);

    // Weak pointers vectors types
    using NodePointerVector = GlobalPointersVector<Node>;
    using ElementPointerVector = GlobalPointersVector<Element>;

    /// Definition of the vector indexes considered
    using VectorIndexType = std::vector<IndexType>;

    /// Definition of the hasher considered
    using VectorIndexHasherType = VectorIndexHasher<VectorIndexType>;

    /// Definition of the key comparor considered
    using VectorIndexComparorType = VectorIndexComparor<VectorIndexType>;

    /// Define the set considered for element pointers
    using HashSetVectorIntType = std::unordered_set<VectorIndexType, VectorIndexHasherType, VectorIndexComparorType>;

    /// Define the HashSetVectorIntTypeIteratorType iterator type
    using HashSetVectorIntTypeIteratorType = HashSetVectorIntType::iterator;

    /// Define the map considered for face ids
    using HashMapVectorIntType = std::unordered_map<VectorIndexType, VectorIndexType, VectorIndexHasherType, VectorIndexComparorType>;

    /// Define the HashMapVectorIntTypeIteratorType iterator type
    using HashMapVectorIntTypeIteratorType = HashMapVectorIntType::iterator;

    /// Define the map considered for properties ids
    using HashMapVectorIntIdsType = std::unordered_map<VectorIndexType, IndexType, VectorIndexHasherType, VectorIndexComparorType>;

    /// Define the HashMapVectorIntIdsType iterator type
    using HashMapVectorIntIdsTypeIteratorType = HashMapVectorIntIdsType::iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the search of neighbours is performed
     * @param ThisParameters The parameters of configuration
     */
    SkinDetectionProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    virtual ~SkinDetectionProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method executes the algorithm that looks for neighbour nodes and elements in a  mesh of prismatic elements
     */
    void Execute() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SkinDetectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SkinDetectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /** @brief Create auxiliary data structures identifying the element faces on the outer boundary.
     *  @param[out] rInverseFaceMap describes the outer faces of the domain.
     *  @param[out] rPropertiesFaceMap identifies the property of the element the face belongs to.
     */
    void GenerateFaceMaps(
        HashMapVectorIntType& rInverseFaceMap,
        HashMapVectorIntIdsType& rPropertiesFaceMap) const;

    /** @brief Create and prepare the SubModelPart containing the new face conditions.
     *  @return A reference to the new SubModelPart.
     */
    ModelPart& SetUpAuxiliaryModelPart();

    /** @brief Assign new conditions to the target SubModelPart.
     *  @param[out] rAuxiliaryModelPart Empty ModelPart to be filled with the new conditions.
     *  @param[in] rInverseFaceMap auxiliary data structure describing the outer faces of the domain.
     *  @param[in] rPropertiesFaceMap auxiliary data structure identifying the property of the element the face belongs to.
     */
    void FillAuxiliaryModelPart(
        ModelPart& rAuxiliaryModelPart,
        HashMapVectorIntType& rInverseFaceMap,
        HashMapVectorIntIdsType& rPropertiesFaceMap);

    /** @brief Create new Conditions based on the results of the face detection algorithm.
     *  @param[in/out] rMainModelPart Complete ModelPart for the domain.
     *  @param[in/out] rSkinModelPart Target ModelPart that will contain the new conditions.
     *  @param[in] rInverseFaceMap auxiliary data structure describing the outer faces of the domain.
     *  @param[in] rPropertiesFaceMap auxiliary data structure identifying the property of the element the face belongs to.
     *  @param[out] rNodesInTheSkin list of all nodes belonging to the model skin.
     *  @param[in] rConditionName base name for the conditions to be created (number of nodes and dimension will be added dynamically).
     */
    virtual void CreateConditions(
        ModelPart& rMainModelPart,
        ModelPart& rSkinModelPart,
        HashMapVectorIntType& rInverseFaceMap,
        HashMapVectorIntIdsType& rPropertiesFaceMap,
        std::unordered_set<IndexType>& rNodesInTheSkin,
        const std::string& rConditionName
        ) const;

    /** @brief Assign new conditions to additional ModelParts (if requested by user).
     *  @param[in] rAuxiliaryModelPart ModelPart containing the new conditions.
     */
    void SetUpAdditionalSubModelParts(const ModelPart& rAuxiliaryModelPart);

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Protected  Access
    ///@{

    ModelPart& GetModelPart() const;

    Parameters GetSettings() const;

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    /// Protected constructor with modified default settings to be defined by derived class.
    SkinDetectionProcess(ModelPart& rModelPart, Parameters Settings, Parameters DefaultSettings);

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;     /// The main model part
    Parameters mThisParameters; /// The parameters (can be used for general pourposes)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method generates a set with the ids of the nodes of the interface
     * @param[in/out]  rSetNodeIdsInterface The set of ids of the nodes of the interface
     */
    void GenerateSetNodeIdsInterface(std::unordered_set<IndexType>& rSetNodeIdsInterface);

    /**
     * @brief This method filters the nodes of the MPI interface from the 
     * @param[in] rSetNodeIdsInterface The set of ids of the nodes of the interface
     * @param[in/out] rInverseFaceMap auxiliary data structure describing the outer faces of the domain.
     */
    void FilterMPIInterfaceNodes(
        const std::unordered_set<IndexType>& rSetNodeIdsInterface,
        HashMapVectorIntType& rInverseFaceMap
        );

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SkinDetectionProcess& operator=(SkinDetectionProcess const& rOther);

    ///@}

}; // Class SkinDetectionProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<SizeType TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  SkinDetectionProcess<TDim>& rThis);

/// output stream function
template<SizeType TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SkinDetectionProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
