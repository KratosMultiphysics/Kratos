//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Marc Nunez
//

#pragma once

// External includes
#include "concurrentqueue/concurrentqueue.h"

// Project includes
#include "includes/kratos_parameters.h"
#include "operations/operation.h"

namespace Kratos
{
///@addtogroup CompressiblePotentialFlowApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class Define3DWakeOperation
 * @ingroup CompressiblePotentialFlowApplication
 * @brief This operation define the wake in 3d problems.
 * @details For a given model, this operation load a modeled stl wake or creates one with minimal 
 * geomtrical data, and uses it to select the wake and kutta elements. 
 * @authors Inigo Lopez and Marc Nunez 
 */
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) Define3DWakeOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Define3DWakeOperation
    KRATOS_CLASS_POINTER_DEFINITION(Define3DWakeOperation);

    typedef Node NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Define3DWakeOperation() : Operation() {};

    /// Constructor

    /**
     * @brief Constructor with Kratos parameters and Model container
     * @param rModel The Model container
     * @param rParameters Kratos parameters encapsulating the settings
     */
    Define3DWakeOperation(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~Define3DWakeOperation() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Define3DWakeOperation& operator=(Define3DWakeOperation const& rOther) = delete;

    /// Copy constructor.
    Define3DWakeOperation(Define3DWakeOperation const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Operation::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<Define3DWakeOperation>(rModel, ThisParameters);
    }

    /**
     * @brief Execute method is used to execute the Operation algorithms.
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;
    
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "Define3DWakeOperation";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Define3DWakeOperation";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{
        
    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics.CompressiblePotentialFlowApplication", Operation, Define3DWakeOperation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", Operation, Define3DWakeOperation)
    
    ///@}
    ///@name Member Variables
    ///@{

    // Model parts
    ModelPart* mpTrailingEdgeModelPart = nullptr;
    ModelPart* mpBodyModelPart = nullptr;
    ModelPart* mpUpperSurfaceModelPart = nullptr;
    ModelPart* mpLowerSurfaceModelPart = nullptr;
    ModelPart* mpRootPointsModelPart = nullptr;
    ModelPart* mpTipPointsModelPart = nullptr;
    ModelPart* mpBluntTESurfaceModelPart = nullptr;
    ModelPart* mpWakeModelPart = nullptr;
    ModelPart* mpRootModelPart = nullptr;

    // Tolerances 
    double mWakeDistanceTolerance;                  // To avoid nodes laying exactly on the wake 

    BoundedVector<double, 3> mWakeNormal;
    BoundedVector<double, 3> mWakeDirection;
    BoundedVector<double, 3> mSpanDirection;
    BoundedVector<double, 3> mWakeTraslationDirection;

    bool mShedWakeFromTrailingEdge;
    bool mVisualizeWakeVTK;

    int mEchoLevel;
    
    double mShedWakeLength;
    double mShedWakeElementSize;
    double mShedGrowFactor;
    double mShedProjectionRootEdge;

    std::unordered_set<IndexType> mBluntIds; 
    
    std::filesystem::path mWakeSTLFileName;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method initializes nodal and elemental variables such as 
     * UPPER_SURFACE, LOWER_SURFACE, TRAILING_EDGE, WAKE, WAKE_DISTANCE,
     * saves the WAKE_NORMAL and computes the ake_direction and span direction.
     */
    void InitializeVariables();

    /**
     * @brief This method initializes the trailing edge submodelpart.
     */
    void InitializeTrailingEdgeSubModelpart() const;

    /**
     * @brief This method initializes the wake submodelpart.
     */
    void InitializeWakeSubModelpart() const;

    /**
     * @brief This method marks the trailing edge nodes and, when defined
     * marks all blunt nodes too, where the airfoil does not end in a thin edge,
     * but in a flat termination.
     */
    void MarkTrailingEdgeAndBluntNodes();

    /**
     * @brief This method adds conditions to the trailing edge model part and marks or
     * finds the tip and root nodes, if necessary.
     */
    void AddTrailingEdgeConditionsAndFindRootAndTipNodes() const;

    /**
     * @brief This method computes the wing lower surface normals and marks the upper
     * and lower surfaces. The wing lower surface normals are used later in
     * RecomputeComputeNodalDistancesToWakeOrWingLowerSurface inside the
     * MarkKuttaElements function to check whether nodes are above or below the wake
     * Upper and Lower surfaces model parts are used, if provided.
     */
    void ComputeWingLowerSurfaceNormals() const;

    /**
     * @brief This method computes the local wake normal at each trailing edge node
     * by avaraging the local wake normals of the surrounding conditions.
     */
    void ComputeAndSaveLocalWakeNormal() const;

    /**
     * @brief This method creates the wake surface automatically by shedding it from the
     * trailing edge/s in the direction of the free stream velocity (mWakeDirection).
     * The user can decide how much distance is to be shedded, the element size
     * of the wake surface in the wake direction and the grow factor of the element size. 
     * Note that the element size in span direction is predetermined by the size of the 
     * conditions constituting the trailing edge.
     */
    void ShedWakeSurfaceFromTheTrailingEdge() const;

    /**
     * @brief This method allows to load a stl model of the wake. 
     */
    void LoadSTL() const;

    /**
     * @brief This method allows to move the wake_model_part. 
     */
    void MoveWakeModelPart() const;

    /**
     * @brief This method creates a vtk output model of the wake. 
     */
    void VisualizeWake() const;

    /**
     * @brief This method creates the wake surface nodes and elements.
     */
    void CreateWakeSurfaceNodesAndElements(ModelPart& ModelPart,
                                           IndexType& rNode_index,
                                           const array_1d<double, 3>& rCoordinates1,
                                           const array_1d<double, 3>& rCoordinates2,
                                           const array_1d<double, 3>& rCoordinates3,
                                           const array_1d<double, 3>& rCoordinates4,
                                           IndexType& rElement_index,
                                           const Properties::Pointer pElemProp) const;

    /**
     * @brief This method creates the wake surface nodes. 
     */
    std::array<ModelPart::IndexType, 4> CreateWakeSurfaceNodes(
        ModelPart& ModelPart,
        IndexType& rNode_index,
        const array_1d<double, 3>& rCoordinates1,
        const array_1d<double, 3>& rCoordinates2,
        const array_1d<double, 3>& rCoordinates3,
        const array_1d<double, 3>& rCoordinates4) const;

    /**
     * @brief This method computes the face normal projected to the wake normal.
     */
    double ComputeFaceNormalProjectionToWakeNormal(const array_1d<double, 3>& rCoordinates1,
                                                   const array_1d<double, 3>& rCoordinates2,
                                                   const array_1d<double, 3>& rCoordinates3,
                                                   const array_1d<double, 3>& rCoordinates4) const;

    /**
     * @brief This method creates wake surface elements. 
     */
    void CreateWakeSurfaceElements(ModelPart& ModelPart,
                                   const double normal_projection,
                                   IndexType& rElement_index,
                                   const std::array<ModelPart::IndexType, 4>& rNodes_ids,
                                   const Properties::Pointer pElemProp) const;

    /**
     * @brief This method checks which elements are cut by the wake and marks them as
     * wake elements.
     */
    void MarkWakeElements() const;

    /**
     * @brief This method checks if the element is touching the trailing edge.
     */
    void CheckIfTrailingEdgeElement(Element& rElement,
                                    const Geometry<NodeType>& rGeometry,
                                    moodycamel::ConcurrentQueue<std::size_t> & rTrailingEdgeElementsOrderedIds) const;

    /**
     * @brief This method adds the trailing edge elements in the trailing_edge_sub_model_part. 
     */
    void AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds,
                                        std::vector<std::size_t>& rTrailingEdgeElementsOrderedIds) const;

    /**
     * @brief This method finds the closest trailing edge node to the given point.
     */
    void FindClosestTrailingEdgeNode(NodeType::Pointer& pClosest_te_node,
        const array_1d<double, 3>& rPoint) const;
        
    /**
     * @brief This method recomputes the wake distances from the nodes belonging to the
     * elements. These distances are used later to decide which elements are KUTTA,
     * which are WAKE (STRUCTURE), and which are NORMAL.
     */
    void RecomputeNodalDistancesToWakeOrWingLowerSurface() const;

    /**
     * @brief This method recomputes the distances from the given node to the wake or to
     * the wing lower surface, depending whether the node is behind or in front of
     * the trailing edge according to the wake direction (free stream velocity
     * direction). 
     */
    void RecomputeDistance(NodeType::Pointer& pClosest_te_node,
                           NodeType& rNode) const;

    /**
     * @brief This method This function selects the kutta elements. Kutta elements are touching the 
     * trailing edge from below. 
     */
    void MarkKuttaElements() const;

    /**
     * @brief This method returns the number of trailing edge nodes.
     */
    unsigned int CountNumberOfTrailindEdgeNodesInElement(const Geometry<NodeType>& rGeometry) const;

    /**
     * @brief This method returns the number of non trailing edge nodes with positive and
     * negative distance.
     */
    void CountNumberOfPositiveAndNegativeDistances(
        const Geometry<NodeType>& rGeometry,
        unsigned int& number_of_nodes_with_negative_distance,
        unsigned int& number_of_nodes_with_positive_distance) const;

    /**
     * @brief This method selects the element type (WAKE, KUTTA or NORMAL)
     */
    void SelectElementType(Element& rElement,
                           const Geometry<NodeType>& rGeometry,
                           const unsigned int number_of_te_nodes,
                           const unsigned int number_of_nodes_with_negative_distance,
                           const unsigned int number_of_nodes_with_positive_distance) const;

    /**
     * @brief This method saves the wake normal in the elements.
     */
    void SaveLocalWakeNormalInElements() const;

    /**
     * @brief This method adds nodes to the wake model part.
     */
    void AddWakeNodesToWakeModelPart() const;

    ///@}

}; // Class Define3DWakeOperation

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Define3DWakeOperation& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Define3DWakeOperation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.