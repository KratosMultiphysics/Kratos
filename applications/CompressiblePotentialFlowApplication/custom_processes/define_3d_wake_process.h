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

#ifndef KRATOS_DEFINE_3D_WAKE_PROCESS_H
#define KRATOS_DEFINE_3D_WAKE_PROCESS_H

// External includes
#include "concurrentqueue/concurrentqueue.h"

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup CompressiblePotentialFlowApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class Define3DWakeProcess
 * @ingroup CompressiblePotentialFlowApplication
 * @brief This process define the wake in 3d problems.
 * @details For a given model, this process load a modeled stl wake or creates one with minimal 
 * geomtrical data, and uses it to select the wake and kutta elements. 
 * @authors Inigo Lopez and Marc Nunez 
 */
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) Define3DWakeProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Define3DWakeProcess
    KRATOS_CLASS_POINTER_DEFINITION(Define3DWakeProcess);

    typedef Node NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Define3DWakeProcess() : Process() {};

    /// Constructor

    /**
     * @brief Constructor with Kratos parameters and Model container
     * @param rModel The Model container
     * @param rParameters Kratos parameters encapsulating the settings
     */
    Define3DWakeProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~Define3DWakeProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Define3DWakeProcess& operator=(Define3DWakeProcess const& rOther) = delete;

    /// Copy constructor.
    Define3DWakeProcess(Define3DWakeProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) override
    {
        return Kratos::make_shared<Define3DWakeProcess>(rModel, ThisParameters);
    }

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;
    
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Define3DWakeProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Define3DWakeProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{
        
    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.KratosMultiphysics.CompressiblePotentialFlowApplication", Process, Define3DWakeProcess)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", Process, Define3DWakeProcess)
    
    ///@}
    ///@name Member Variables
    ///@{

    // Model parts
    ModelPart* mrTrailingEdgeModelPart = nullptr;
    ModelPart* mrBodyModelPart = nullptr;
    ModelPart* mrUpperSurfaceModelPart = nullptr;
    ModelPart* mrLowerSurfaceModelPart = nullptr;
    ModelPart* mrRootPointsModelPart = nullptr;
    ModelPart* mrTipPointsModelPart = nullptr;
    ModelPart* mrBluntTESurfaceModelPart = nullptr;
    ModelPart* mrWakeModelPart = nullptr;
    ModelPart* mrRootModelPart = nullptr;

    // Tolerances 
    double mWakeDistanceTolerance;                  // To avoid nodes laying exactly on the wake 
    double mCheckWakeConditionTolerance;            // To check if wake conditions are fulfilled

    BoundedVector<double, 3> mWakeNormal;
    BoundedVector<double, 3> mWakeDirection;
    BoundedVector<double, 3> mSpanDirection;
    BoundedVector<double, 3> mWakedrTraslation;

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

    void InitializeVariables();

    void InitializeTrailingEdgeSubModelpart() const;

    void InitializeWakeSubModelpart() const;

    void MarkTrailingEdgeAndBluntNodes();

    void AddTrailingEdgeConditionsAndFindRootAndTipNodes() const;

    void ComputeWingLowerSurfaceNormals() const;

    void ComputeAndSaveLocalWakeNormal() const;

    void ShedWakeSurfaceFromTheTrailingEdge() const;

    void LoadSTL() const;

    void MoveWakeModelPart() const;

    void VisualizeWake() const;

    void DecreaseWakeWidthAtTheWingTips(array_1d<double, 3>& rPoint1,
                                        const array_1d<double, 3>& rPoint2) const;

    void CreateWakeSurfaceNodesAndElements(ModelPart& ModelPart,
                                           IndexType& rNode_index,
                                           const array_1d<double, 3>& rCoordinates1,
                                           const array_1d<double, 3>& rCoordinates2,
                                           const array_1d<double, 3>& rCoordinates3,
                                           const array_1d<double, 3>& rCoordinates4,
                                           IndexType& rElement_index,
                                           const Properties::Pointer pElemProp) const;

    std::array<ModelPart::IndexType, 4> CreateWakeSurfaceNodes(
        ModelPart& ModelPart,
        IndexType& rNode_index,
        const array_1d<double, 3>& rCoordinates1,
        const array_1d<double, 3>& rCoordinates2,
        const array_1d<double, 3>& rCoordinates3,
        const array_1d<double, 3>& rCoordinates4) const;

    double ComputeFaceNormalProjectionToWakeNormal(const array_1d<double, 3>& rCoordinates1,
                                                   const array_1d<double, 3>& rCoordinates2,
                                                   const array_1d<double, 3>& rCoordinates3,
                                                   const array_1d<double, 3>& rCoordinates4) const;

    void CreateWakeSurfaceElements(ModelPart& ModelPart,
                                   const double normal_projection,
                                   IndexType& rElement_index,
                                   const std::array<ModelPart::IndexType, 4>& rNodes_ids,
                                   const Properties::Pointer pElemProp) const;

    void MarkWakeElements() const;

    void CheckIfTrailingEdgeElement(Element& rElement,
                                    const Geometry<NodeType>& rGeometry,
                                    moodycamel::ConcurrentQueue<std::size_t> & rTrailingEdgeElementsOrderedIds) const;

    void AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds,
                                        std::vector<std::size_t>& rTrailingEdgeElementsOrderedIds) const;

    void FindClosestTrailingEdgeNode(NodeType::Pointer& pClosest_te_node,
                                     const array_1d<double, 3>& rPoint) const;

    void RecomputeNodalDistancesToWakeOrWingLowerSurface() const;

    void RecomputeDistance(NodeType::Pointer& pClosest_te_node,
                           NodeType& rNode) const;

    void MarkKuttaElements() const;

    unsigned int CountNumberOfTrailindEdgeNodesInElement(const Geometry<NodeType>& rGeometry) const;

    void CountNumberOfPositiveAndNegativeDistances(
        const Geometry<NodeType>& rGeometry,
        unsigned int& number_of_nodes_with_negative_distance,
        unsigned int& number_of_nodes_with_positive_distance) const;

    void SelectElementType(Element& rElement,
                           const Geometry<NodeType>& rGeometry,
                           const unsigned int number_of_te_nodes,
                           const unsigned int number_of_nodes_with_negative_distance,
                           const unsigned int number_of_nodes_with_positive_distance) const;

    void SaveLocalWakeNormalInElements() const;

    void AddWakeNodesToWakeModelPart() const;

    ///@}

}; // Class Define3DWakeProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Define3DWakeProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Define3DWakeProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DEFINE_3D_WAKE_PROCESS_H