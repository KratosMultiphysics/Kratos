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

#if !defined(KRATOS_DEFINE_3D_WAKE_PROCESS_H_INCLUDED )
#define  KRATOS_DEFINE_3D_WAKE_PROCESS_H_INCLUDED

#include "concurrentqueue/concurrentqueue.h"

namespace Kratos
{

///@name Kratos Classes
///@{

  /// Auxiliary process to define the wake in 2 dimensional problems.
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) Define3DWakeProcess
    : public Process
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

    /// Constructor
    Define3DWakeProcess(ModelPart& rTrailingEdgeModelPart,
                        ModelPart& rBodyModelPart,
                        ModelPart& rStlWakeModelPart,
                        Parameters ThisParameters);

    /// Copy constructor.
    Define3DWakeProcess(Define3DWakeProcess const& rOther) = delete;

    /// Destructor.
    ~Define3DWakeProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Define3DWakeProcess& operator=(Define3DWakeProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// ExecuteInitialize method is used to execute the Define3DWakeProcess algorithms.
    void ExecuteInitialize() override;

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
    ///@name Member Variables
    ///@{

    // The airfoil model part containing the trailing edge
    ModelPart& mrTrailingEdgeModelPart;
    ModelPart& mrBodyModelPart;
    ModelPart& mrStlWakeModelPart;
    // Tolerance to avoid nodes laying exactly on the wake
    double mTolerance;
    BoundedVector<double, 3> mWakeNormal;
    BoundedVector<double, 3> mWakeDirection;
    BoundedVector<double, 3> mSpanDirection;

    bool mSwitchWakeDirection;
    bool mCountElementsNumber;
    bool mWriteElementsIdsToFile;
    bool mShedWakeFromTrailingEdge;
    bool mDecreaseWakeWidthAtTheWingTips;
    int mEchoLevel;

    double mSheddedWakeDistance;
    double mSheddedWakeElementSize;

    BoundedVector<double, 3> mWakeNormalOld;

    ///@}
    ///@name Private Operators
    ///@{
    void InitializeTrailingEdgeSubModelpart() const;

    void InitializeWakeSubModelpart() const;

    void MarkTrailingEdgeNodesAndFindWingtipNodes();

    void ComputeWingLowerSurfaceNormals() const;

    void ComputeAndSaveLocalWakeNormal() const;

    void ShedWakeSurfaceFromTheTrailingEdge() const;

    void DecreaseWakeWidthAtTheWingTips(array_1d<double, 3>& rPoint1,
                                        const array_1d<double, 3>& rPoint2) const;

    void CreateWakeSurfaceNodesAndElements(IndexType& rNode_index,
                                           const array_1d<double, 3>& rCoordinates1,
                                           const array_1d<double, 3>& rCoordinates2,
                                           const array_1d<double, 3>& rCoordinates3,
                                           const array_1d<double, 3>& rCoordinates4,
                                           IndexType& rElement_index,
                                           const Properties::Pointer pElemProp) const;

    std::array<ModelPart::IndexType, 4> CreateWakeSurfaceNodes(
        IndexType& rNode_index,
        const array_1d<double, 3>& rCoordinates1,
        const array_1d<double, 3>& rCoordinates2,
        const array_1d<double, 3>& rCoordinates3,
        const array_1d<double, 3>& rCoordinates4) const;

    double ComputeFaceNormalProjectionToWakeNormal(const array_1d<double, 3>& rCoordinates1,
                                                   const array_1d<double, 3>& rCoordinates2,
                                                   const array_1d<double, 3>& rCoordinates3,
                                                   const array_1d<double, 3>& rCoordinates4) const;

    void CreateWakeSurfaceElements(const double normal_projection,
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

    void CountElementsNumber() const;

    void WriteElementIdsToFile() const;
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

}  // namespace Kratos.

#endif // KRATOS_DEFINE_3D_WAKE_PROCESS_H_INCLUDED  defined