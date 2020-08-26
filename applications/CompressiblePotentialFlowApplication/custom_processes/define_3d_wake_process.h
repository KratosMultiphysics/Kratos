//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

#if !defined(KRATOS_DEFINE_3D_WAKE_PROCESS_H_INCLUDED )
#define  KRATOS_DEFINE_3D_WAKE_PROCESS_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "compressible_potential_flow_application_variables.h"

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

    typedef Node <3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Define3DWakeProcess(ModelPart& rTrailingEdgeModelPart,
                        ModelPart& rBodyModelPart,
                        ModelPart& rStlWakeModelPart,
                        const double Tolerance,
                        const Vector& rWakeNormal);

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

    // The airfoil model part conatining the trailing edge
    ModelPart& mrTrailingEdgeModelPart;
    ModelPart& mrBodyModelPart;
    ModelPart& mrStlWakeModelPart;
    // Tolerance to avoid nodes laying exactly on the wake
    const double mTolerance;
    const Vector& mWakeNormal;
    BoundedVector<double, 3> mWakeDirection;
    BoundedVector<double, 3> mSpanDirection;

    NodeType* mpRightWingTipNode;
    NodeType* mpLeftWingTipNode;

    // Vector to store the trailing edge elements ids
    std::vector<std::size_t> mTrailingEdgeElementsOrderedIds;

    NodeType* mpTrailingEdgeNode;
    BoundedVector<double, 3> mWakeNormalOld;

    ///@}
    ///@name Private Operators
    ///@{
    void InitializeTrailingEdgeSubModelpart() const;

    void InitializeWakeSubModelpart() const;

    void SetWakeAndSpanDirections();

    void MarkTrailingEdgeNodes();

    void ComputeLowerSurfaceNormals() const;

    void MarkWakeElements();

    void CheckIfTrailingEdgeElement(Element& rElement);

    void AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds);

    void MarkKuttaElements();

    void ComputeNodalDistancesToWakeAndLowerSurface(const Element::GeometryType& rGeom, NodeType::Pointer pTrailingEdgeNode, Vector& rNodalDistancesToTe) const;

    void AddWakeNodes() const;

    void SelectWakeTipNodes() const;

    void SelectDoubleTrailingEdgeElements() const;

    void ApplyMPConstraints() const;

    void Print() const;
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