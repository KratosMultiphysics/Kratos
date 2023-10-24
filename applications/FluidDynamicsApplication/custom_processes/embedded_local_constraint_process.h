//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//
//

#ifndef KRATOS_EMBEDDED_LOCAL_CONSTRAINT_PROCESS_H
#define  KRATOS_EMBEDDED_LOCAL_CONSTRAINT_PROCESS_H

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Utility to add a constraint to negative nodes of intersected elements
/// in order to prevent small cut instabilities.
/// By default, the process also deactivates full negative distance elements.
class EmbeddedLocalConstraintProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    using DofType = Dof<double>;

    using DofPointerVectorType = std::vector<DofType::Pointer>;

    using MatrixType = Matrix;

    using VectorType = Vector;

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using CloudDataVectorType = DenseVector<std::pair<NodeType::Pointer, double>>;

    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using NodesOffsetMapType = std::unordered_map<NodeType::Pointer, double, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of EmbeddedLocalConstraintProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedLocalConstraintProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedLocalConstraintProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    virtual ~EmbeddedLocalConstraintProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name"                   : "",
            "check_at_each_time_step"           : true,
            "apply_to_all_negative_cut_nodes"   : true,
            "include_intersection_points"       : true,
            "avoid_zero_distances"              : true,
            "deactivate_full_negative_elements" : true,
            "slip_length"                       : 0.0
        })" );

        return default_parameters;
    }

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
    virtual std::string Info() const override
    {
        return "EmbeddedLocalConstraintProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmbeddedLocalConstraintProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}

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

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;
    const std::array<std::string,4> mComponents = {"PRESSURE","VELOCITY_X","VELOCITY_Y","VELOCITY_Z"};

    bool mConstraintsAreCalculated;
    bool mCheckAtEachStep;

    bool mApplyToAllNegativeCutNodes;
    bool mIncludeIntersectionPoints;
    double mSlipLength;

    bool mAvoidZeroDistances;

    bool mNegElemDeactivation;

    // Count number of positive small cuts
    std::size_t mNSmallCutPos;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateNodeClouds(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap);

    void AddAveragedNodeClouds(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap,
        std::vector<NodeType::Pointer> neg_nodes_element,
        std::vector<NodeType::Pointer> pos_nodes_element);

    void AddAveragedNodeCloudsIncludingBC(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap,
        std::vector<NodeType::Pointer> neg_nodes_element,
        std::vector<NodeType::Pointer> pos_nodes_element);

    void ApplyConstraints(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap);

    void DeactivateFullNegativeElements();

    void RecoverDeactivationPreviousState();

    void ModifyDistances();

    bool IsSplit(const GeometryType& rGeometry);

    bool IsSmallCut(const GeometryType& rGeometry);

    bool IsNegative(const GeometryType& rGeometry);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    EmbeddedLocalConstraintProcess() = delete;

    /// Assignment operator.
    // EmbeddedLocalConstraintProcess& operator=(EmbeddedLocalConstraintProcess const& rOther);

    /// Copy constructor.
    //EmbeddedLocalConstraintProcess(EmbeddedLocalConstraintProcess const& rOther);

    ///@}

}; // Class EmbeddedLocalConstraintProcess

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_LOCAL_CONSTRAINT_PROCESS_H  defined
