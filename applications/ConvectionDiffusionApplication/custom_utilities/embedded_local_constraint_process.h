// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

#ifndef KRATOS_EMBEDDED_LOCAL_CONSTRAINT_H
#define  KRATOS_EMBEDDED_LOCAL_CONSTRAINT_H

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "utilities/mls_shape_functions_utility.h"

// Application includes

namespace Kratos
{

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

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using CloudDataVectorType = DenseVector<std::pair<NodeType::Pointer, double>>;

    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using NodesOffsetMapType = std::unordered_map<NodeType::Pointer, double, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

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
        Parameters ThisParameters);

    /// Destructor.
    virtual ~EmbeddedLocalConstraintProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "apply_to_all_negative_cut_nodes" : false,
            "use_mls_shape_functions" : true,
            "include_intersection_points" : true,
            "avoid_zero_distances" : true,
            "deactivate_negative_elements" : true,
            "deactivate_intersected_elements" : false
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
    virtual void PrintData(std::ostream& rOStream) const override
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

    bool mApplyToAllNegativeCutNodes;
    bool mUseMLSShapeFunctions;
    bool mIncludeIntersectionPoints;

    bool mAvoidZeroDistances;

    bool mDeactivateNegativeElements;
    bool mDeactivateIntersectedElements;

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

    void AddMLSNodeClouds(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap,
        std::vector<NodeType::Pointer> neg_nodes_element,
        std::vector<NodeType::Pointer> pos_nodes_element);

    void AddMLSNodeCloudsIncludingBC(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap,
        std::vector<NodeType::Pointer> neg_nodes_element,
        std::vector<NodeType::Pointer> pos_nodes_element);

    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin);

    void ApplyConstraints(
        NodesCloudMapType& rCloudsMap,
        NodesOffsetMapType& rOffsetsMap);

    void DeactivateElementsAndNodes();

    void ModifyDistances();

    bool IsSplit(const GeometryType& rGeometry);

    bool IsSmallCut(const GeometryType& rGeometry);

    bool IsNegative(const GeometryType& rGeometry);

    MLSShapeFunctionsFunctionType GetMLSShapeFunctionsFunction();

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

#endif // KRATOS_EMBEDDED_LOCAL_CONSTRAINT_H  defined
