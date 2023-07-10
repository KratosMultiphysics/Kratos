// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Franziska Wahl
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
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

/// Utility to add MLS-based master-slave-constraints to negative nodes of intersected elements
/// in order to prevent small cut instabilities.
/// By default, the process also deactivates full negative distance elements.
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) EmbeddedMLSConstraintProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

    using NodesCloudSetType = std::unordered_set<NodeType::Pointer, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using CloudDataVectorType = DenseVector<std::pair<NodeType::Pointer, double>>;

    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of EmbeddedMLSConstraintProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedMLSConstraintProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EmbeddedMLSConstraintProcess() = delete;

    /// Constructor.
    EmbeddedMLSConstraintProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    virtual ~EmbeddedMLSConstraintProcess() = default;

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
            "unknown_variable" : "TEMPERATURE",
            "mls_extension_operator_order" : 1,
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
        return "EmbeddedMLSConstraintProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmbeddedMLSConstraintProcess";
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

    std::string mUnknownVariable;

    std::size_t mMLSExtensionOperatorOrder;

    bool mDeactivateNegativeElements;
    bool mDeactivateIntersectedElements;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateConformingExtensionBasis(
        NodesCloudMapType& rExtensionOperatorMap);

    void ApplyExtensionConstraints(
        NodesCloudMapType& rExtensionOperatorMap);

    void SetInterfaceFlags();

    void ReactivateElementsAndNodes();

    void ModifyDistances();

    bool IsSplit(const GeometryType& rGeometry);

    bool IsNegative(const GeometryType& rGeometry);

    MLSShapeFunctionsFunctionType GetMLSShapeFunctionsFunction();

    void SetNegativeNodeSupportCloud(
        const NodeType& rNegativeNode,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates);

    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin);

    std::size_t GetRequiredNumberOfPoints();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; // Class EmbeddedMLSConstraintProcess

} // namespace Kratos.

