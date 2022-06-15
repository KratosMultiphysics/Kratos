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

#ifndef KRATOS_EMBEDDED_CONSTRAINT_H
#define  KRATOS_EMBEDDED_CONSTRAINT_H

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

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

    bool mAvoidZeroDistances;

    bool mDeactivateNegativeElements;
    bool mDeactivateIntersectedElements;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ApplyConstraints();

    void DeactivateElementsAndNodes();

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

#endif // KRATOS_EMBEDDED_CONSTRAINT_UTILITIES_H  defined
