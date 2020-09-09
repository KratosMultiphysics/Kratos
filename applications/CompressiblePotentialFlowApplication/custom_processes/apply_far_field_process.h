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

#if !defined(KRATOS_APPLY_FAR_FIELD_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_FAR_FIELD_PROCESS_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

  /// Auxiliary process to define the wake in 2 dimensional problems.
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) ApplyFarFieldProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyFarFieldProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyFarFieldProcess);

    typedef Node <3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ApplyFarFieldProcess(ModelPart& rModelPart,
                         const double ReferencePotential,
                         const bool InitializeFlowField,
                         const bool PerturbationField);

    /// Copy constructor.
    ApplyFarFieldProcess(ApplyFarFieldProcess const& rOther) = delete;

    /// Destructor.
    ~ApplyFarFieldProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ApplyFarFieldProcess& operator=(ApplyFarFieldProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the ApplyFarFieldProcess algorithms.
    void Execute() override;
    void InitializeFlowField();
    void FindFarthestUpstreamBoundaryNode();
    void AssignFarFieldBoundaryConditions();
    void AssignDirichletFarFieldBoundaryCondition(Geometry<NodeType>& rGeometry);
    void AssignNeumannFarFieldBoundaryCondition(Condition& rCondition);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyFarFieldProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyFarFieldProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;
    NodeType::Pointer mpReferenceNode;
    const double mReferencePotential = 1.0;
    const bool mInitializeFlowField;
    const bool mPerturbationField;
    const array_1d<double,3> mFreeStreamVelocity;
    ///@}
    ///@name Private Operators
    ///@{

    ///@}

}; // Class ApplyFarFieldProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyFarFieldProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyFarFieldProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_DEFINE_2D_WAKE_PROCESS_H_INCLUDED  defined