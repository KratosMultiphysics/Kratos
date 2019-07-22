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

#if !defined(KRATOS_CHECK_WAKE_CONDITION_PROCESS_H_INCLUDED )
#define  KRATOS_CHECK_WAKE_CONDITION_PROCESS_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

  /// Auxiliary process to define the wake in 2 dimensional problems.
template<int Dim = 3>
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) CheckWakeConditionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CheckWakeConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(CheckWakeConditionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    CheckWakeConditionProcess(ModelPart& rWakeModelPart, const double Tolerance, const int EchoLevel)
        : Process(),
          mrWakeModelPart(rWakeModelPart),
          mTolerance(Tolerance),
          mEchoLevel(EchoLevel){}

    /// Copy constructor.
    CheckWakeConditionProcess(CheckWakeConditionProcess const& rOther) = delete;

    /// Destructor.
    ~CheckWakeConditionProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CheckWakeConditionProcess& operator=(CheckWakeConditionProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// ExecuteInitialize method is used to execute the CheckWakeConditionProcess algorithms.
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CheckWakeConditionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CheckWakeConditionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    ModelPart& mrWakeModelPart;
    const double mTolerance;
    const int mEchoLevel;

    ///@}

}; // Class CheckWakeConditionProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CheckWakeConditionProcess<>& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CheckWakeConditionProcess<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_CHECK_WAKE_CONDITION_PROCESS_H_INCLUDED  defined