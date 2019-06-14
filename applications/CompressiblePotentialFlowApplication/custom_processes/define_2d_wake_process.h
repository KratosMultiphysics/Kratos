//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

#if !defined(KRATOS_DEFINE_2D_WAKE_PROCESS_H_INCLUDED )
#define  KRATOS_DEFINE_2D_WAKE_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

  /// Auxiliary process to define the wake in 2 dimensional problems.
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) Define2DWakeProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Define2DWakeProcess
    KRATOS_CLASS_POINTER_DEFINITION(Define2DWakeProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Define2DWakeProcess(
        ModelPart& rModelPart
        ) : Process() ,
            mrModelPart(rModelPart)
    {
        std::cout << "Entering Constructor" << std::endl;
    }

    /// Copy constructor.
    Define2DWakeProcess(Define2DWakeProcess const& rOther) = delete;

    /// Destructor.
    ~Define2DWakeProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Define2DWakeProcess& operator=(Define2DWakeProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the Define2DWakeProcess algorithms.
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Define2DWakeProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Define2DWakeProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The main model part where the wake elements will be defined

    ///@}

}; // Class Define2DWakeProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Define2DWakeProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Define2DWakeProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_DEFINE_2D_WAKE_PROCESS_H_INCLUDED  defined
