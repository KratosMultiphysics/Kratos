//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_FORMULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_FORMULATION_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief This class is extending standard Process interface
 *
 * This class extends standard Process interface to have methods to be
 * called before and after coupling step solve in Formulations. This is inherited
 * from Process so, same processes can be used as normal Processes, or if required
 * they can be used in Formulation Process lists.
 *
 */
class RansFormulationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansFormulationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansFormulationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RansFormulationProcess() : Process()
    {
    }

    /// Destructor.
    ~RansFormulationProcess() override = default;

    RansFormulationProcess& operator=(
        RansFormulationProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief The method is called before executing a coupling solve step
     *
     * This method is called before executing Formulation.SolveCouplingStep()
     *
     */
    virtual void ExecuteBeforeCouplingSolveStep()
    {
    }

    /**
     * @brief The method is called after executing a coupling solve step
     *
     * This method is called after executing Formulation.SolveCouplingStep()
     *
     */
    virtual void ExecuteAfterCouplingSolveStep()
    {
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "RansFormulationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansFormulationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

}; // Class RansFormulationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(
    std::istream& rIStream,
    RansFormulationProcess& rThis);

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansFormulationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_RANS_FORMULATION_PROCESS_H_INCLUDED  defined
