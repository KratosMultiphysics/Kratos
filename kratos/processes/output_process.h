//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


#if !defined(KRATOS_OUTPUT_PROCESS_H_INCLUDED )
#define  KRATOS_OUTPUT_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class OutputProcess
 * @ingroup KratosCore
 * @brief The base class for output processes
 * @details The different output processes must be derived from here
 * @author Philipp Bucher
 */
class OutputProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(OutputProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OutputProcess() : Process() {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    OutputProcess& operator=(OutputProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    virtual bool IsOutputStep()
    {
        return true;
    }

    virtual void PrintOutput()
    {
        KRATOS_ERROR << "Calling the base OutputProcess class PrintOutput. Please implement the PrintOutput in your derived OutputProcess class." << std::endl;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "OutputProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "OutputProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.KratosMultiphysics", Process, OutputProcess)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", Process, OutputProcess)

    ///@}

}; // Class OutputProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const OutputProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_OUTPUT_PROCESS_H_INCLUDED defined
