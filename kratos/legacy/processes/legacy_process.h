//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"

namespace Kratos::Legacy
{

/**
 * @class Process
 * @ingroup KratosCore
 * @brief The base class for all processes in Kratos.
 * @details This is a dummy class to test the legacy namespace and should never be used. Please do not remove this class.
*/
class Process : public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(Process);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Process() : Flags() {}
    explicit Process(const Flags options) : Flags( options ) {}

    /// Destructor.
    ~Process() override {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    virtual void Execute();

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
    std::string Info() const override
    {
        return "Process";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.KratosMultiphysics.Legacy", Process, Process)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All.Legacy", Process, Process)

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Process& operator=(Process const& rOther);

    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Process& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Process& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos::Legacy.
