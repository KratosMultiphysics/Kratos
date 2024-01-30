// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Luis Antonio Goncalves Junior
//                   Alejandro Cornejo
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"

namespace Kratos
{

/**
 * @class SetAutomatedInitialVariableProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class automotes the creation of the variables INITIAL_STRAIN_VECTOR and INITIAL_STRESS_VECTOR using tables imported from csv files
 * @author Luis Antonio Goncalves Junior
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetAutomatedInitialVariableProcess
    : public Process
{

public:

    static constexpr double tolerance         = 1.0e-6;
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_CLASS_POINTER_DEFINITION(SetAutomatedInitialVariableProcess);


    /// Constructor
    SetAutomatedInitialVariableProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~SetAutomatedInitialVariableProcess() override = default;

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SetAutomatedInitialVariableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetAutomatedInitialVariableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


protected:

    /// Member Variables

    ModelPart& mrThisModelPart;
    Parameters mThisParameters;

private:

    /// Assignment operator.
    SetAutomatedInitialVariableProcess& operator=(SetAutomatedInitialVariableProcess const& rOther);

    /// Copy constructor.
    //SetAutomatedInitialVariableProcess(SetAutomatedInitialVariableProcess const& rOther);

}; // Class SetAutomatedInitialVariableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetAutomatedInitialVariableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetAutomatedInitialVariableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

