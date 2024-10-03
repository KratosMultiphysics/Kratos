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
 * @class SetAutomatedInitialVariableProcessCartesian
 * @ingroup StructuralMechanicsApplication
 * @brief This class automotes the creation of the variables INITIAL_STRAIN_VECTOR and INITIAL_STRESS_VECTOR using tables imported from csv files
 * @author Luis Antonio Goncalves Junior
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetAutomatedInitialVariableProcessCartesian
    : public Process
{

public:

    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_CLASS_POINTER_DEFINITION(SetAutomatedInitialVariableProcessCartesian);


    /// Constructor
    SetAutomatedInitialVariableProcessCartesian(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~SetAutomatedInitialVariableProcessCartesian() override = default;

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
        return "SetAutomatedInitialVariableProcessCartesian";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetAutomatedInitialVariableProcessCartesian";
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
    SetAutomatedInitialVariableProcessCartesian& operator=(SetAutomatedInitialVariableProcessCartesian const& rOther);

    /// Copy constructor.
    //SetAutomatedInitialVariableProcessCartesian(SetAutomatedInitialVariableProcessCartesian const& rOther);

}; // Class SetAutomatedInitialVariableProcessCartesian

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetAutomatedInitialVariableProcessCartesian& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetAutomatedInitialVariableProcessCartesian& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

