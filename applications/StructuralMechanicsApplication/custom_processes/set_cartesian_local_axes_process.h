// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

#include "processes/process.h"

namespace Kratos
{

/**
 * @class SetCartesianLocalAxesProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class set the local axes of the elements according to a given set of cartesian axes
 * @author Alejandro Cornejo
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetCartesianLocalAxesProcess
    : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SetCartesianLocalAxesProcess);


    /// Constructor
    SetCartesianLocalAxesProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~SetCartesianLocalAxesProcess() override = default;

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function is designed for being called at the beginning each time step
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SetCartesianLocalAxesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetCartesianLocalAxesProcess";
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
    SetCartesianLocalAxesProcess& operator=(SetCartesianLocalAxesProcess const& rOther);

    /// Copy constructor.
    //SetCartesianLocalAxesProcess(SetCartesianLocalAxesProcess const& rOther);

}; // Class SetCartesianLocalAxesProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetCartesianLocalAxesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetCartesianLocalAxesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.
