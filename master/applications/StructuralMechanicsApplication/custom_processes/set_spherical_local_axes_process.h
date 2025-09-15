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

// System includes

// External includes

// Project includes
#include "processes/process.h"

namespace Kratos
{

/**
 * @class SetSphericalLocalAxesProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class set the local axes of the elements according to a SPHERICAL coordinates
 * @author Alejandro Cornejo
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetSphericalLocalAxesProcess
    : public Process
{

public:


    KRATOS_CLASS_POINTER_DEFINITION(SetSphericalLocalAxesProcess);


    /// Constructor
    SetSphericalLocalAxesProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~SetSphericalLocalAxesProcess() override = default;

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
        return "SetSphericalLocalAxesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetSphericalLocalAxesProcess";
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
    SetSphericalLocalAxesProcess& operator=(SetSphericalLocalAxesProcess const& rOther);

    /// Copy constructor.
    //SetSphericalLocalAxesProcess(SetSphericalLocalAxesProcess const& rOther);

}; // Class SetSphericalLocalAxesProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetSphericalLocalAxesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetSphericalLocalAxesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.
