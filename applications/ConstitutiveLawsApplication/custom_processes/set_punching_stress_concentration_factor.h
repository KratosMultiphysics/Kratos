// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: kratos/license.txt
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
 * @class SetPunchingStressConcentrationFactor
 * @ingroup ConstitutiveLawsApplication
 * @brief This class sets the hole diameter d and the specimen width H to compute the punching stress concentration factor K_t
 * @author Luis Antonio Goncalves Junior
*/

class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SetPunchingStressConcentrationFactor
    : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SetPunchingStressConcentrationFactor);


    /// Constructor
    SetPunchingStressConcentrationFactor(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~SetPunchingStressConcentrationFactor() override = default;

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SetPunchingStressConcentrationFactor";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetPunchingStressConcentrationFactor";
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
    SetPunchingStressConcentrationFactor& operator=(SetPunchingStressConcentrationFactor const& rOther);

    /// Copy constructor.
    //SetPunchingStressConcentrationFactor(SetPunchingStressConcentrationFactor const& rOther);

}; // Class SetPunchingStressConcentrationFactor

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetPunchingStressConcentrationFactor& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetPunchingStressConcentrationFactor& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.