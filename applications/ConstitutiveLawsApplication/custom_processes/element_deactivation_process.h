// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"

namespace Kratos
{

/**
 * @class ElementDeactivationProcess
 * @ingroup ConstitutiveLawsApplication
 * @brief This class process deactivates elements according to a certain variable threshold
 * We currently suport double and Vector variable types
 * @author Alejandro Cornejo
*/

class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ElementDeactivationProcess
    : public Process
{

    using IndexType = std::size_t;

public:

    KRATOS_CLASS_POINTER_DEFINITION(ElementDeactivationProcess);


    /// Constructor
    ElementDeactivationProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~ElementDeactivationProcess() override = default;

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ElementDeactivationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ElementDeactivationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


protected:

    /// Member Variables

    ModelPart& mrThisModelPart;
    Parameters mThisParameters;
    std::string mVariableName;
    double mThreshold;
    bool mAverageOverIP = true;

private:

    /// Assignment operator.
    ElementDeactivationProcess& operator=(ElementDeactivationProcess const& rOther);

    /// Copy constructor.
    //ElementDeactivationProcess(ElementDeactivationProcess const& rOther);

}; // Class ElementDeactivationProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ElementDeactivationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ElementDeactivationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.