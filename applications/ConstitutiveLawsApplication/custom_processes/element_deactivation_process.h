// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
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
 * @class ElementDeactivationProcess
 * @ingroup ConstitutiveLawsApplication
 * @brief This class process deactivates elements according to a certain variable threshold
 * We currently suport double and Vector variable types
 * @author Alejandro Cornejo
*/
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ElementDeactivationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{
    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ElementDeactivationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ElementDeactivationProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor
    ~ElementDeactivationProcess() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

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

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mrThisModelPart;
    Parameters mThisParameters;
    std::string mVariableName;
    double mThreshold;
    bool mAverageOverIP = true;

    ///@}
private:
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ElementDeactivationProcess& operator=(ElementDeactivationProcess const& rOther);

    /// Copy constructor.
    //ElementDeactivationProcess(ElementDeactivationProcess const& rOther);
    
    ///@}
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