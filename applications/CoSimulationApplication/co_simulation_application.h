// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

#if !defined(KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED )
#define  KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CO_SIMULATION_APPLICATION) KratosCoSimulationApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosCoSimulationApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosCoSimulationApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosCoSimulationApplication();

    /// Copy constructor.
    KratosCoSimulationApplication(KratosCoSimulationApplication const& rOther) = delete;

    /// Destructor.
    ~KratosCoSimulationApplication() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    KratosCoSimulationApplication& operator=(KratosCoSimulationApplication const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosCoSimulationApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

}; // Class KratosCoSimulationApplication

///@}

}  // namespace Kratos.

#endif // KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED  defined
