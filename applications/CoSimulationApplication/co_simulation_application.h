//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
class KratosCoSimulationApplication : public KratosApplication {
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
