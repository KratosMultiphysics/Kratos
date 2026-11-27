//  ______  _                _            ______        ______                           _ _                 _
// (_____ \| |              (_)          |  ___ \      |  ___ \         /\              | (_)           _   (_)
//  _____) ) | _  _   _  ___ _  ____  ___| |   | | ____| | _ | | ___   /  \  ____  ____ | |_  ____ ____| |_  _  ___  ____
// |  ____/| || \| | | |/___) |/ ___)/___) |   | |/ _  ) || || |/ _ \ / /\ \|  _ \|  _ \| | |/ ___) _  |  _)| |/ _ \|  _ \
// | |     | | | | |_| |___ | ( (___|___ | |   | ( (/ /| || || | |_| | |__| | | | | | | | | ( (__( ( | | |__| | |_| | | | |
// |_|     |_| |_|\__  (___/|_|\____|___/|_|   |_|\____)_||_||_|\___/|______| ||_/| ||_/|_|_|\____)_||_|\___)_|\___/|_| |_|
//               (____/                                                     |_|   |_|
//
//  License:         BSD License
//                   license: PhysicsNeMoApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class KratosPhysicsNeMoApplication
 * @ingroup PhysicsNeMoApplication
 * @brief Bridge between Kratos Multiphysics and NVIDIA PhysicsNeMo.
 * @details The C++ core of this application is intentionally minimal: all
 * torch/physicsnemo functionality lives in lazily-imported Python modules
 * (python_scripts/), built on top of the core tensor adaptors.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(PHYSICS_NEMO_APPLICATION) KratosPhysicsNeMoApplication
    : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosPhysicsNeMoApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPhysicsNeMoApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosPhysicsNeMoApplication();

    /// Destructor.
    ~KratosPhysicsNeMoApplication() override = default;

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
        return "KratosPhysicsNeMoApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "In KratosPhysicsNeMoApplication" << std::endl;
        KratosApplication::PrintData(rOStream);
    }

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Copy constructor.
    KratosPhysicsNeMoApplication(KratosPhysicsNeMoApplication const& rOther) = delete;

    /// Assignment operator.
    KratosPhysicsNeMoApplication& operator=(KratosPhysicsNeMoApplication const& rOther) = delete;

    ///@}
};

///@}

} // namespace Kratos
