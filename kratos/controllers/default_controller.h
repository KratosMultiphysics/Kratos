//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "controllers/controller.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class Default controller
 * @ingroup KratosCore
 * @brief A controller which is always evaluated to be true.
 * @details This controller always evaluates to true.
 * @author Suneth Warnakulasuriya
*/
class DefaultController: public Controller
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DefaultController
    KRATOS_CLASS_POINTER_DEFINITION(DefaultController);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DefaultController() noexcept = default;

    /// Destructor.
    ~DefaultController() override = default;

    /// Copy constructor.
    //TODO: Check. It is required by the registry
    DefaultController(DefaultController const& rOther) = default;

    /// Move constructor
    DefaultController(DefaultController&& rOther) noexcept = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    DefaultController& operator=(DefaultController const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the controller
     * @details We consider as input, a Model and a set of Parameters for the sake of generality
     * @warning Must be overrided in controller implementation
     * @param rModel The model to be considered
     * @param ThisParameters The configuration parameters
     */
    Controller::Pointer Create(
        Model&,
        Parameters) const override
    {
        return Kratos::make_shared<DefaultController>();
    }

    /**
     * @brief Using input data, returns bool.
     */
    bool Evaluate() const override
    {
        return true;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DefaultController";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DefaultController";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
}; // Class DefaultController

}  // namespace Kratos.
