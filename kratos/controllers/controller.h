//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/registry.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class Controller
 * @ingroup KratosCore
 * @brief The base  class for all Controllers in Kratos.
 * @details Controller is the base class for all controllers and defines the interface for them.
 * Each class must define their own Evaluate method.
 * Controller parameters must be passed at construction time.
  @author Daniel Diez
*/
class Controller
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Controller
    KRATOS_CLASS_POINTER_DEFINITION(Controller);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Controller() noexcept = default;

    /// Destructor.
    virtual ~Controller() = default;

    /// Copy constructor.
    //TODO: Check. It is required by the registry
    Controller(Controller const& rOther) {}

    /// Move constructor
    Controller(Controller&& rOther) noexcept = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Controller& operator=(Controller const& rOther) = delete;

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
    virtual Controller::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const
    {
        KRATOS_ERROR << "Calling base class Create. Please override this method in the corresonding Controller" << std::endl;
        return nullptr;
    }

    /**
     * @brief Checks that input conditions are correct. 
     */
    virtual int Check() const
    {
        return 0;
    }

    /**
     * @brief Using input data, returns bool.
     */
    virtual bool Evaluate()
    {
        KRATOS_ERROR << "Calling base class Evaluate. Please override this method in the corresonding Controller" << std::endl;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    virtual Parameters GetDefaultParameters() const
    {
        KRATOS_ERROR << "Calling the base Controller class GetDefaultParameters. Please implement the GetDefaultParameters in your derived process class." << std::endl;
        const Parameters default_parameters = Parameters(R"({})" );

        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "Controller";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Controller";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
}; // Class Controller

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
std::istream& operator >> (std::istream& rIStream, Controller& rThis);

/// output stream function
std::ostream& operator << (std::ostream& rOStream, const Controller& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
}  // namespace Kratos.
