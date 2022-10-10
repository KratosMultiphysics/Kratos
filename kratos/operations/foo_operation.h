//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "operations/operation.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class FooOperation
 * @ingroup KratosCore
 * @brief Foo operation to be used while developing
 * @details This fake operation deriving from the base Operation must be removed after development
  @author Ruben Zorrilla
*/
class FooOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FooOperation
    KRATOS_CLASS_POINTER_DEFINITION(FooOperation);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit FooOperation() {}

    /// Model and settings constructor
    explicit FooOperation(
        Model& rModel,
        Parameters ThisParameters) {}

    /// Destructor.
    virtual ~FooOperation() {}

    /// Copy constructor.
    //TODO: Check. It is required by the registry
    FooOperation(FooOperation const& rOther) {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Mmodel and a set of Parameters for the sake of generality
     * @warning Must be overrided in each process implementation
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Operation::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const override
    {
        KRATOS_INFO("FooOperation") << "Create." <<  std::endl;
        return Kratos::make_shared<FooOperation>(rModel, ThisParameters);
    }

    /**
     * @brief Execute method is used to execute the Operation algorithms.
     */
    void Execute() override
    {
        KRATOS_INFO("FooOperation") << "Execute." <<  std::endl;
    }
    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        KRATOS_INFO("FooOperation") << "GetDefaultParameters()." <<  std::endl;
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
        return "FooOperation";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FooOperation";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FooOperation& operator=(FooOperation const& rOther) = delete;

    //TODO: Check. It is required by the registry
    // FooOperation(FooOperation const& rOther) = delete;

    ///@}
}; // Class FooOperation

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, FooOperation& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FooOperation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
}  // namespace Kratos.
