//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define_registry.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/registry.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class Operation
 * @ingroup KratosCore
 * @brief The base class for all operations in Kratos.
 * @details The operation is the base class for all operations and defines the interface for them.
 * Execute method is used to execute all the Operation algorithms. Operation parameters must be passed at construction time.
  @author Carlos Roig
  @author Ruben Zorrilla
*/
class Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Operation
    KRATOS_CLASS_POINTER_DEFINITION(Operation);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit Operation() = default;

    /// Destructor.
    virtual ~Operation() {}

    /// Copy constructor.
    //TODO: Check. It is required by the registry
    Operation(Operation const& rOther) {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Operation& operator=(Operation const& rOther) = delete;

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
    virtual Operation::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const
    {
        KRATOS_ERROR << "Calling base class Create. Please override this method in the corresonding Operation" << std::endl;
        return nullptr;
    }

    /**
     * @brief Execute method is used to execute the Operation algorithms.
     */
    virtual void Execute()
    {
        KRATOS_ERROR << "Calling base class Execute. Please override this method in the corresonding Operation" << std::endl;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    virtual const Parameters GetDefaultParameters() const
    {
        KRATOS_ERROR << "Calling the base Operation class GetDefaultParameters. Please implement the GetDefaultParameters in your derived process class." << std::endl;
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
        return "Operation";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Operation";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics", Operation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", Operation)

    ///@}
}; // Class Operation

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, Operation& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const Operation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
}  // namespace Kratos.
