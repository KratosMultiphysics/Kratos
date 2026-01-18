//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//                   

#pragma once

// System includes

// External includes

// Project includes
#include "operations/operation.h"

namespace Kratos
{
  ///@addtogroup CompressiblePotentialFlowApplication
  ///@{

///@name Kratos Classes
///@{

/**
 * @author Marco Antonio Zuñiga Perez
 * @class PrimalToAdjointOperation
 * @ingroup CompressiblePotentialFlowApplication
 * @brief This operation passes the selected variables from an origin model part 
 * to a destination model part.
*/
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) PrimalToAdjointOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Operation
    KRATOS_CLASS_POINTER_DEFINITION(PrimalToAdjointOperation);

    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics.CompressiblePotentialFlowApplication", Operation, PrimalToAdjointOperation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", Operation, PrimalToAdjointOperation)

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PrimalToAdjointOperation() : Operation() {}

    PrimalToAdjointOperation(
        Model& rModel,
        Parameters OperationParameters);

    /// Destructor.
    ~PrimalToAdjointOperation() override = default;

    /// Copy constructor.
    PrimalToAdjointOperation(PrimalToAdjointOperation const& rOther);

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PrimalToAdjointOperation& operator=(PrimalToAdjointOperation const& rOther) = delete;


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the operation
     * @details We consider as input a Mmodel and a set of Parameters for the sake of generality
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Operation::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /**
     * @brief Execute method is used to execute the Operation algorithms.
     */
    void Execute() override;

    ///@}
private:
    ///@}
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;

    ///@}
}; // Class PrimalToAdjointOperation

    ///@} addtogroup block

}  // namespace Kratos.
