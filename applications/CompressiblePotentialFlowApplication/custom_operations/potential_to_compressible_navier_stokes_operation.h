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
 * @class PotentialToCompressibleNavierStokesOperation
 * @ingroup CompressiblePotentialFlowApplication
 * @brief This operation pass the nodal velocities from Potential model part as initial condition of a compressible 
 * Navier Stokes model part in conservative variables form.
*/
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) PotentialToCompressibleNavierStokesOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Operation
    KRATOS_CLASS_POINTER_DEFINITION(PotentialToCompressibleNavierStokesOperation);

    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics.CompressiblePotentialFlowApplication", PotentialToCompressibleNavierStokesOperation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", PotentialToCompressibleNavierStokesOperation)

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PotentialToCompressibleNavierStokesOperation() : Operation() {}

    PotentialToCompressibleNavierStokesOperation(
        Model& rModel,
        Parameters OperationParameters);

    /// Destructor.
    ~PotentialToCompressibleNavierStokesOperation() override = default;

    /// Copy constructor.
    PotentialToCompressibleNavierStokesOperation(PotentialToCompressibleNavierStokesOperation const& rOther);

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PotentialToCompressibleNavierStokesOperation& operator=(PotentialToCompressibleNavierStokesOperation const& rOther) = delete;


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
}; // Class PotentialToCompressibleNavierStokesOperation

    ///@} addtogroup block

}  // namespace Kratos.
