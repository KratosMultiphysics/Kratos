// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Lorenzo Gracia
//                   Aron Noordam
//                   Vahid Galavi
//                   Marjan Fathian
//                   Ruben Zorrilla
//

#pragma once

// Project includes
#include "operations/operation.h"

namespace Kratos
{

///@addtogroup GeoMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @author Lorenzo Gracia
 * @author Aron Noordam
 * @author Vahid Galavi
 * @author Marjan Fathian
 * @author Ruben Zorrilla
 * @class DeactivateModelPartOperation
 * @ingroup GeoMechanicsApplication
 * @brief This operation emulate the effect of an excavation by deactivating some parts of the computational domain
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) DeactivateModelPartOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DeactivateModelPartOperation
    KRATOS_CLASS_POINTER_DEFINITION(DeactivateModelPartOperation);

    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics.GeoMechanicsApplication", Operation, DeactivateModelPartOperation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", Operation, DeactivateModelPartOperation)

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DeactivateModelPartOperation() : Operation() {}

    /// @brief Model-parameters constructor.
    /// @param rModel Reference to the model container
    /// @param rSettings Input settings
    DeactivateModelPartOperation(
        Model& rModel,
        const Parameters rSettings);

    /// Destructor
    ~DeactivateModelPartOperation() override = default;

    /// Copy constructor
    DeactivateModelPartOperation(const DeactivateModelPartOperation&) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    DeactivateModelPartOperation& operator=(const DeactivateModelPartOperation&) = delete;

    ///@}
    ///@name Operations
    ///@{

    Operation::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const override;

    void Execute() override;

    ///@}
private:
    ///@}
    ///@name Member Variables
    ///@{

    const ModelPart* mpModelPart = nullptr; // Reference to the model part to which the operation is applied

    ///@}
}; // Class DeactivateModelPartOperation

///@} addtogroup block

} // namespace Kratos.