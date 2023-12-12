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
 * @class ApplyExcavationOperation
 * @ingroup GeoMechanicsApplication
 * @brief This operation emulate the effect of an excavation by deactivating some parts of the computational domain
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyExcavationOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyExcavationOperation
    KRATOS_CLASS_POINTER_DEFINITION(ApplyExcavationOperation);

    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics.GeoMechanicsApplication", Operation, ApplyExcavationOperation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", Operation, ApplyExcavationOperation)

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ApplyExcavationOperation() : Operation() {}

    /// @brief Model-parameters constructor.
    /// @param rModel Reference to the model container
    /// @param rSettings Input settings
    ApplyExcavationOperation(
        Model& rModel,
        const Parameters rSettings);

    /// Destructor
    ~ApplyExcavationOperation() override = default;

    /// Copy constructor
    ApplyExcavationOperation(const ApplyExcavationOperation&) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ApplyExcavationOperation& operator=(const ApplyExcavationOperation&) = delete;

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

    const bool mDeactivateSoilPart = true; // Auxiliary deactivation flag
    const ModelPart* mpModelPart = nullptr; // Reference to the model part to which the operation is applied

    ///@}
}; // Class ApplyExcavationOperation

///@} addtogroup block

} // namespace Kratos.