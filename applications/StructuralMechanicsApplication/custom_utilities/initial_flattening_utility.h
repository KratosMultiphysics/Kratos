// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class InitialFlatteningUtility
 * @brief Computes an initial guess for a cutting-pattern flattening problem by projecting the
 * (reference) nodal positions of a membrane model part onto a plane
 * @details
 * @note 
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InitialFlatteningUtility
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(InitialFlatteningUtility);

    typedef array_1d<double, 3> Vector3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor, not accessible
    InitialFlatteningUtility() = delete;

    /// Destructor.
    virtual ~InitialFlatteningUtility() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the initial flattening projection and sets it as the initial DISPLACEMENT.
     */
    static void Execute(ModelPart& rModelPart, Parameters ThisParameters);

    ///@}

private:
    ///@name Private Operations
    ///@{

    static Vector3 GetFixedDirection(Parameters ThisParameters);

    static Vector3 ComputeMeanSurfaceNormal(const ModelPart& rModelPart);

    static void ProjectNodesOntoPlane(ModelPart& rModelPart, const Vector3& rNormal, const int EchoLevel);

    ///@}

}; // Class InitialFlatteningUtility

///@}

///@} addtogroup block

}  // namespace Kratos.
