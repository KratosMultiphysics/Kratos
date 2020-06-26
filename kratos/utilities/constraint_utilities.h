//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Klaus B Sautter
//

#if !defined(KRATOS_CONSTRAINT_UTILITIES)
#define KRATOS_CONSTRAINT_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @namespace ConstraintUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for the computation of the MPC
 * @author Vicente Mataix Ferrandiz
 */
namespace ConstraintUtilities
{
    /**
     * @brief This method computes the active dofs
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rActiveDofs Vector containing the active dofs
     * @param rDofSet The whole set of dofs
     */
    void KRATOS_API(KRATOS_CORE) ComputeActiveDofs(
        ModelPart& rModelPart,
        std::vector<int>& rActiveDofs,
        const ModelPart::DofsArrayType& rDofSet
        );

    /**
     * @brief This method resets the values of the slave dofs
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) ResetSlaveDofs(ModelPart& rModelPart);

    /**
     * @brief This method resets the values of the slave dofs
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) ApplyConstraints(ModelPart& rModelPart);

    /**
     * @brief This method precomputes the contribution of the explicit MPC over nodal residual forces
     * @param rModelPart The model of the problem to solve
     * @param rDofVariableNames The name of the Dof variables to check
     * @param rResidualDofVariableNames The name name of the corresponding residual variables
     */
    void KRATOS_API(KRATOS_CORE) PreComputeExplicitConstraintConstribution(
        ModelPart& rModelPart,
        const std::vector<std::string>& rDofVariableNames,
        const std::vector<std::string>& rResidualDofVariableNames
        );

    /**
     * @brief This method precomputes the contribution of the explicit MPC over nodal masses and inertias
     * @todo The inertia must be computed using the Steiner theorem https://en.wikipedia.org/wiki/Parallel_axis_theorem
     * @param rModelPart The model of the problem to solve
     * @param DofDisplacementVariableName The name of the Dof variables to check
     * @param MassVariableName The name of the variable of the nodal mass
     * @param DofRotationVariableName The name of the rotational variable name
     * @param InertiaVariableName The inertia variable to be considered
     */
    void KRATOS_API(KRATOS_CORE) PreComputeExplicitConstraintMassAndInertia(
        ModelPart& rModelPart,
        const std::string DofDisplacementVariableName = "DISPLACEMENT",
        const std::string MassVariableName = "NODAL_MASS",
        const std::string DofRotationVariableName = "ROTATION",
        const std::string InertiaVariableName = "NODAL_INERTIA_TENSOR"
        );

}; // namespace ConstraintUtilities
}  // namespace Kratos
#endif /* KRATOS_CONSTRAINT_UTILITIES defined */
