//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:             BSD License
//                            Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_CFD_UTILITIES_H_INCLUDED)
#define KRATOS_CFD_UTILITIES_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace CFDUtilities
{
using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;
using GeometryType = Geometry<NodeType>;

/**
 * @brief Calculate gauss point data for condition
 *
 * The method calculates and retrieves gauss point data namely
 * gauss_weight, and gauss point shape function values at each gauss point.
 * rNContainer matrix rows corresponds to gauss point index
 *
 * @param rGeometry
 * @param rIntegrationMethod
 * @param rGaussWeights
 * @param rNContainer
 */
void CalculateConditionGeometryData(const GeometryType& rGeometry,
                                    const GeometryData::IntegrationMethod& rIntegrationMethod,
                                    Vector& rGaussWeights,
                                    Matrix& rNContainer);

/**
 * @brief Calculates outwards pointing normal for a condition
 *
 * This method calculates outward pointing normal. The geometry should be
 * oriented properly to have outwards pointing normal. Use
 * TetrahedralMeshOrientationCheck before using this method to calculate normal
 * since this process checks and corrects orientation of geometry nodes.
 * Magnitude of normal contains area(in 3D), length (in 2D) of the condition.
 *
 * @tparam TDim       Dimensionality of condition. Either 2 or 3
 * @param rNormal     Calculated normal
 * @param rCondition  Condition
 *
 * @see TetrahedralMeshOrientationCheck
 */

template <unsigned int TDim>
void CalculateConditionNormal(array_1d<double, 3>& rNormal, const ConditionType& rCondition);

/**
 * @brief Calculates wall height for condition
 *
 * Wall height is for given condition is calculated. Wall height is the wall normal distance
 * between center of condition and center of parent element for that particular condition.
 * TetrahedralMeshOrientationCheck needs to be used with ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS flag
 * to assign parent element for each condition to use this method.
 *
 * @param rCondition  Condition for which wall height is calculated
 * @param rNormal     Wall normal (outwards pointing)
 * @return double     Wall height
 *
 * @see TetrahedralMeshOrientationCheck
 */
double CalculateConditionWallHeight(const ConditionType& rCondition,
                                    const array_1d<double, 3>& rNormal);

/**
 * @brief Evaluates variable value at gauss point
 *
 * This method evaluates a variable for a given geometry at a given gauss point.
 * The gauss point is chosen by providing proper shape function values at that
 * gauss point via rShapeFunction.
 *
 * @tparam TDataType        Data type
 * @param rGeometry         Geometry for gauss point evaluation
 * @param rVariable         Variable
 * @param rShapeFunction    Shape function values evaluated at gauss point
 * @param Step              Step
 * @return TDataType        variable value evaluated at gauss point
 */
template <typename TDataType>
TDataType EvaluateInPoint(const GeometryType& rGeometry,
                          const Variable<TDataType>& rVariable,
                          const Vector& rShapeFunction,
                          const int Step = 0)
{
    const unsigned int number_of_nodes = rGeometry.PointsNumber();
    TDataType value =
        rGeometry[0].FastGetSolutionStepValue(rVariable, Step) * rShapeFunction[0];
    for (unsigned int c = 1; c < number_of_nodes; ++c)
    {
        value += rGeometry[c].FastGetSolutionStepValue(rVariable, Step) *
                 rShapeFunction[c];
    }

    return value;
}

/**
 * @brief Calculates neighbour conditions
 *
 * @param rModelPart
 */
void KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateNumberOfNeighbourConditions(ModelPart& rModelPart);

/**
 * @brief Calculates $y^+$ limit for linear-logarithmic wall law
 *
 * Linear law
 * \[
 *      y^+ = \frac{u_\tau y}{\nu}
 * \]
 *
 * Log law
 * \[
 *      u^+ = \frac{1}{\kappa}log\left(y^+\right) + \beta
 * \]
 *
 * Where $\kappa$ is the VonKarman constant, and $\beta$ is the wall smoothness parameter
 * $\nu$ is the kinematic viscosity, $u_\tau$ is the wall friction velocity, $y$ is the wall height.
 *
 * At the limit following relation is assumed.
 *
 * \[
 *      y^+_{limit} =  \frac{1}{\kappa}log\left(y^+_{limit}\right) + \beta
 * \]
 *
 * Therefore, this method tries to find limiting $y^+$ where linear and log law meets.
 *
 * @param VonKarman             Von Karman constant ($\kappa$)
 * @param WallSmoothness        Wall smoothness parameter ($\beta$)
 * @param MaxIterations         Number of maximum iterations to use to identify limit
 * @param Tolerance             Tolerance for convergence
 * @return double               $y^+_{limit}$ limit at the boundary
 */
double KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateLinearLogarithmicWallFunctionBasedYPlusLimit(const double VonKarman = 0.41,
                                                          const double WallSmoothness = 5.2,
                                                          const int MaxIterations = 20,
                                                          const double Tolerance = 1e-6);

/**
 * @brief Calculate $y^+$ and $u_\tau$ based on log and linear laws
 *
 * $y^+$ and $u_\tau$ is calculated using either linear or log wall laws.
 *
 * In the linear region where $y^+ \le y^+_{limit}$
 * \[
 *      y^+ = u^+
 * \]
 * where $y^+=\frac{u_\tau y}{\nu}$ and $u^+=\frac{|\mathbf{u}|}{u_\tau}$. In here $u_\tau$ is the
 * friction velocity. $y$ is the wall height calculated from CalculateConditionWallHeight method.
 * $\nu$ is the kinematic viscosity. $\mathbf{u}$ is the wall velocity, where tangential component
 * is calculated inside this method and used in the wall law.
 *
 * If $y^+>y^+_{limit}$ then log law is applied where
 * \[
 *      u^+ = \frac{1}{\kappa}log\left(y^+\right) + \beta
 * \]
 * Where $y^+$ and $u^+$ definitions remains the same as in linear region. $\kappa$ is the Von Karman
 * constant, and $\beta$ is the wall smoothness parameter. Iterative Newton-Raphson method is used
 * to calculate proper $y^+$ and $u_\tau$.
 *
 * Finally calculated $u_\tau$ is presented in the opposite direction to tangential velocity at wall.
 *
 * It is advised to use this method with slip boundary conditions, where wall laws are applied.
 *
 * @param rFrictionVelocity         Output friction velocity ($u_\tau$)
 * @param rWallVelocity             Wall velocity (tangential component will be derrived within this method)
 * @param rNormal                   Wall outwards pointing normal for condition
 * @param KinematicViscosity        Kinematic viscosity ($\nu$)
 * @param WallHeight                Wall height ($y$)
 * @param VonKarman                 Von Karman constant ($\kappa$)
 * @param WallSmoothness            Wall smoothness parameter ($\beta$)
 * @param MaxIterations             Max iterations used in Newton-Raphson solver
 * @param Tolerance                 Tolerance to be converged for Newton-Raphson solver
 * @return double                   Output $y^+$ value
 */
double KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau(
        array_1d<double, 3>& rFrictionVelocity,
        const array_1d<double, 3>& rWallVelocity,
        const array_1d<double, 3>& rNormal,
        const double KinematicViscosity,
        const double WallHeight,
        const double VonKarman = 0.41,
        const double WallSmoothness = 5.2,
        const int MaxIterations = 20,
        const double Tolerance = 1e-6);

/**
 * @brief Calculates $y^+$ at the cell center (no wall laws assumed)
 *
 * This method calculates $y^+$ at center of the parent element of the condition.
 * No wall laws are assumed.
 *
 * \[
 *      u_\tau = \sqrt{\frac{\tau}{\rho}}
 * \]
 *
 * Where $u_\tau$ is the wall friction velocity, $\tau$ is the shear stress at walls.
 * $\rho$ is the density.
 *
 * Total force acting on condition surface ($\mathbf{R}$) is calculate using,
 * \[
 *      \mathbf{R} = \sum^M_i \frac{\mathbf{R}_i}{N_i}
 * \]
 * Where $M$ is number of nodes in the condition. $N_i$ is number of neighbour conditions at $i^{th}$ node.
 * \mathbf{R}_i is the reaction force at $i^{th}$ node.
 *
 * $\tau$ is calculated using
 * \[
 *      \tau = \frac{\mathbf{R} - \left(\mathbf{R}\cdot\mathbf{n}\right)\mathbf{n}}{A}
 * \]
 * Where $A$ is the surface area of the condition (i.e. line length in 2D). $\mathbf{n}$ is the
 * outward pointing unit normal. $\mathbf{R}$ is the reaction acting upon fluid flow. Finally direction
 * of $u_\tau$ is kept as the tangential direction of $\mathbf{R}$.
 *
 * It is advised to use this method with no slip boundary conditions.
 *
 * @param rFrictionVelocity     Output friction velocity $u_\tau$
 * @param rReaction             Reaction force vector (tangential component is calculated within the method)
 * @param rNormal               Outwards pointing condition normal (Magnitude contains surface area)
 * @param Density               Density of the fluid
 * @param KinematicViscosity    Kinematic viscosity of fluid
 * @param WallHeight            Wall height at center of the parent element
 * @return double               Output $y^+$ at center of the parent element
 */
double KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateReactionBasedYPlusUTau(array_1d<double, 3>& rFrictionVelocity,
                                    const array_1d<double, 3>& rReaction,
                                    const array_1d<double, 3>& rNormal,
                                    const double Density,
                                    const double KinematicViscosity,
                                    const double WallHeight);

/**
 * @brief Calculates $y^+$ and $u_\tau$ for given method
 *
 * This holds common steps required to calculate and store $y^+$ and $u_\tau$ for given model part's conditions.
 *
 * @param rModelPart                        Model part
 * @param rKinematicViscosityVariable       Variable for kinematic viscosity
 * @param rYPlusAndUTauCalculationMethod    $y^+$ and $u_\tau$ calculation method.
 *
 * @see CalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau
 * @see CalculateReactionBasedYPlusUTau
 */
void CalculateYPlusAndUTauForConditions(
    ModelPart& rModelPart,
    const Variable<double>& rKinematicViscosityVariable,
    const std::function<double(
        array_1d<double, 3>&, const GeometryType&, const array_1d<double, 3>&, const Vector&, const double, const double, const double)>&
        rYPlusAndUTauCalculationMethod);
/**
 * @brief Calculates $y^+$ and $u_\tau$ based on reaction for conditions
 *
 * This method calculates $y^+$ and $u_\tau$ for conditions in given model part based on reaction.
 * Calculated $y^+$ and $u_\tau$ is stored at Y_PLUS, and FRICTION_VELOCITY variables in condition
 * data value container.
 *
 * It is advised to use model parts where only no slip boundary condition is used for conditions
 *
 * @param rModelPart                        Model part
 * @param rKinematicViscosityVariable       Kinematic viscosity variable
 * @param rReactionVariable                 Reaction variable
 *
 * @see CalculateReactionBasedYPlusUTau
 * @see CalculateYPlusAndUTauForConditions
 */
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateYPlusAndUTauForConditionsBasedOnReaction(
    ModelPart& rModelPart,
    const Variable<double>& rKinematicViscosityVariable,
    const Variable<array_1d<double, 3>>& rReactionVariable);

/**
 * @brief Calculates $y^+$ and $u_\tau$ based on linear-log wall laws for conditions
 *
 * This method calculates $y^+$ and $u_\tau$ for conditions in given model part based on
 * linear and log wall laws. Calculated $y^+$ and $u_\tau$ is stored at Y_PLUS, and FRICTION_VELOCITY variables in condition
 * data value container.
 *
 * It is advied to use model parts where only slip boundary with wall functions are used for conditions
 *
 * @param rModelPart                    Model part where $y^+$ and $u_\tau$ is calculated for conditions
 * @param rKinematicViscosityVariable   Kinematic viscosity variable
 * @param VonKarman                     Von Karman constant
 * @param WallSmoothness                Wall smoothness constant
 * @param MaxIterations                 Max iterations used in Newton-Raphson solver
 * @param Tolerance                     Tolerance used in Newton-Raphson solver
 *
 * @see CalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau
 * @see CalculateYPlusAndUTauForConditions
 */
void KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
        ModelPart& rModelPart,
        const Variable<double>& rKinematicViscosityVariable,
        const double VonKarman = 0.41,
        const double WallSmoothness = 5.2,
        const int MaxIterations = 20,
        const double Tolerance = 1e-6);

/**
 * @brief Distributes variable values in conditions to nodes
 *
 * This method distributes variables values stored in condition data value container to nodes.
 * Constant weighting is used. The distributed variable value is stored in nodal non-historical
 * data value container under the same variable.
 *
 * @tparam TDataType    Data type
 * @param rModelPart    Model part
 * @param rVariable     Variable to be distributed
 */
template <typename TDataType>
void DistributeConditionVariableToNodes(ModelPart& rModelPart,
                                        const Variable<TDataType>& rVariable);

} // namespace CFDUtilities

} // namespace Kratos.

#endif // KRATOS_CFD_UTILITIES_H_INCLUDED  defined
