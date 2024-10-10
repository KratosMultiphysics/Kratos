//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLUID_CALCULATION_UTILITIES_H)
#define KRATOS_FLUID_CALCULATION_UTILITIES_H

// system includes
#include <tuple>
#include <type_traits>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidCalculationUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    template<class TDataType>
    using ComponentDataType = std::tuple<TDataType&, const TDataType&>;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Evaluates given list of variable pairs at gauss point locations at step
     *
     * Example:
     *      double density;
     *      array_1d<double, 3> velocity
     *      EvaluateInPoint(rGeometry, rShapeFunction, Step,
     *                      std::tie(density, DENSITY),
     *                      std::tie(velocity, VELOCITY));
     *
     *      The above evaluation will fill density, and velocity variables with gauss point
     *      evaluated DENSITY and VELOCITY at gauss point with shape function values rShapeFunction
     *      in the geometry named rGeometry.
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunction            Shape function values at gauss point
     * @param[in] Step                      Step to be evaluated
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType, const Variable<TDataType>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void EvaluateInPoint(
        const GeometryType& rGeometry,
        const Vector& rShapeFunction,
        const int Step,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        KRATOS_TRY

        const auto& r_node = rGeometry[0];
        const double shape_function_value = rShapeFunction[0];

        (std::apply([&](auto&&... args) {((std::get<0>(args) = std::get<1>(args) * shape_function_value), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step))), ...);

        for (IndexType c = 1; c < rGeometry.PointsNumber(); ++c) {
            const auto& r_node = rGeometry[c];
            const double shape_function_value = rShapeFunction[c];

            (std::apply([&](auto&&... args) {((std::get<0>(args) += std::get<1>(args) * shape_function_value), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step))), ...);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Evaluates given list of variable pairs at gauss point locations at current step
     *
     * Example:
     *      double density;
     *      array_1d<double, 3> velocity
     *      EvaluateInPoint(rGeometry, rShapeFunction,
     *                      std::tie(density, DENSITY),
     *                      std::tie(velocity, VELOCITY));
     *
     *      The above evaluation will fill density, and velocity variables with gauss point
     *      evaluated DENSITY and VELOCITY at gauss point with shape function values rShapeFunction
     *      in the geometry named rGeometry.
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunction            Shape function values at gauss point
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType, const Variable<TDataType>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void inline EvaluateInPoint(
        const GeometryType& rGeometry,
        const Vector& rShapeFunction,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        EvaluateInPoint<TRefVariableValuePairArgs...>(
            rGeometry, rShapeFunction, 0, rValueVariablePairs...);
    }

    /**
     * @brief Evaluates given list of non historical variable pairs at gauss point locations
     *
     * Example:
     *      double density;
     *      array_1d<double, 3> velocity
     *      EvaluateInPoint(rGeometry, rShapeFunction,
     *                      std::tie(density, DENSITY),
     *                      std::tie(velocity, VELOCITY));
     *
     *      The above evaluation will fill density, and velocity variables with gauss point
     *      evaluated DENSITY and VELOCITY at gauss point with shape function values rShapeFunction
     *      in the geometry named rGeometry.
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunction            Shape function values at gauss point
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType, const Variable<TDataType>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void EvaluateNonHistoricalInPoint(
        const GeometryType& rGeometry,
        const Vector& rShapeFunction,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        KRATOS_TRY

        const auto& r_node = rGeometry[0];
        const double shape_function_value = rShapeFunction[0];

        (std::apply([&](auto&&... args) {((std::get<0>(args) = std::get<1>(args) * shape_function_value), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.GetValue(std::get<1>(rValueVariablePairs)))), ...);

        for (IndexType c = 1; c < rGeometry.PointsNumber(); ++c) {
            const auto& r_node = rGeometry[c];
            const double shape_function_value = rShapeFunction[c];

            (std::apply([&](auto&&... args) {((std::get<0>(args) += std::get<1>(args) * shape_function_value), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.GetValue(std::get<1>(rValueVariablePairs)))), ...);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Evaluates gradients of given list of variable pairs at gauss point locations at step
     *
     * Example:
     *      array_1d<double, 3> density_gradient;
     *      BoundedMatrix<double, 3, 3> velocity_gradient
     *      EvaluateGradientInPoint(rGeometry, rShapeFunctionDerivatives, Step,
     *                              std::tie(density_gradient, DENSITY),
     *                              std::tie(velocity_gradient, VELOCITY));
     *
     *      The above evaluation will fill density_gradient, and velocity_gradient variables with gauss point
     *      evaluated gradient of DENSITY and VELOCITY at gauss point with shape function derivative values rShapeFunctionDerivatives
     *      in the geometry named rGeometry.
     *
     *      The output gradient vectors will be having variable direction components in rows in the case where
     *      variable is of type scalar
     *
     *      The output gradient matrices will be having variable components in rows, direction components in columns in the case where
     *      variable is of type vector
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunctionDerivatives Shape function  derivative values at gauss point
     * @param[in] Step                      Step to be evaluated
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType1, const Variable<TDataType2>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void EvaluateGradientInPoint(
        const GeometryType& rGeometry,
        const Matrix& rShapeFunctionDerivatives,
        const int Step,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        KRATOS_TRY

        const auto& r_node = rGeometry[0];
        const Vector& shape_function_derivative = row(rShapeFunctionDerivatives, 0);

        for (IndexType i = 0; i < rShapeFunctionDerivatives.size2(); ++i) {
            (std::apply([&](auto&&... args) {((std::get<0>(args) = std::get<1>(args) * shape_function_derivative[i]), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step), i)), ...);
        }

        for (IndexType c = 1; c < rGeometry.PointsNumber(); ++c) {
            const auto& r_node = rGeometry[c];
            const Vector& shape_function_derivative = row(rShapeFunctionDerivatives, c);

            for (IndexType i = 0; i < rShapeFunctionDerivatives.size2(); ++i) {
                (std::apply([&](auto&&... args) {((std::get<0>(args) += std::get<1>(args) * shape_function_derivative[i]), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step), i)), ...);
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Evaluates gradients of given list of variable pairs at gauss point locations at current step
     *
     * Example:
     *      array_1d<double, 3> density_gradient;
     *      BoundedMatrix<double, 3, 3> velocity_gradient
     *      EvaluateGradientInPoint(rGeometry, rShapeFunctionDerivatives, Step,
     *                              std::tie(density_gradient, DENSITY),
     *                              std::tie(velocity_gradient, VELOCITY));
     *
     *      The above evaluation will fill density_gradient, and velocity_gradient variables with gauss point
     *      evaluated gradient of DENSITY and VELOCITY at gauss point with shape function derivative values rShapeFunctionDerivatives
     *      in the geometry named rGeometry.
     *
     *      The output gradient vectors will be having variable direction components in rows in the case where
     *      variable is of type scalar
     *
     *      The output gradient matrices will be having variable components in rows, direction components in columns in the case where
     *      variable is of type vector
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunctionDerivatives Shape function  derivative values at gauss point
     * @param[in] Step                      Step to be evaluated
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType1, const Variable<TDataType2>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void inline EvaluateGradientInPoint(
        const GeometryType& rGeometry,
        const Matrix& rShapeFunctionDerivatives,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        EvaluateGradientInPoint<TRefVariableValuePairArgs...>(
            rGeometry, rShapeFunctionDerivatives, 0, rValueVariablePairs...);
    }

    /**
     * @brief Evaluates non historical gradients of given list of variable pairs at gauss point locations
     *
     * Example:
     *      array_1d<double, 3> density_gradient;
     *      BoundedMatrix<double, 3, 3> velocity_gradient
     *      EvaluateGradientInPoint(rGeometry, rShapeFunctionDerivatives, Step,
     *                              std::tie(density_gradient, DENSITY),
     *                              std::tie(velocity_gradient, VELOCITY));
     *
     *      The above evaluation will fill density_gradient, and velocity_gradient variables with gauss point
     *      evaluated gradient of DENSITY and VELOCITY at gauss point with shape function derivative values rShapeFunctionDerivatives
     *      in the geometry named rGeometry.
     *
     *      The output gradient vectors will be having variable direction components in rows in the case where
     *      variable is of type scalar
     *
     *      The output gradient matrices will be having variable components in rows, direction components in columns in the case where
     *      variable is of type vector
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunctionDerivatives Shape function  derivative values at gauss point
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType1, const Variable<TDataType2>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void EvaluateNonHistoricalGradientInPoint(
        const GeometryType& rGeometry,
        const Matrix& rShapeFunctionDerivatives,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        KRATOS_TRY

        const auto& r_node = rGeometry[0];
        const Vector& shape_function_derivative = row(rShapeFunctionDerivatives, 0);

        for (IndexType i = 0; i < rShapeFunctionDerivatives.size2(); ++i) {
            (std::apply([&](auto&&... args) {((std::get<0>(args) = std::get<1>(args) * shape_function_derivative[i]), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.GetValue(std::get<1>(rValueVariablePairs)), i)), ...);
        }

        for (IndexType c = 1; c < rGeometry.PointsNumber(); ++c) {
            const auto& r_node = rGeometry[c];
            const Vector& shape_function_derivative = row(rShapeFunctionDerivatives, c);

            for (IndexType i = 0; i < rShapeFunctionDerivatives.size2(); ++i) {
                (std::apply([&](auto&&... args) {((std::get<0>(args) += std::get<1>(args) * shape_function_derivative[i]), ...);}, GetComponentsArray(std::get<0>(rValueVariablePairs), r_node.GetValue(std::get<1>(rValueVariablePairs)), i)), ...);
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Get a sub vector from a vector
     *
     * This method returns a sub vector from a vector with size
     * TSize, starting from Position index.
     *
     * @tparam TSize            Size of the output vector
     * @param rOutput           Output vector
     * @param rInput            Input vector
     * @param Position          Starting index in input vector
     */
    template <unsigned int TSize>
    static void ReadSubVector(
        BoundedVector<double, TSize>& rOutput,
        const Vector& rInput,
        const IndexType Position)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rInput.size() < Position + TSize)
            << "rInput size is not enough for sub vector retireval. [ "
               "rInput.size() = "
            << rInput.size() << ", Requested Position = " << Position
            << ", Requested vector size = " << TSize << " ].\n";

        for (IndexType i = 0; i < TSize; ++i) {
            rOutput[i] = rInput[Position + i];
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Adds values of a vector to a given matrix
     *
     * This method adds values of a vector to a matrix where
     * matrix row is given by RowIndex, and matrix column starting
     * index is given by ColumnIndex. The vector values are added along
     * the row of the given matrix.
     *
     * @tparam TSize                Input vector size
     * @param rOutput               Output matrix
     * @param rInput                Input vector
     * @param RowIndex              Row index
     * @param ColumnIndex           Column index
     */
    template<unsigned int TSize>
    static void AddSubVector(
        Matrix& rOutput,
        const BoundedVector<double, TSize>& rInput,
        const unsigned int RowIndex,
        const unsigned int ColumnIndex)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(RowIndex >= rOutput.size1())
            << "RowIndex is larger than or equal to rOutput number of rows. [ "
               "rOutput.size1() = "
            << rOutput.size1() << ", RowIndex = " << RowIndex << " ].\n";

        KRATOS_DEBUG_ERROR_IF(ColumnIndex + TSize > rOutput.size2())
            << "rOutput matrix does not have sufficient number of columns [ "
               "rOutput.size2() = "
            << rOutput.size2() << ", ColumnIndex = " << ColumnIndex
            << ", TSize = " << TSize << " ].\n";

        for (unsigned int i = 0; i < TSize; ++i) {
            rOutput(RowIndex, ColumnIndex + i) += rInput[i];
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates logarithmic and linear wall law region limit
     *
     * This method calculates logarithmic and linear wall law region limit according to following
     * formulare
     *
     * \[
     *      y^+ = \frac{1}{\kappa}ln\left(y^+\right) + \beta
     * \]
     *
     * @param Kappa             $\kappa$ value
     * @param Beta              $\beta$ value
     * @param MaxIterations     Number of max iterations to consider
     * @param Tolerance         Tolerance for the iterative solve
     * @return double           The limiting $y^+$ value
     */
    static double CalculateLogarithmicYPlusLimit(
        const double Kappa,
        const double Beta,
        const int MaxIterations = 20,
        const double Tolerance = 1e-6);

    /**
     * @brief Calculates linear or logarithmic y_plus value
     *
     * This method calculates $y^+$ values based on the following formulare.
     *
     * if $y^+ \geq y^+_{limit}$ using logarithmic wall law
     *
     * \[
     *      u^+ = \frac{1}{\kappa}ln\left(y^+\right) + \beta
     * \]
     *
     * else using linear wall law
     *
     * \[
     *      u^+ = y^+
     * \]
     *
     * Where
     * \[
     *      u^+ = \frac{|\underline{u}|}{u_\tau}
     * \]
     * \[
     *      y^+ = \frac{u_\tau y}{\nu}
     * \]
     *
     * $u_\tau$ is the friction velocity, $y$ is the wall distance, and $\nu$ is the fluid kinematic viscosity.
     *
     * @param WallVelocityMagnitude         $|\underline{u}|$
     * @param WallHeight                    $y$
     * @param KinematicViscosity            $\nu$
     * @param Kappa                         $\kappa$
     * @param Beta                          $\beta$
     * @param YPlusLimit                    $y^+_{limit}$
     * @param MaxIterations                 Number of max iterations to solve for $y^+$
     * @param Tolerance                     Tolerance to be used for solving $y^+$
     * @return double                       The $y^+$
     */
    static double CalculateLogarithmicYPlus(
        const double WallVelocityMagnitude,
        const double WallHeight,
        const double KinematicViscosity,
        const double Kappa,
        const double Beta,
        const double YPlusLimit,
        const int MaxIterations = 20,
        const double Tolerance = 1e-6);

    ///@}

private:
    ///@name Private Operations
    ///@{

    // gauss point evaluation templates
    static auto inline GetComponentsArray(double& rOutput, const double& rInput) { return std::array<ComponentDataType<double>, 1>{std::tie(rOutput, rInput)};}
    static auto inline GetComponentsArray(array_1d<double, 2>& rOutput, const array_1d<double, 3>& rInput) { return std::array<ComponentDataType<double>, 2>{std::tie(rOutput[0], rInput[0]), std::tie(rOutput[1], rInput[1])};}
    static auto inline GetComponentsArray(array_1d<double, 3>& rOutput, const array_1d<double, 3>& rInput) { return std::array<ComponentDataType<double>, 3>{std::tie(rOutput[0], rInput[0]), std::tie(rOutput[1], rInput[1]), std::tie(rOutput[2], rInput[2])};}
    static auto inline GetComponentsArray(Vector& rOutput, const Vector& rInput) { return std::array<ComponentDataType<Vector>, 1>{std::tie(rOutput, rInput)};}
    static auto inline GetComponentsArray(Matrix& rOutput, const Matrix& rInput) { return std::array<ComponentDataType<Matrix>, 1>{std::tie(rOutput, rInput)};}

    // gauss point gradient evaluation templates
    static auto inline GetComponentsArray(array_1d<double, 2>& rOutput, const double& rInput, const IndexType DerivativeIndex) { return std::array<ComponentDataType<double>, 1>{std::tie(rOutput[DerivativeIndex], rInput)};}
    static auto inline GetComponentsArray(array_1d<double, 3>& rOutput, const double& rInput, const IndexType DerivativeIndex) { return std::array<ComponentDataType<double>, 1>{std::tie(rOutput[DerivativeIndex], rInput)};}
    static auto inline GetComponentsArray(BoundedMatrix<double, 2, 2>& rOutput, const array_1d<double, 3>& rInput, const IndexType DerivativeIndex) { return std::array<ComponentDataType<double>, 2>{std::tie(rOutput(0, DerivativeIndex), rInput[0]), std::tie(rOutput(1, DerivativeIndex), rInput[1])};}
    static auto inline GetComponentsArray(BoundedMatrix<double, 3, 3>& rOutput, const array_1d<double, 3>& rInput, const IndexType DerivativeIndex) { return std::array<ComponentDataType<double>, 3>{std::tie(rOutput(0, DerivativeIndex), rInput[0]), std::tie(rOutput(1, DerivativeIndex), rInput[1]), std::tie(rOutput(2, DerivativeIndex), rInput[2])};}

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_CALCULATION_UTILITIES_H
