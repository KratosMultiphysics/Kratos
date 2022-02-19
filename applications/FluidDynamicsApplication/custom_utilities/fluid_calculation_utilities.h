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

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

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

        int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
            AssignValue<
                typename std::remove_reference<typename std::tuple_element<0, TRefVariableValuePairArgs>::type>::type,
                typename std::remove_reference<typename std::tuple_element<1, TRefVariableValuePairArgs>::type>::type::Type
                >
                (
                    r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value,
                    std::get<0>(rValueVariablePairs)
                ),
            0)...};

        // this can be removed with fold expressions in c++17
        *dummy = 0;

        for (IndexType c = 1; c < rGeometry.PointsNumber(); ++c) {
            const auto& r_node = rGeometry[c];
            const double shape_function_value = rShapeFunction[c];

            int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
                UpdateValue<
                    typename std::remove_reference<typename std::tuple_element<0, TRefVariableValuePairArgs>::type>::type,
                    typename std::remove_reference<typename std::tuple_element<1, TRefVariableValuePairArgs>::type>::type::Type
                    >
                    (
                        r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value,
                        std::get<0>(rValueVariablePairs)
                    ),
                0)...};

            // this can be removed with fold expressions in c++17
            *dummy = 0;
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

    template<class TOutputDataType, class TInputDataType = TOutputDataType>
    static void AssignValue(
        const TInputDataType& rInput,
        TOutputDataType& rOutput);

    template<class TOutputDataType, class TInputDataType = TOutputDataType>
    static void UpdateValue(
        const TInputDataType& rInput,
        TOutputDataType& rOutput);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_CALCULATION_UTILITIES_H
