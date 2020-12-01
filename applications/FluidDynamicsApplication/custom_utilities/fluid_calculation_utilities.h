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
