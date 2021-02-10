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

#if !defined(KRATOS_FLUID_ADJOINT_VARIABLE_INFORMATION_H)
#define KRATOS_FLUID_ADJOINT_VARIABLE_INFORMATION_H

// System includes
#include <array>

// External includes

// Project includes
#include "includes/variables.h"

// Application includes
#include "custom_utilities/fluid_adjoint_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

namespace AdjointVariableInformation
{
///@name Type Definitions
///@{

using IndexType = std::size_t;

template<unsigned int TDim>
using ScalarVariableGradientComponents = std::array<const Variable<double>*, TDim>;

///@}
///@name Kratos classes
///@{

template <unsigned int TDim>
class VariableInformation
{
public:
    ///@name Life Cycle
    ///{

    VariableInformation(
        const Variable<double>& rVariable,
        const ScalarVariableGradientComponents<TDim>& rVariableGradientComponents)
        : mrVariable(rVariable),
          mVariableGradients(rVariableGradientComponents)
    {}

    ///@}
    ///@name Operations
    ///@{

    const Variable<double>& GetVariable() const { return mrVariable; }

    const ScalarVariableGradientComponents<TDim> GetVariableGradientComponents() const { return mVariableGradients; }

    ///@}
private:
    ///@name Private Members
    ///@{

    const Variable<double>& mrVariable;

    const ScalarVariableGradientComponents<TDim> mVariableGradients;

    ///@}
};

template <unsigned int TDim>
class ScalarVariableInformation : public VariableInformation<TDim>
{
public:
    ///@name Life Cycle
    ///@{

    ScalarVariableInformation(
        const Variable<double>& rVariable,
        const ScalarVariableGradientComponents<3>& rVariableGradientComponents)
        : VariableInformation<TDim>(
            rVariable,
            FluidAdjointUtilities<TDim>::GetRelevantGradientVariableComponentList(0, rVariable, rVariableGradientComponents))
    {}

    ///@}
};

template <unsigned int TDim>
class VectorVariableInformation : public VariableInformation<TDim>
{
public:
    ///@name Life Cycle
    ///@{

    VectorVariableInformation(
        const IndexType DirectionIndex,
        const Variable<array_1d<double, 3>>& rVariable,
        const ScalarVariableGradientComponents<3>& rVariableComponents,
        const ScalarVariableGradientComponents<9>& rVariableGradientComponents)
        : VariableInformation<TDim>(
            FluidAdjointUtilities<TDim>::GetRelevantVariable(DirectionIndex, rVariable, rVariableComponents),
            FluidAdjointUtilities<TDim>::GetRelevantGradientVariableComponentList(DirectionIndex, rVariable, rVariableGradientComponents)
        )
    {}

    ///@}
};

template <unsigned int TDim>
class VelocityInformation : public VectorVariableInformation<TDim>
{
public:
    ///@name Life Cycle
    ///@{

    VelocityInformation(const IndexType DirectionIndex)
        : VectorVariableInformation<TDim>(
            DirectionIndex,
            VELOCITY,
            {
                &VELOCITY_X,
                &VELOCITY_Y,
                &VELOCITY_Z
            },
            {
                &VELOCITY_GRADIENT_TENSOR_XX,
                &VELOCITY_GRADIENT_TENSOR_XY,
                &VELOCITY_GRADIENT_TENSOR_XZ,
                &VELOCITY_GRADIENT_TENSOR_YX,
                &VELOCITY_GRADIENT_TENSOR_YY,
                &VELOCITY_GRADIENT_TENSOR_YZ,
                &VELOCITY_GRADIENT_TENSOR_ZX,
                &VELOCITY_GRADIENT_TENSOR_ZY,
                &VELOCITY_GRADIENT_TENSOR_ZZ
            })
    {}

    ///@}
};

template <unsigned int TDim>
class PressureInformation : public VariableInformation<TDim>
{
public:
    ///@name Life Cycle
    ///@{

    PressureInformation()
    : VariableInformation<TDim>(
        PRESSURE,
        {
            &Variable<double>::StaticObject(),
            &Variable<double>::StaticObject()
        })
    {}

    ///@}
};

template <unsigned int TDim>
class ShapeSensitivityInformation : public VectorVariableInformation<TDim>
{
public:
    ///@name Life Cycle
    ///@{

    ShapeSensitivityInformation(const IndexType DirectionIndex)
        : VectorVariableInformation<TDim>(
            DirectionIndex,
            SHAPE_SENSITIVITY,
            {
                &SHAPE_SENSITIVITY_X,
                &SHAPE_SENSITIVITY_Y,
                &SHAPE_SENSITIVITY_Z
            },
            {
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_XX,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_XY,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_XZ,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_YX,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_YY,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_YZ,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_ZX,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_ZY,
                &SHAPE_SENSITIVITY_GRADIENT_TENSOR_ZZ
            })
    {}

    ///@}
};

///@}

} // namespace AdjointVariableInformation

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_VARIABLE_INFORMATION_H