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

#pragma once

// System includes
#include <array>

// External includes

// Project includes
#include "containers/variable.h"
#include "fluid_dynamics_application_variables.h"
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_utilities/fluid_adjoint_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief These are helper classes to define derivatives of coefficients of Navier-Stokes.
 *
 * This class defines CalculationContainerHelper classes which can be used to compute
 * analytical Navier-Stokes coefficient derivatives.
 *
 * Coefficients being:
 *      EffectiveVelocity
 *      ElementLength
 *      StrainRate
 *
 * @tparam TDim
 */
template <unsigned int TDim>
class QSVMSDerivativeUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    using DependentVariablesListType = std::vector<
                                            std::tuple<
                                                const Variable<double>&,
                                                std::vector<const Variable<double>*>
                                            >
                                        >;

    using DerivativeGradientsArray = std::array<const Variable<double>*, 9>;

    constexpr static IndexType TStrainSize = (TDim - 1) * 3; // 3 in 2D, 6 in 3D

    ///@}
    ///@name Static Operations
    ///@{

    // Ask Ruben: This is redundant code, can we agree to
    // create qvms_utilities.h, and store this there to be used in primal and adjoint analysis?
    static void CalculateStrainRate(
        Vector& rOutput,
        const Matrix& rNodalVelocity,
        const Matrix& rdNdX);

    static const std::array<const Variable<double>*, TStrainSize> GetStrainRateVariables();

    ///@}
    ///@name Classes
    ///@{

    /**
     * @brief Base class for Derivatives
     *
     * This implements some checks which are required for all the derivatives
     *
     * @tparam TComponentIndex
     */
    template<unsigned int TComponentIndex = 0>
    class Derivative
    {
    public:
        ///@name Type definitions
        ///@{

        static constexpr IndexType ComponentIndex = TComponentIndex;

        ///@}
        ///@name Life Cycle
        ///@{

        Derivative(
            const IndexType NodeIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        ///@}
        ///@name Operations
        ///@{

        DependentVariablesListType GetEffectiveViscosityDependentVariables() const { return DependentVariablesListType({}); }

        ///@}

    protected:
        ///@name Private Members
        ///@{

        const IndexType mNodeIndex;
        const GeometryType& mrGeometry;
        const double mW;
        const Vector& mrN;
        const Matrix& mrdNdX;
        const double mWDerivative;
        const double mDetJDerivative;
        const Matrix& mrdNdXDerivative;

        ///@}
    };

    /**
     * @brief Velocity derivative computation container
     *
     * This class is used to compute derivatives of effective velocity, length
     * and strain rate with respect velocity TComponentIndex component.
     *
     * @tparam TNumNodes        Number of nodes in the element.
     * @tparam TComponentIndex  Component index of the velocity.
     */
    template<unsigned int TNumNodes, unsigned int TComponentIndex>
    class VelocityDerivative : public Derivative<TComponentIndex>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = Derivative<TComponentIndex>;

        static constexpr IndexType ComponentIndex = BaseType::ComponentIndex;

        static constexpr double VelocityDerivativeFactor = 1.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        static constexpr unsigned int TDerivativeDimension = TDim;

        ///@}
        ///@name Life Cycle
        ///@{

        VelocityDerivative(
            const IndexType NodeIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
            : BaseType(NodeIndex, rGeometry, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative)
        {
        }

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        array_1d<double, TDim> CalculateEffectiveVelocityDerivative(const array_1d<double, TDim>& rVelocity) const;

        double CalculateElementLengthDerivative(const double ElementLength) const;

        void CalculateStrainRateDerivative(
            Vector& rOutput,
            const Matrix& rNodalVelocity) const;

        ///@}
    };

    /**
     * @brief Pressure derivative computation container
     *
     * This class is used to compute derivatives of effective velocity, length
     * and strain rate with respect pressure.
     *
     * @tparam TNumNodes        Number of nodes in the element.
     */
    template<unsigned int TNumNodes>
    class PressureDerivative : public Derivative<0>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = Derivative<0>;

        static constexpr IndexType ComponentIndex = BaseType::ComponentIndex;

        static constexpr double VelocityDerivativeFactor = 0.0;

        static constexpr double PressureDerivativeFactor = 1.0;

        static constexpr unsigned int TDerivativeDimension = 1;

        ///@}
        ///@name Life Cycle
        ///@{

        PressureDerivative(
            const IndexType NodeIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
            : BaseType(NodeIndex, rGeometry, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative)
        {
        }

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const { return PRESSURE; }

        array_1d<double, TDim> CalculateEffectiveVelocityDerivative(const array_1d<double, TDim>& rVelocity) const;

        double CalculateElementLengthDerivative(const double ElementLength) const;

        void CalculateStrainRateDerivative(
            Vector& rOutput,
            const Matrix& rNodalVelocity) const;

        ///@}
    };

    /**
     * @brief Shape derivative computation container
     *
     * This class is used to compute derivatives of effective velocity, length
     * and strain rate with respect shape TComponentIndex component.
     *
     * @tparam TNumNodes        Number of nodes in the element.
     * @tparam TComponentIndex  Component index of the shape variable (nodal coordinates component).
     */
    template<unsigned int TNumNodes, unsigned int TComponentIndex>
    class ShapeDerivative : public Derivative<TComponentIndex>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = Derivative<TComponentIndex>;

        static constexpr IndexType ComponentIndex = BaseType::ComponentIndex;

        static constexpr double VelocityDerivativeFactor = 0.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        static constexpr unsigned int TDerivativeDimension = TDim;

        ///@}
        ///@name Life Cycle
        ///@{

        ShapeDerivative(
            const IndexType NodeIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
            : BaseType(NodeIndex, rGeometry, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative)
        {
        }

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        array_1d<double, TDim> CalculateEffectiveVelocityDerivative(const array_1d<double, TDim>& rVelocity) const;

        double CalculateElementLengthDerivative(const double ElementLength) const;

        void CalculateStrainRateDerivative(
            Vector& rOutput,
            const Matrix& rNodalVelocity) const;

        ///@}
    };

    ///@}
private:
    ///@name Private Static Operations
    ///@{

    static void CalculateStrainRateVelocityDerivative(
        Vector& rOutput,
        const IndexType DerivativeNodeIndex,
        const IndexType DerivativeDirectionIndex,
        const Matrix& rdNdX);

    ///@}
};

///@}

} // namespace Kratos