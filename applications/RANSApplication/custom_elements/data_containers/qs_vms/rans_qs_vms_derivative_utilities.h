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

#if !defined(KRATOS_RANS_QS_VMS_DERIVATIVE_UTILITIES_H)
#define KRATOS_RANS_QS_VMS_DERIVATIVE_UTILITIES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim>
class RansQSVMSDerivativeUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    using DependentVariablesListType = typename QSVMSDerivativeUtilities<TDim>::DependentVariablesListType;

    ///@}
    ///@name Classes
    ///@{

    template<unsigned int TNumNodes, class TElementData>
    class TurbulenceVariableDerivative : public QSVMSDerivativeUtilities<TDim>::template Derivative<0>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = typename QSVMSDerivativeUtilities<TDim>::template Derivative<0>;

        using ElementDataType = TElementData;

        static constexpr double VelocityDerivativeFactor = 0.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        static constexpr unsigned int TDerivativeDimension = 1;

        ///@}
        ///@name Life Cycle
        ///@{

        TurbulenceVariableDerivative(
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
     * @brief This class is used with k-omega-sst adjoints
     *
     * k-omega-sst turbulence model calculate nu_t using tke, omega, and velocity_gradient.
     * Therefore, in the case of the derivative w.r.t. VELOCITY, we need to calculate
     * velocity_gradient derivatives as well. Then QSVMSDerivatives will apply chain rule
     * to computed velocity_gradient derivatives to convert them to velocity derivatives.
     *
     * @tparam TNumNodes
     */
    template<unsigned int TNumNodes, unsigned int TComponentIndex>
    class KOmegaSSTVelocityDerivative : public QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes, TComponentIndex>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes, TComponentIndex>;

        static constexpr double VelocityDerivativeFactor = 1.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        static constexpr unsigned int TDerivativeDimension = TDim;

        ///@}
        ///@name Life Cycle
        ///@{

        KOmegaSSTVelocityDerivative(
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

        DependentVariablesListType GetEffectiveViscosityDependentVariables() const;

        ///@}
    };

    /**
     * @brief This class is used with k-omega-sst adjoints
     *
     * k-omega-sst turbulence model calculate nu_t using tke, omega, and velocity_gradient.
     * Therefore, in the case of the derivative w.r.t. VELOCITY, we need to calculate
     * velocity_gradient derivatives as well. Then QSVMSDerivatives will apply chain rule
     * to computed velocity_gradient derivatives to convert them to velocity derivatives.
     *
     * @tparam TNumNodes
     */
    template<unsigned int TNumNodes, unsigned int TComponentIndex>
    class KOmegaSSTShapeDerivative : public QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes, TComponentIndex>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = typename QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes, TComponentIndex>;

        static constexpr double VelocityDerivativeFactor = 0.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        static constexpr unsigned int TDerivativeDimension = TDim;

        ///@}
        ///@name Life Cycle
        ///@{

        KOmegaSSTShapeDerivative(
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

        DependentVariablesListType GetEffectiveViscosityDependentVariables() const;

        ///@}
    };

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_RANS_QS_VMS_DERIVATIVE_UTILITIES_H