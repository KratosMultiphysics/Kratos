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

#if !defined(KRATOS_QS_VMS_DERIVATIVE_UTILITIES_H)
#define KRATOS_QS_VMS_DERIVATIVE_UTILITIES_H

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

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim>
class QSVMSDerivativeUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    ///@}
    ///@name Static Operations
    ///@{

    // Ask Ruben: This is redundant code, can we agree to
    // create qvms_utilities.h, and store this there to be used in primal and adjoint analysis?
    static void CalculateStrainRate(
        Vector& rOutput,
        const Matrix& rNodalVelocity,
        const Matrix& rdNdX);

    static void CalculateGradient(
        array_1d<double, 3>& rOutput,
        const Variable<double>& rVariable,
        const GeometryType& rGeometry,
        const Matrix& rdNdX);

    static void CalculateGradient(
        BoundedMatrix<double, TDim, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rVariable,
        const GeometryType& rGeometry,
        const Matrix& rdNdX);

    ///@}
    ///@name Classes
    ///@{

    class Derivative
    {
    public:
        ///@name Life Cycle
        ///@{

        Derivative(
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        virtual ~Derivative() = default;

        ///@}

    protected:
        ///@name Private Members
        ///@{

        const IndexType mNodeIndex;
        const IndexType mDirectionIndex;
        const GeometryType& mrGeometry;
        const double mW;
        const Vector& mrN;
        const Matrix& mrdNdX;
        const double mWDerivative;
        const double mDetJDerivative;
        const Matrix& mrdNdXDerivative;

        ///@}
    };

    class VelocityDerivative : public Derivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = Derivative;

        static constexpr double VelocityDerivativeFactor = 1.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        VelocityDerivative(
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
            : BaseType(NodeIndex, DirectionIndex, rGeometry, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative)
        {
        }

        ~VelocityDerivative() override = default;

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        array_1d<double, 3> CalculateEffectiveVelocityDerivative(const array_1d<double, 3>& rVelocity) const;

        double CalculateElementLengthDerivative(const double ElementLength) const;

        void CalculateStrainRateDerivative(
            Vector& rOutput,
            const Matrix& rNodalVelocity) const;

        ///@}
    };

    class PressureDerivative : public Derivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = Derivative;

        static constexpr double VelocityDerivativeFactor = 0.0;

        static constexpr double PressureDerivativeFactor = 1.0;

        ///@}
        ///@name Life Cycle
        ///@{

        PressureDerivative(
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
            : BaseType(NodeIndex, DirectionIndex, rGeometry, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative)
        {
        }

        ~PressureDerivative() override = default;

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        array_1d<double, 3> CalculateEffectiveVelocityDerivative(const array_1d<double, 3>& rVelocity) const;

        double CalculateElementLengthDerivative(const double ElementLength) const;

        void CalculateStrainRateDerivative(
            Vector& rOutput,
            const Matrix& rNodalVelocity) const;

        ///@}
    };

    class ShapeDerivative : public Derivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = Derivative;

        static constexpr double VelocityDerivativeFactor = 0.0;

        static constexpr double PressureDerivativeFactor = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        ShapeDerivative(
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const GeometryType& rGeometry,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
            : BaseType(NodeIndex, DirectionIndex, rGeometry, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative)
        {
        }

        ~ShapeDerivative() override = default;

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        array_1d<double, 3> CalculateEffectiveVelocityDerivative(const array_1d<double, 3>& rVelocity) const;

        double CalculateElementLengthDerivative(const double ElementLength) const;

        void CalculateStrainRateDerivative(
            Vector& rOutput,
            const Matrix& rNodalVelocity) const;

        ///@}
    };

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_QS_VMS_DERIVATIVE_UTILITIES_H