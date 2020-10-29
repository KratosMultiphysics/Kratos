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

#if !defined(KRATOS_K_EPSILON_HIGH_RE_DERIVATIVE_ELEMENT_DATA_EPSILON_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_K_EPSILON_HIGH_RE_DERIVATIVE_ELEMENT_DATA_EPSILON_ELEMENT_DATA_H_INCLUDED

// System includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_adjoint_element_data.h"
#include "custom_elements/convection_diffusion_reaction_sensitivity_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"

namespace Kratos
{
///@name  Functions
///@{

namespace KEpsilonAdjointElementData
{
template <unsigned int TDim, unsigned int TNumNodes>
class EpsilonAdjointElementData
{
public:
    ///@name Public forward declarations
    ///@{

    class AdjointElementData;

    class SensitivityElementData;

    ///@}
    ///@name Public classes
    ///@{

    class VelocityDerivatives
        : public ConvectionDiffusionReactionAdjointElementData<TDim, TNumNodes, AdjointElementData>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = ConvectionDiffusionReactionAdjointElementData<TDim, TNumNodes, AdjointElementData>;
        static constexpr unsigned int TDerivativesSize = BaseType::TDerivativesSize;

        ///@}
        ///@name Life cycle
        ///@{

        VelocityDerivatives(const AdjointElementData& rElementData);

        ///@}
        ///@name Public operations
        ///@{

        static const Variable<array_1d<double, 3>>& GetDerivativeVariable();

        void CalculateEffectiveVelocityDerivatives(
            BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateEffectiveKinematicViscosityDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateReactionTermDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateSourceTermDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        ///@}
    };

    class KDerivatives
        : public ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, AdjointElementData>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, AdjointElementData>;
        static constexpr unsigned int TDerivativesSize = BaseType::TDerivativesSize;

        ///@}
        ///@name Life cycle
        ///@{

        KDerivatives(const AdjointElementData& rElementData);

        ///@}
        ///@name Public operations
        ///@{

        static const Variable<double>& GetDerivativeVariable();

        void CalculateEffectiveVelocityDerivatives(
            BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateEffectiveKinematicViscosityDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateReactionTermDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateSourceTermDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        ///@}
    };

    class EpsilonDerivatives
        : public ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, AdjointElementData>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, AdjointElementData>;
        static constexpr unsigned int TDerivativesSize = BaseType::TDerivativesSize;

        ///@}
        ///@name Life cycle
        ///@{

        EpsilonDerivatives(const AdjointElementData& rElementData);

        ///@}
        ///@name Public operations
        ///@{

        static const Variable<double>& GetDerivativeVariable();

        void CalculateEffectiveVelocityDerivatives(
            BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateEffectiveKinematicViscosityDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateReactionTermDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateSourceTermDerivatives(
            BoundedVector<double, TDerivativesSize>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const override;

        ///@}
    };

    class ShapeSensitivity
        : public ConvectionDiffusionReactionSensitivityElementData<TDim, TNumNodes, SensitivityElementData>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = ConvectionDiffusionReactionSensitivityElementData<TDim, TNumNodes, SensitivityElementData>;
        static constexpr unsigned int TDerivativesSize = TDim * TNumNodes;

        ///@}
        ///@name Life cycle
        ///@{

        ShapeSensitivity(const SensitivityElementData& rElementData);

        ///@}
        ///@name Public operations
        ///@{

        static const Variable<array_1d<double, 3>>& GetDerivativeVariable();

        array_1d<double, 3> CalculateEffectiveVelocityDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

        double CalculateEffectiveKinematicViscosityDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

        double CalculateReactionTermDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

        double CalculateSourceTermDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

        ///@}
    };

    /**
     * @brief Data container for adjoint elements
     *
     * This class holds data for adjoint elements. Only one object
     * of this class is designed to be used with all the derivative classes
     * making CalculateGaussPointData, CalculateContants method to be called
     * only once per gauss point, and once per all derivatives
     *
     */
    class AdjointElementData : public KEpsilonElementData::EpsilonElementData<TDim>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = KEpsilonElementData::EpsilonElementData<TDim>;
        using GeometryType = typename BaseType::GeometryType;

        ///@}
        ///@name Life cycle
        ///@{

        AdjointElementData(const GeometryType& rGeometry);

        ///@}
        ///@name Public operations
        ///@{

        void CalculateGaussPointData(
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const int Step = 0) override;

        static void Check(
            const GeometryType& rGeometry,
            const ProcessInfo& rCurrentProcessInfo);

        static const Variable<double>& GetAdjointScalarVariable();

        ///@}

    private:
        ///@name Private members
        ///{

        double mProductionTerm;
        double mReactionTerm;
        double mEffectiveKinematicViscosity;
        BoundedMatrix<double, TNumNodes, TDim> mNodalVelocity;
        BoundedVector<double, TNumNodes> mGaussTurbulentKinematicViscositySensitivitiesK;
        BoundedVector<double, TNumNodes> mGaussTurbulentKinematicViscositySensitivitiesEpsilon;

        ///@}
        ///@name Private friend class definitions
        ///@{

        // following classes are made friends in order to access private members of this class
        // as well as the base class. (To avoid having lots of setters and getters)
        friend class VelocityDerivatives;
        friend class KDerivatives;
        friend class EpsilonDerivatives;

        ///@}
    };

    class SensitivityElementData : public KEpsilonElementData::EpsilonElementData<TDim>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = KEpsilonElementData::EpsilonElementData<TDim>;
        using GeometryType = typename BaseType::GeometryType;

        ///@}
        ///@name Life cycle
        ///@{

        SensitivityElementData(const GeometryType& rGeometry);

        ///@}
        ///@name Public operations
        ///@{

        void CalculateGaussPointData(
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const int Step = 0) override;

        static void Check(
            const GeometryType& rGeometry,
            const ProcessInfo& rCurrentProcessInfo);

        static const Variable<double>& GetAdjointScalarVariable();

        ///@}

    private:
        ///@name Private members
        ///{

        double mProductionTerm;
        double mReactionTerm;
        BoundedMatrix<double, TNumNodes, TDim> mNodalVelocity;

        ///@}
        ///@name Private friend class definitions
        ///@{

        // following classes are made friends in order to access private members of this class
        // as well as the base class. (To avoid having lots of setters and getters)
        friend class ShapeSensitivity;

        ///@}
    };

    ///@}
};
} // namespace KEpsilonAdjointElementData
} // namespace Kratos

#endif // KRATOS_K_EPSILON_HIGH_RE_DERIVATIVE_ELEMENT_DATA_EPSILON_ELEMENT_DATA_H_INCLUDED