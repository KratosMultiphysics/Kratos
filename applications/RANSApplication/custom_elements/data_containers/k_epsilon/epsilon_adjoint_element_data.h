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
class EpsilonAdjointStateDerivatives
{
public:
    ///@name Public forward declarations
    ///@{

    class Data;

    ///@}
    ///@name Public classes
    ///@{

    class VelocityDerivatives
        : public ConvectionDiffusionReactionAdjointElementData<TDim, TNumNodes, Data>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType =
            ConvectionDiffusionReactionAdjointElementData<TDim, TNumNodes, Data>;
        static constexpr unsigned int TDerivativesSize = BaseType::TDerivativesSize;

        ///@}
        ///@name Life cycle
        ///@{

        VelocityDerivatives(const Data& rElementData);

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

        void CalculateReactionTermDerivatives(BoundedVector<double, TDerivativesSize>& rOutput,
                                              const Vector& rShapeFunctions,
                                              const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateSourceTermDerivatives(BoundedVector<double, TDerivativesSize>& rOutput,
                                            const Vector& rShapeFunctions,
                                            const Matrix& rShapeFunctionDerivatives) const override;

        ///@}
    };

    class KDerivatives
        : public ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, Data>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, Data>;
        static constexpr unsigned int TDerivativesSize = BaseType::TDerivativesSize;

        ///@}
        ///@name Life cycle
        ///@{

        KDerivatives(const Data& rElementData);

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

        void CalculateReactionTermDerivatives(BoundedVector<double, TDerivativesSize>& rOutput,
                                              const Vector& rShapeFunctions,
                                              const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateSourceTermDerivatives(BoundedVector<double, TDerivativesSize>& rOutput,
                                            const Vector& rShapeFunctions,
                                            const Matrix& rShapeFunctionDerivatives) const override;

        ///@}
    };

    class EpsilonDerivatives
        : public ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, Data>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = ConvectionDiffusionReactionAdjointElementData<1, TNumNodes, Data>;
        static constexpr unsigned int TDerivativesSize = BaseType::TDerivativesSize;

        ///@}
        ///@name Life cycle
        ///@{

        EpsilonDerivatives(const Data& rElementData);

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

        void CalculateReactionTermDerivatives(BoundedVector<double, TDerivativesSize>& rOutput,
                                              const Vector& rShapeFunctions,
                                              const Matrix& rShapeFunctionDerivatives) const override;

        void CalculateSourceTermDerivatives(BoundedVector<double, TDerivativesSize>& rOutput,
                                            const Vector& rShapeFunctions,
                                            const Matrix& rShapeFunctionDerivatives) const override;

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
    class Data : public KEpsilonElementData::EpsilonElementData<TDim>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = KEpsilonElementData::EpsilonElementData<TDim>;
        using GeometryType = typename BaseType::GeometryType;

        ///@}
        ///@name Life cycle
        ///@{

        Data(const GeometryType& rGeometry);

        ///@}
        ///@name Public operations
        ///@{

        void CalculateGaussPointData(const Vector& rShapeFunctions,
                                     const Matrix& rShapeFunctionDerivatives,
                                     const int Step = 0) override;

        static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

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

    ///@}
};

template <unsigned int TDim, unsigned int TNumNodes>
class EpsilonAdjointShapeDerivatives
{
public:
    ///@name Public forward declarations
    ///@{

    class Data;

    ///@}
    ///@name Public classes
    ///@{

    class ShapeDerivatives
        : public ConvectionDiffusionReactionSensitivityElementData<TDim, TNumNodes, Data>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType =
            ConvectionDiffusionReactionSensitivityElementData<TDim, TNumNodes, Data>;
        static constexpr unsigned int TDerivativesSize = TDim * TNumNodes;

        ///@}
        ///@name Life cycle
        ///@{

        ShapeDerivatives(const Data& rElementData);

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

    class Data : public KEpsilonElementData::EpsilonElementData<TDim>
    {
    public:
        ///@name Public type definitions
        ///@{

        using BaseType = KEpsilonElementData::EpsilonElementData<TDim>;
        using GeometryType = typename BaseType::GeometryType;

        ///@}
        ///@name Life cycle
        ///@{

        Data(const GeometryType& rGeometry);

        ///@}
        ///@name Public operations
        ///@{

        void CalculateGaussPointData(const Vector& rShapeFunctions,
                                     const Matrix& rShapeFunctionDerivatives,
                                     const int Step = 0) override;

        static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

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
        friend class ShapeDerivatives;

        ///@}
    };

    ///@}
};

///@}

} // namespace KEpsilonAdjointElementData
} // namespace Kratos

#endif // KRATOS_K_EPSILON_HIGH_RE_DERIVATIVE_ELEMENT_DATA_EPSILON_ELEMENT_DATA_H_INCLUDED