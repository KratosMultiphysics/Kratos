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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <
    unsigned int TDim,
    unsigned int TNumNodes,
    class TElementDataType>
class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using PropertiesType = typename Element::PropertiesType;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

    using MatrixNN = BoundedMatrix<double, TNumNodes, TNumNodes>;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Element& rElement,
        const ProcessInfo& rProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    ///@}
    ///@name Class Declarations
    ///@{

    class Data;

    ///@}
    ///@name Classes
    ///@{

    class ResidualsContributions
    {
    public:
        ///@name Life Cycle
        ///@{

        ResidualsContributions(Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void AddGaussPointResidualsContributions(
            VectorN& rResidual,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        void AddResidualsContributionsAfterGaussPointLoop(
            VectorN& rResidual);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        MatrixNN mPrimalMatrix;

        ///@}
    };

    template<class TDerivativesType>
    class VariableDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        static constexpr double TSelfWeight = TDerivativesType::TSelfWeight;

        static constexpr IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        static constexpr IndexType TDerivativesSize = TDerivativeDimension * TNumNodes;

        ///@}
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorN& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative,
            const double MassTermsDerivativesWeight = 1.0);

        void CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(
            VectorN& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        ResidualsContributions mResidualWeightDerivativeContributions;
        BoundedVector<MatrixNN, TDerivativesSize> mPrimalDampingMatrixDerivatives;
        BoundedVector<double, TDerivativesSize> mScalarMultiplierDerivatives;

        ///@}
    };

    class SecondDerivatives
    {
    public:
        ///@name Life Cycle
        ///@{

        SecondDerivatives(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorN& rResidualDerivative,
            const int NodeIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        void CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(
            VectorN& rResidualDerivative,
            const int NodeIndex);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;

        BoundedVector<double, TNumNodes> mScalarMultiplierDerivatives;

        ///@}
    };

    class Data
    {
    public:
        ///@name Life Cycle
        ///@{

        Data(
            const Element& rElement,
            const ProcessInfo& rProcessInfo);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointData(
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const int Step = 0);

        void CalculateDataAfterGaussPointPointLoop();

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element& mrElement;
        TElementDataType mrElementData;
        double mNumberOfGaussPoints;

        // Primal data

        double mStabilizationTau;
        double mAbsoluteReactionTerm;
        double mVelocityMagnitude;
        double mElementLength;
        double mDeltaTime;
        double mBossakAlpha;
        double mBossakGamma;
        double mDynamicTau;
        double mPrimalRelaxedVariableRateValue;
        double mPrimalVariableValue;
        double mPrimalVariableGradientDotVelocity;
        double mScalarMultiplier;
        double mResidual;
        double mAbsoluteResidual;
        double mDiagonalCoefficient;

        ArrayD mPrimalVariableGradient;

        VectorN mVelocityConvectiveTerms;
        VectorN mPrimalVariableGradientDotDnDx;
        VectorN mNodalPrimalVariableValues;
        VectorN mNodalPrimalRelaxedVariableRateValues;

        MatrixNN mdNadNb;
        MatrixNN mPrimalDampingMatrix;
        MatrixNN mDiscreteDiffusionMatrix;

        double mDiscreteUpwindOperatorCoefficient;
        double mDiagonalPositivityPreservingCoefficient;

        ///@}
        ///@name Private Operations
        ///@{

        void AddDampingMatrixContributions(
            MatrixNN& rDampingMatrix,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX) const;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType>
        friend class VariableDerivatives;
        friend class SecondDerivatives;
        friend class ResidualsContributions;

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H
