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

#if !defined(KRATOS_K_OMEGA_SST_OMEGA_ELEMENT_DATA_DERIVATIVES_H)
#define KRATOS_K_OMEGA_SST_OMEGA_ELEMENT_DATA_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_stabilization_adjoint_utilities.h"
#include "custom_elements/data_containers/element_data_derivative.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"

namespace Kratos
{
namespace KOmegaSSTElementData
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class OmegaElementDataDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using UtilityDerivatives = typename ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>::Derivatives;

    ///@}
    ///@name Classes
    ///@{

    // forward declaration of the data gauss point data holder
    class Data;

    class UDerivative : public ElementDataDerivative<TDim>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = ElementDataDerivative<TDim>;

        using DataType = Data;

        using UtilitiesDerivative = typename UtilityDerivatives::NonRelatedVariable;

        static constexpr unsigned int TDerivativeDimension = TDim;

        static constexpr double TSelfWeight = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        UDerivative(
            const DataType& rData,
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        ArrayD CalculateEffectiveVelocityDerivative() const;

        double CalculateEffectiveKinematicViscosityDerivative() const;

        double CalculateReactionTermDerivative() const;

        double CalculateSourceTermDerivative() const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        double mGaussTurbulentKinematicViscosityDerivative;
        double mBlendedGammaDerivative;
        double mF1Derivative;
        double mCrossDiffusionDerivative;

        ///@}
    };

    class KDerivative : public ElementDataDerivative<TDim>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = ElementDataDerivative<TDim>;

        using DataType = Data;

        using UtilitiesDerivative = typename UtilityDerivatives::NonRelatedVariable;

        static constexpr unsigned int TDerivativeDimension = 1;

        static constexpr double TSelfWeight = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        KDerivative(
            const DataType& rData,
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        ArrayD CalculateEffectiveVelocityDerivative() const;

        double CalculateEffectiveKinematicViscosityDerivative() const;

        double CalculateReactionTermDerivative() const;

        double CalculateSourceTermDerivative() const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        double mGaussTurbulentKinematicViscosityDerivative;
        double mBlendedGammaDerivative;
        double mF1Derivative;
        double mCrossDiffusionDerivative;

        ///@}
    };

    class OmegaDerivative : public ElementDataDerivative<TDim>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = ElementDataDerivative<TDim>;

        using DataType = Data;

        using UtilitiesDerivative = typename UtilityDerivatives::NonRelatedVariable;

        static constexpr unsigned int TDerivativeDimension = 1;

        static constexpr double TSelfWeight = 1.0;

        ///@}
        ///@name Life Cycle
        ///@{

        OmegaDerivative(
            const DataType& rData,
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        ArrayD CalculateEffectiveVelocityDerivative() const;

        double CalculateEffectiveKinematicViscosityDerivative() const;

        double CalculateReactionTermDerivative() const;

        double CalculateSourceTermDerivative() const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        double mGaussTurbulentKinematicViscosityDerivative;
        double mBlendedGammaDerivative;
        double mF1Derivative;
        double mCrossDiffusionDerivative;

        ///@}
    };

    class ShapeDerivative : public ElementDataDerivative<TDim>
    {
    public:
        /// name@ Type Definitions
        ///@{

        using BaseType = ElementDataDerivative<TDim>;

        using DataType = Data;

        using UtilitiesDerivative = typename UtilityDerivatives::Shape;

        static constexpr unsigned int TDerivativeDimension = TDim;

        static constexpr double TSelfWeight = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        ShapeDerivative(
            const DataType& rData,
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        ///@}
        ///@name Operations
        ///@{

        const Variable<double>& GetDerivativeVariable() const;

        ArrayD CalculateEffectiveVelocityDerivative() const;

        double CalculateEffectiveKinematicViscosityDerivative() const;

        double CalculateReactionTermDerivative() const;

        double CalculateSourceTermDerivative() const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        double mGaussTurbulentKinematicViscosityDerivative;
        double mBlendedGammaDerivative;
        double mF1Derivative;
        double mCrossDiffusionDerivative;

        MatrixDD mVelocityGradientDerivative;
        ArrayD mTurbulentKineticEnergyGradientDerivative;
        ArrayD mTurbulentSpecificEnergyDissipationRateGradientDerivative;

        ///@}
    };

    class Data : public OmegaElementData<TDim>
    {
    public:
        ///@name Type Definitions
        ///@{

        using BaseType = OmegaElementData<TDim>;

        ///@}
        ///@name Life Cycle
        ///@{

        Data(
            const GeometryType& rGeometry,
            const Properties& rProperties,
            const ProcessInfo& rProcessInfo)
            : BaseType(rGeometry, rProperties, rProcessInfo)
        {}

        ///@}
        ///@name Static Operations
        ///@{

        using BaseType::GetScalarVariable;

        static const Variable<double>& GetAdjointScalarVariable();

        static void Check(
            const Element& rElement,
            const ProcessInfo& rCurrentProcessInfo);

        static const std::string GetName() { return "KOmegaSSTOmegaAdjointElementData"; }

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointData(
            const Vector& rN,
            const Matrix& rdNdX,
            const int Step = 0);

        ///@}

    protected:
        ///@name Protected Members
        ///@{

        using BaseType::mEffectiveVelocity;
        using BaseType::mEffectiveKinematicViscosity;
        using BaseType::mReactionTerm;
        using BaseType::mSourceTerm;
        using BaseType::mVelocityGradient;
        using BaseType::mTurbulentKineticEnergyGradient;
        using BaseType::mTurbulentSpecificEnergyDissipationRateGradient;
        using BaseType::mBeta1;
        using BaseType::mBeta2;
        using BaseType::mSigmaOmega1;
        using BaseType::mSigmaOmega2;
        using BaseType::mBetaStar;
        using BaseType::mTurbulentKineticEnergy;
        using BaseType::mTurbulentSpecificEnergyDissipationRate;
        using BaseType::mKinematicViscosity;
        using BaseType::mWallDistance;
        using BaseType::mCrossDiffusion;
        using BaseType::mF1;
        using BaseType::mBlendedSigmaOmega;
        using BaseType::mBlendedBeta;
        using BaseType::mBlendedGamma;
        using BaseType::mVelocityDivergence;
        using BaseType::mTurbulentKinematicViscosity;
        using BaseType::mKappa;
        using BaseType::mDensity;
        using BaseType::mA1;

        double mProductionTerm;
        double mF2;
        double mT;
        double mGamma1;
        double mGamma2;
        double mNonZeroOmega;

        MatrixDD mSymmetricVelocityGradient;
        MatrixND mNodalVelocity;


        ///@}
        ///@name Private Friends
        ///@{

        friend class UDerivative;
        friend class KDerivative;
        friend class OmegaDerivative;
        friend class ShapeDerivative;

        ///@}
    };

    ///@}
private:
    ///@name Private Static Operations
    ///@{

    ///@}
};

///@}

} // namespace KOmegaSSTElementData

} // namespace Kratos

#endif // KRATOS_K_OMEGA_SST_OMEGA_ELEMENT_DATA_DERIVATIVES_H