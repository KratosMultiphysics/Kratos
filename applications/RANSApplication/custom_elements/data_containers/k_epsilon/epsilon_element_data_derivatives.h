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

#if !defined(KRATOS_EPSILON_ELEMENT_DATA_DERIVATIVES_H)
#define KRATOS_EPSILON_ELEMENT_DATA_DERIVATIVES_H

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
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"

namespace Kratos
{
namespace KEpsilonElementData
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class EpsilonElementDataDerivatives
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
            const Matrix& rdNdXDerivative)
            : BaseType(
                NodeIndex,
                DirectionIndex,
                rData.GetGeometry(),
                W,
                rN,
                rdNdX,
                WDerivative,
                DetJDerivative,
                rdNdXDerivative),
              mrData(rData)
        {
        }

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
        double mGammaDerivative;

        ///@}
    };

    class EpsilonDerivative : public ElementDataDerivative<TDim>
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

        EpsilonDerivative(
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
        double mGammaDerivative;

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
            const Matrix& rdNdXDerivative)
            : BaseType(
                NodeIndex,
                DirectionIndex,
                rData.GetGeometry(),
                W,
                rN,
                rdNdX,
                WDerivative,
                DetJDerivative,
                rdNdXDerivative),
              mrData(rData)
        {
        }

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

        ///@}
    };

    class Data : public EpsilonElementData<TDim>
    {
    public:
        ///@name Type Definitions
        ///@{

        using BaseType = EpsilonElementData<TDim>;

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

        static const std::string GetName() { return "KEpsilonEpsilonAdjointElementData"; }

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
        using BaseType::mC1;
        using BaseType::mC2;
        using BaseType::mCmu;
        using BaseType::mGamma;
        using BaseType::mTurbulentKineticEnergy;
        using BaseType::mTurbulentKinematicViscosity;
        using BaseType::mKinematicViscosity;
        using BaseType::mVelocityDivergence;
        using BaseType::mInvEpsilonSigma;
        using BaseType::mDensity;

        double mTurbulentEnergyDissipationRate;
        double mProductionTerm;

        MatrixND mNodalVelocity;

        ///@}
        ///@name Private Friends
        ///@{

        friend class UDerivative;
        friend class KDerivative;
        friend class EpsilonDerivative;
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

} // namespace KEpsilonElementData

} // namespace Kratos

#endif // KRATOS_EPSILON_ELEMENT_DATA_DERIVATIVES_H