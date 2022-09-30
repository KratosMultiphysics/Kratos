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

#if !defined(KRATOS_K_OMEGA_SST_OMEGA_U_BASED_WALL_CONDITION_DATA_DERIVATIVES_H)
#define KRATOS_K_OMEGA_SST_OMEGA_U_BASED_WALL_CONDITION_DATA_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_conditions/data_containers/k_omega_sst/omega_u_based_wall_condition_data.h"
#include "custom_conditions/scalar_wall_flux_condition_data.h"

namespace Kratos
{
namespace KOmegaSSTWallConditionData
{
///@name Kratos Classes
///@{

template <unsigned int TDim>
class OmegaUBasedWallConditionDataDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    ///@}
    ///@name Classes
    ///@{

    // forward declaration of the data gauss point data holder
    class Data;

    class UDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = TDim;

        ///@}
        ///@name Life Cycle
        ///@{

        UDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };

    class KDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = 1;

        ///@}
        ///@name Life Cycle
        ///@{

        KDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };

    class OmegaDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = 1;

        ///@}
        ///@name Life Cycle
        ///@{

        OmegaDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };

    class ShapeDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = TDim;

        ///@}
        ///@name Life Cycle
        ///@{

        ShapeDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };



    class Data : public OmegaUBasedWallConditionData
    {
    public:
        ///@name Type Definitions
        ///@{

        using BaseType = OmegaUBasedWallConditionData;

        ///@}
        ///@name Life Cycle
        ///@{

        Data(
            const GeometryType& rGeometry,
            const Properties& rProperties,
            const ProcessInfo& rProcessInfo,
            const ScalarWallFluxConditionData::Parameters& rParameters);

        ///@}
        ///@name Static Operations
        ///@{

        using BaseType::GetScalarVariable;

        static const Variable<double>& GetAdjointScalarVariable();

        static void Check(
            const Condition& rCondition,
            const ProcessInfo& rCurrentProcessInfo);

        static const std::string GetName() { return "KOmegaSSTOmegaKBasedWallConditionAdjointData"; }

        static void InitializeCondition(
            Condition& rCondition,
            const ProcessInfo& rProcessInfo) {}

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointData(
            const Vector& rN,
            const int Step = 0);

        double GetWallFlux() const { return mWallFlux; }

        ///@}

    protected:
        ///@name Protected Members
        ///@{

        const GeometryType& mrParentElementGeometry;

        const ScalarWallFluxConditionData::Parameters& mrParameters;

        using BaseType::mBeta;
        using BaseType::mCmu25;
        using BaseType::mSigmaOmega1;
        using BaseType::mSigmaOmega2;
        using BaseType::mBetaStar;
        using BaseType::mWallHeight;
        using BaseType::mWallCorrectionFactor;

        double mUTau;
        double mTurbulentKineticEnergy;
        double mTurbulentSpecificEnergyDissipationRate;
        double mWallVelocityMagnitude;
        double mNormalMagnitude;
        double mYPlusLimit;
        double mWallFlux;
        double mF1;
        double mBlendedSigma;

        array_1d<double, 3> mUnitNormal;
        array_1d<double, 3> mWallVelocity;

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
};

///@}

} // namespace KOmegaSSTElementData

} // namespace Kratos

#endif // KRATOS_K_OMEGA_SST_OMEGA_U_BASED_WALL_CONDITION_DATA_DERIVATIVES_H